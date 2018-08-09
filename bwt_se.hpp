/******************************************************************************
 * bwt_se.hpp
 *
 * Parallel computation of the Burrows-Wheeler Transform
 *
 ******************************************************************************
 * Copyright (C) 2018 José Fuentes-Sepúlveda <jfuentess@dcc.uchile.cl>, Gonzalo
 * Navarro and Yakov Nekrich
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *****************************************************************************/

#ifndef BWT_SE_HPP
#define BWT_SE_HPP

#include "utils.hpp"
#include "lookup_table.hpp"
#include <parallel-range.h>

#include "external/sortPBBS/blockRadixSort.hpp"

/*
  Answer a batch of rank queries.

  Input:
  - B: Partial BWT
  - tbl: Lookup table
  - Q: Array of paris. Q[i].first represents the index of the rank query and
  Q[i].second represents the index of the metasymbol associated to the i-th
  query
  - m: Number of rank queries
  - size_B: Number of symbols in B
  - out: Answer of the rank queries
  

  The sequence B is (logically) divided into num_threads subsequences. To
  support rank queries, we will traverse the sequence B left-to-right using an
  universal table. We will take w consecutive symbols as the key the universal
  table, where w is the size of the machine word (16 bits in our case).

  Each rank query with index i will be answered in two steps. In the first step,
  the sequence B is traversed left-to-right until j=floor(i/w)*w, the closest
  multiple of w of i. In the second step, the remaining i-j bits are traversed,
  in order to complete the query.

  To support parallel processing, the first step is performed using domain
  decomposition, meanwhile the second step is performed naively, processing all
  the partial queries concurrently (each partial query takes O(1) by using
  universal tables)


 Q       :                 1 2           3 4    5              6    7 8     9
                           | |           | |    |              |    | |     |
 B       : ======|======|======|======|======|======|======|======|======|======

 Threads :|---------------------------|---------------------------|-------------|
 division          thread 1                     thread 2             thread 3
*/
template<typename T_in>
void batch_rank(std::vector<T_in> &B, lookup_table tbl,
		      std::pair<uint32_t, uint32_t> *Q, uint32_t m,
		      uint32_t size_B, uint32_t *out) {

  // Number of subsequences of tbl.npw consecutive symbols. It represents an
  // upper limit in the traversal of B
  uint32_t lm = Q[m-1].first/tbl.npw;
  uint32_t chunk = (uint32_t)ceil((double)lm/num_threads);
  uint32_t epch = chunk*tbl.npw; // Elements per chunk

  cilk_for(uint32_t i = 0; i < m; i++)
    out[i] = 0;
  
  uint32_t **partial = (uint32_t **)calloc(num_threads, sizeof(uint32_t *));
  for(uint32_t i = 0; i < num_threads; i++)
    partial[i] = (uint32_t *)calloc(tbl.sigma, sizeof(uint32_t));

  // Initial query of each chunk
  uint32_t *init_q = (uint32_t *)calloc(num_threads, sizeof(uint32_t));

  // 1. Detect the first query (leftmost) in each chunk
  cilk_for(uint32_t i=1; i < m; i++) {
    uint32_t chk_prev = Q[i-1].first/epch;
    uint32_t chk_curr = Q[i].first/epch;
    if(chk_curr < num_threads && chk_prev != chk_curr) {
      init_q[chk_curr] = i;
    }
  }

  // The first chunk containing at least one query
  uint32_t chk_first = Q[0].first/epch;  
  init_q[chk_first] = 0;

  // 2. Compute the rank of the last bits
  cilk_for(uint32_t k=0; k < m; k++) {
    uint32_t rem = Q[k].first%tbl.npw;
    uint32_t chk = Q[k].first-rem;
    lt_w w = 0, s = B[Q[k].first];
    if(rem) {
      // Packing symbols into a lt_w
      for(uint32_t i=chk, j=0; i < Q[k].first; i++, j++)
	w += ((lt_w)B[i] << tbl.bpe*j);
      out[k] += (uint32_t)tbl.popcount_i(w, s, rem);
    }
  }

  // 3. Compute partial ranks in each chunk. The final rank of each chunk is
  // stored in the array 'partial'
  cilk_for(uint32_t i=0; i < num_threads; i++) {
    uint32_t ll = i*epch;
    uint32_t ul = (i == num_threads-1)? std::min(Q[m-1].first+1, size_B): ll + epch;
    uint32_t curr_q = init_q[i], curr_i = Q[curr_q].first;
    uint32_t end_q = (i == num_threads-1)? m: init_q[i+1];

    for(uint j = ll; j < ul; j+=tbl.npw) {
      lt_w w = 0;

      while(curr_q < end_q && curr_i <= j+tbl.npw-1) {
      	lt_w s = B[Q[curr_q].first];
      	out[curr_q] += partial[i][s];
      	curr_q++;	

      	if(curr_q < end_q)
      	  curr_i = Q[curr_q].first;
      }
      
      for(uint32_t k=0; k < tbl.npw; k++)
      	w += ((lt_w)B[j+k] << tbl.bpe*k);
      uint8_t* pc = tbl.popcount(w);
      
      for(uint32_t j=0; j < tbl.sigma; j++)
      	  partial[i][j] += pc[j];      
    }
  }
  
  // 4. Prefix sum over the final ranks of each chunk
  for(uint32_t i=1; i < num_threads; i++)
    for(uint32_t j=0; j < tbl.sigma; j++)
      partial[i][j] += partial[i-1][j];
  
  // 5. Transform local ranks into global ranks
  cilk_for(uint32_t i=1; i < num_threads; i++) {
    uint32_t ll = init_q[i];
    uint32_t ul = (i == num_threads-1)? m: init_q[i+1];

    for(uint32_t j=ll; j < ul; j++) {
      lt_w s = B[Q[j].first];
      out[j] += partial[i-1][s];
    }
  }

  cilk_for(uint32_t i=0; i < m; i++) {
    Q[i].first = out[i];
  }

  for(uint32_t i = 0; i < num_threads; i++)
    free(partial[i]);
  free(partial);
  free(init_q);
}

/* Main function for the computation of the BWT

   Input:
   - input_file: Name of the file containing the input sequence
   - delta: Width (in bits) of the metasymbols
*/
template<typename T_in = uint8_t> // Data type of the input sequence
std::vector<T_in> bwt_se(const string& input_file, uint32_t delta) {

  T_in* map; // Table to convert symbols of the input sequence to its compact form
  uint32_t n; // Length of the BWT
  uint32_t sigma; // Alphabet size
  int32_t pad; // Number of $-symbols added to the input sequence
  
  // Load the input sequence
  int_vector<> infile;
  load_vector_from_file(infile, input_file, sizeof(T_in));

  /* [INIT] Measuring running time */
  struct timespec stime, etime; // To measure total running time
  clock_gettime(CLOCK_MONOTONIC , &stime);

  // Compute the compact form of the input sequence using 64-bits words
  uint64_t* seq = pack_symbols_par<T_in>(infile, delta, map, &n,
					 &sigma, &n_m, &pad);

  /***********************
   **** Steps 1 and 2 ****
   ***********************/

  // Compute the lookup table
  lookup_table tbl(sigma);

  bpe = ceil(log(sigma)/log(2)); // Bits per elements of input symbol
  epw = delta/bpe; // Elements per packed word

  // Define masks and threshold for parallel recursion
  mask_get_sym = ((uint64_t)1 << bpe)-1;
  mask_get_32 = ((uint64_t)1 << 32)-1;
  spawn_th = n_m/(5*delta);

  // Array to store the position of the metasymbols in the input sequence
  uint32_t *rotation = (uint32_t*)calloc(n_m+1, sizeof(uint32_t));
  // Representation of rotation of two symbols to the right of the input sequence
  cilk_for(int32_t i=0; i<n_m; i++)
    rotation[n_m-i-1] = i*epw+2;
  rotation[n_m] = n_m*epw;

  n_c = 2*n_m+1;

  // Convert the concatenation of two rotations of the input sequence to a
  // sequence with contiguous alphabet
  uint32_t *ISA = contig_alphabet(seq, rotation, 2, delta);

  // Computation of the suffix array of the ISA array
  uint32_t *SA = new uint32_t[n_c];

  cilk_for (uint32_t i = 0; i < n_c; ++i)
    SA[i] = i;
  parallelrangelite(ISA, SA, n_c);
  delete[] ISA;

  int32_t size_B = n_c-1; // Initial size of B. It will be increase iteratively  
  std::vector<T_in> B(n_m*epw); // Array to store the BWT

  int32_t tmp = 0, accum = 0;

  // Insert the first 2*n_m symbols into the BWT
  cilk_for(int32_t i=1; i < n_c; i++) {
    int32_t idx;
    if(SA[i] >= n_m) idx = rotation[SA[i]-n_m]+epw;
    else idx = rotation[SA[i]]+epw-1; // Rotation 1
    
    T_in s = get_sym<T_in>(seq, idx);
    B[i-1] = s;
  }

  // Array to store the position in B of the last inserted symbols
  pair<uint32_t, uint32_t> *W = new pair<uint32_t, uint32_t>[n_m];

  // Array to store the accumulated frequency of the symbols in B
  uint32_t * Acc = (uint32_t *)calloc(sigma+1, sizeof(uint32_t));

  cilk_for(int32_t i=1; i < n_c; i++) {
    if(SA[i] >= n_m) {
      W[SA[i] - n_m].first = i-1;
      W[SA[i] - n_m].second = SA[i] - n_m;      
    }
  }

  for(int32_t i=0; i < n_c-1; i++)
    Acc[B[i]]++;

  // In general, sigma is small enough to use a sequential algorithm to perform
  // the prefix sum
  for(int32_t i=0; i < sigma+1; i++) {
    tmp = Acc[i];
    Acc[i] = accum;
    accum += tmp;
  }

  delete[] SA;

  /************************
   **** Steps 3 to epw ****
   ************************/
  for(int32_t i=3; i <= epw; i++) {
    uint32_t * Acc_tmp = (uint32_t *)calloc(sigma+1, sizeof(uint32_t));

    // Compute a new rotation to the right
    cilk_for(int32_t i=0; i<n_m; i++)
      rotation[i]++;  

    // Restore the order in W. Notice that, at this point, the first and the
    // second components of the W entries are different, therefore, we could use
    // any sorting algorithm
    integerSort(W, n_m);
    
    // Convert the concatenation of two rotations of the input sequence to a
    // sequence with contiguous alphabet
    ISA = contig_alphabet(seq, rotation, i, delta);
    
    // Computation of the suffix array of the ISA array
    SA = new uint32_t[n_c];

    cilk_for (uint32_t i = 0; i < n_c; ++i)
      SA[i] = i;
    parallelrangelite(ISA, SA, n_c);

    // Answer the batch of rank queries. The queries are stored in W
    batch_rank<T_in>(B, tbl, W, n_m, size_B, ISA);

    // Add the accumulated frequency of the symbols to the answer of the rank
    // queries
    cilk_for(int32_t j=0; j < n_m; j++) {
      int32_t idx = rotation[W[j].second]+epw-1;
      T_in s = get_sym<T_in>(seq, idx);
      W[j].first += Acc[s];
    }

    uint chunk = (uint)ceil((double)n_c/num_threads);
    int32_t* partial_scan = new int32_t[num_threads];
    partial_scan[0] = 0;

    // Count suffixes of the first rotation to the right
    cilk_for(uint j = 0; j < num_threads-1; j++) {
      uint32_t ll = j*chunk;
      uint32_t ul = (ll+chunk > n_c)? n_c: ll+chunk;
      if(!j)
    	ll++;

      int32_t in_zero = 0;
      for(uint k = ll; k < ul; k++)
    	if(SA[k] < n_m)
    	  in_zero++;
      partial_scan[j+1] = in_zero;      
    }
    
    for(uint j = 1; j < num_threads; j++)
      partial_scan[j] += partial_scan[j-1];

    cilk_for(uint j = 0; j < num_threads; j++) {
      uint32_t ll = j*chunk;
      uint32_t ul = (ll+chunk > n_c)? n_c: ll+chunk;
      if(!j)
    	ll++;
      
      int32_t acc = partial_scan[j];
      int32_t in_zero = 0;
      for(uint k = ll; k < ul; k++)
    	if(SA[k] >= n_m)
    	  ISA[SA[k]-n_m] = (acc + in_zero);
    	else
    	  in_zero++;
    }

    // Compute the final position in B of the new symbols
    cilk_for(int32_t j=0; j < n_m; j++)
      W[j].first += ISA[W[j].second];

    // The arrays SA and ISA must be released here in order to reduce the memory
    // peak
    delete[] SA;
    delete[] ISA;

    // Sort W by the first component. At this point, W contains the position in
    // B of the new symbols
    integerSort(W, n_m);
    
    // Insert the new symbol in B
    // NOTE: We try to duplicate B temporarily in order to parallelize the
    // insertion of new symbols, but the extra space increases the total peak of
    // memory consumption. Additionally, the improvement in running time was
    // negligible. Therefore, how decided to use a sequential version
    int32_t cnt = size_B-1;
    size_B += n_m;
    int32_t ll2 = W[n_m-1].first+1, ul2 = size_B-1;

    for(int32_t k=ul2; k >= ll2; k--, cnt--)
      B[k] = B[cnt];
    
    for(int32_t j=n_m-1; j > 0 ; j--) {
      ul2 = W[j].first-j-1; ll2 = W[j-1].first-j+1; cnt = W[j].first-1;
      for(int32_t k=ul2; k >= ll2; k--, cnt--)
    	B[cnt] = B[k];
    }
    
    ul2 = W[0].first-1, ll2 = 0, cnt = W[0].first-1;

    for(int32_t k=ul2; k >= ll2; k--, cnt--)
      B[cnt] = B[k];

    integerSort2(W, n_m);

    for(int32_t j=0; j < n_m; j++) {
      B[W[j].first] = get_sym<T_in>(seq, rotation[j]+epw);
      Acc_tmp[B[W[j].first]]++;
    }

    // Update the array Acc
    tmp = 0, accum = 0;
    for(int32_t i=0; i < sigma+1; i++) {
      tmp = Acc_tmp[i];
      Acc_tmp[i] = accum;
      accum += tmp;
    }
    
    for(int32_t i=0; i < sigma+1; i++)
      Acc[i] += Acc_tmp[i];    

    free(Acc_tmp);
  }

  delete[] W;
  free(seq);
  free(Acc);

  // Prepare the final BWT deleting the added $-symbols
  prepare_bwt<T_in>(B, pad-1, map);

  /* [END] Measuring running time */
  clock_gettime(CLOCK_MONOTONIC , &etime); // Jose
  double et = (etime.tv_sec - stime.tv_sec) + (etime.tv_nsec - stime.tv_nsec) / 1000000000.0;
  printf(";%lf", et);

  return B;
}

#endif
