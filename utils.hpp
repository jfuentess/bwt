/******************************************************************************
 * utils.hpp
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

#ifndef UTILS_HPP
#define UTILS_HPP

#ifdef NOPARALLEL
#define __cilkrts_get_nworkers() 1
#define cilk_for for
#define cilk_spawn 
#define cilk_sync 
#else
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/common.h>
#endif

#define num_threads __cilkrts_get_nworkers()

#include <sdsl/int_vector.hpp>

#include "defs.hpp"

using namespace std;
using namespace sdsl;

/* Parallel exclusive prefix sum

   Input:
   - A: Input array
   - length: Number of elements
*/
template<typename T_in>
void scan(T_in *A, uint32_t length) {
  uint chunk = (uint)ceil((double)length/num_threads);
  T_in* partial_scan = new T_in[num_threads];
  cilk_for(uint i = 0; i < num_threads; i++) {
    uint32_t ll = i*chunk;
    uint32_t ul = ((ll + chunk) > length)? length: ll + chunk;

    T_in acc = 0, tmp = 0;
    for(uint j = ll; j < ul; j++) {
      tmp = A[j];
      A[j] = acc;
      acc += tmp;
    }
    partial_scan[i] = acc;
  }

  for(uint i = 1; i < num_threads; i++)
    partial_scan[i] += partial_scan[i-1];
    
  cilk_for(uint i = 1; i < num_threads; i++) {
    uint32_t ll = i*chunk;
    uint32_t ul = ((ll + chunk) > length)? length: ll + chunk;

    T_in acc = partial_scan[i-1];
    for(uint j = ll; j < ul; j++) {
      A[j] += acc;
    }
  }

  delete[] partial_scan;
}

/* Parallel exclusive prefix sum of an array of pairs. Only the first component
   of each entry is considered

   Input:
   - A: Input array
   - length: Number of elements
*/
template<typename T_in>
void scan(pair<T_in,T_in> *A, uint32_t length) {
  uint chunk = (uint)ceil((double)length/num_threads);
  T_in* partial_scan = new T_in[num_threads];
  cilk_for(uint i = 0; i < num_threads; i++) {
    uint32_t ll = i*chunk;
    uint32_t ul = ((ll + chunk) > length)? length: ll + chunk;

    T_in acc = 0, tmp = 0;
    for(uint j = ll; j < ul; j++) {
      tmp = A[j].first;
      A[j].first = acc;
      acc += tmp;
    }
    partial_scan[i] = acc;
  }

  for(uint i = 1; i < num_threads; i++)
    partial_scan[i] += partial_scan[i-1];
    
  cilk_for(uint i = 1; i < num_threads; i++) {
    uint32_t ll = i*chunk;
    uint32_t ul = ((ll + chunk) > length)? length: ll + chunk;

    T_in acc = partial_scan[i-1];
    for(uint j = ll; j < ul; j++) {
      A[j].first += acc;
    }
  }

  delete[] partial_scan;
}


/* Parallel computation of the maximum value in the vector A

   Input:
   - A: Input vector
*/
template<typename T_in>
T_in get_max(int_vector<>& A) {
  uint32_t length = A.size();
  
  uint chunk = (uint)ceil((double)length/num_threads);

  T_in* local_max = new T_in[num_threads];
  T_in global_max = 0;
  cilk_for(int i=0; i<num_threads; i++) {
    uint32_t ll = i*chunk;
    uint32_t ul = ((ll + chunk) > length)? length: ll + chunk;

    T_in out = 0;
    for(int j=ll; j<ul; j++) {
      if (A[j] > out)
	out = A[j];
    }
    local_max[i] = out;
  }

  for(int i=0; i<num_threads; i++) {
    if (local_max[i] > global_max)
      global_max = local_max[i];
  }

  return global_max;
}


/* Parallel computation of the maximum value in the array A

   Input:
   - A: Input array
   - length: Number of elements
*/
template<typename T_in>
uint32_t get_max(uint32_t *A, uint32_t length) {
  
  uint chunk = (uint)ceil((double)length/num_threads);

  T_in* local_max = new T_in[num_threads];
  T_in global_max = 0;
  cilk_for(int i=0; i<num_threads; i++) {
    uint32_t ll = i*chunk;
    uint32_t ul = ((ll + chunk) > length)? length: ll + chunk;

    T_in out = 0;
    for(int j=ll; j<ul; j++) {
      if (A[j] > out)
	out = A[j];
    }
    local_max[i] = out;
  }

  for(int i=0; i<num_threads; i++) {
    if (local_max[i] > global_max)
      global_max = local_max[i];
  }

  return global_max;
}

/* Parallel computation of the sum of all the elements in the array A

   Input:
   - A: Input array
   - length: Number of elements
*/
template<typename T_in>
T_in get_total(T_in *A, uint32_t length) {
  uint chunk = (uint)ceil((double)length/num_threads);

  T_in *local = (T_in *)calloc(num_threads, sizeof(T_in));
  T_in total = 0;
  
  cilk_for(uint i = 0; i < num_threads; i++) {
    uint32_t ll = i*chunk;
    uint32_t ul = ((ll + chunk) > length)? length: ll + chunk;

    for(uint j = ll; j < ul; j++) {
      if(A[j])
      	local[i] += A[j];
    }
  }

  for(int i=0; i<num_threads; i++) {
    total += local[i];
  }

  return total;
}

/* Compute the compact form of the sequence 'seq'
   
   Input
   - seq: Input sequence. Each symbol uses 8*T_in bits
   - delta: The width (in bits) of the metasymbols
   - table_map: Table to map symbols from the original input sequence to symbols
   in the compact form
   - n: Number of symbols
   - sigma: Alphabet size
   - n_m: Number of metasymbols
   - pad: Number of $-symbols add to the input sequence
*/
template<typename T_in>
uint64_t *pack_symbols_par(int_vector<>& seq, uint32_t delta, T_in * &table_map,
			   uint32_t *n, uint32_t *sigma, int32_t *n_m, int32_t
			   *pad) {
#if defined(MEM_DBG) || defined(MEM)
  malloc_count_reset_peak();
#endif  

  /* Computing a contiguous alphabet from the input sequence */
  *n = seq.size();
  
  T_in max_sym = get_max<T_in>(seq);
  uint32_t map_size = max_sym+1;

  T_in *local_map = (T_in*)calloc(map_size, sizeof(T_in));

  // Assumption: There is not 0-symbols in the input sequence
  local_map[0] = 1; // symbol $
  cilk_for (size_t i=0; i < *n; i++) {
    T_in sym = (T_in)seq[i];
    if(!local_map[sym])
      local_map[sym] = 1;
  }

  T_in local_sigma = get_total<T_in>(local_map, map_size);  
  *sigma = local_sigma;

  T_in flag = local_map[map_size-1];
  scan<T_in>(local_map, map_size);
   
  table_map = (T_in*)calloc(*sigma, sizeof(T_in));
  cilk_for(size_t i=0; i < map_size-1; i++)
    if((local_map[i+1]-local_map[i]))
      table_map[local_map[i]] = i;

  if(flag)
    table_map[local_map[map_size-1]] = map_size-1;

  /* Adding padding $ symbols */
  uint32_t wb = sizeof(T_in)*8;
  uint32_t bpe = ceil(log(*sigma)/log(2));
  uint32_t epw = delta/bpe; // Elements per packed word
  uint64_t new_n = ceil((double)(*n+1)/epw)*epw; // Adding padding ($)
                                                 // symbols to the input
                                                 // sequence

  new_n += epw; // We add a terminal metasymbol to avoid problems with periodicities

  *n +=1; // count the first $
  *pad = new_n-*n+1;
  *n_m = ceil((double)new_n/epw); // Number of metasymbols
  
  seq.resize(new_n);
  for(uint32_t i=*n; i < new_n; i++)
    seq[i] = 0; // $ symbol  

  // Reverse the input sequence
  uint32_t h = new_n/2;
  cilk_for(uint32_t i=0; i < h; i++) {
    int32_t j = new_n-i-1;
    T_in t = seq[j];
    seq[j] = seq[i];
    seq[i] = t;
  }

  // Copy the symbols of the first delta bits to avoid special cases for the
  // rotations
  new_n += epw;

  seq.resize(new_n+2*epw);
  for(uint32_t i=new_n-epw, j=0; i < new_n+epw; i++, j++)
    seq[i] = seq[j];

  double val = (double)new_n/64;
  int64_t n64 = ceil((double)val*bpe); // Number of 64-bits words to
				       // store the packed sequence
  uint64_t *pack_seq = (uint64_t *)calloc(n64, sizeof(uint64_t));

  /*
    The ith 64-bits word will contains completely the symbols in the range
    [(64i+bpe-1)/bpe , (64(i+1)-bpe)/bpe]
   */
  uint32_t op1 = bpe-1;
  uint32_t op2 = 64-bpe;
  cilk_for (int64_t i = 0; i < n64; i++) {
    int64_t v1 = 64*i+op1;
    int64_t v2 = 64*i+op2;
    int64_t l = (int64_t)v1/bpe;
    int64_t u = (int64_t)v2/bpe;

    if(u >= new_n)
      u = new_n-1;
    int32_t d = l*bpe- 64*i;

    // There exist a symbol splitted into two 64-bit words
    if(d && d < bpe) {
      uint64_t sym = (uint64_t)local_map[seq[l-1]];
      sym >>= (bpe-d);
      pack_seq[i] += sym;
    }
    
    for (int64_t j=l; j <= u; j++, d+=bpe) {
      uint64_t sym = (uint64_t)local_map[seq[j]];
      sym <<= d;
      pack_seq[i] += sym;
    }

    //    Again, there exist a symbol splitted into two 64-bit words
    if(d < 64) {
      uint64_t sym = (uint64_t)local_map[seq[u+1]];
      sym <<= d;
      pack_seq[i] += sym;
    }
  }

  util::clear(seq);

  return pack_seq;
}

/* Recover the i-th symbol of the sequence 'seq' */
template<typename T_in>
T_in get_sym(uint64_t *seq, uint32_t i) {
  uint32_t w2 = i*bpe/64;
  uint32_t d2 = i*bpe - w2*64;
  uint64_t idx = (uint64_t)i*bpe;
  uint32_t w = idx >> 6;
  uint32_t d = idx - (w << 6);
  uint64_t mask = mask_get_sym;

  mask <<= d;
  uint64_t sym = seq[w] & mask;
  sym >>= d;

  // The bits of the symbol belong to two consecutive 64-bit words
  d = 64-d;
  if(d < bpe) {
    uint64_t mask = ((uint64_t)1 << (bpe-d))-1;
    uint64_t sym2 = seq[w+1] & mask;
    sym2 <<= d;

    sym += sym2;
  }

  return (T_in)sym;
}

/* Recover 32 consecutive bits of 'seq' starting in the position 'i' */
uint32_t get_32b(uint64_t *seq, uint32_t i, uint32_t offset) {
  uint64_t idx = (uint64_t)i*bpe+offset;
  uint32_t w = idx >> 6;
  uint32_t d = idx - (w << 6);

  uint64_t mask = mask_get_32;
  mask <<= d;
  uint64_t sym = seq[w] & mask;
  sym >>= d;

  // The bits of the symbol belong to two consecutive 64-bit words
  d = 64-d;
  if(d < 32) {
    uint64_t mask = ((uint64_t)1 << (32-d))-1;
    uint64_t sym2 = seq[w+1] & mask;
    sym2 <<= d;

    sym += sym2;
  }
  return (uint32_t)sym;
}

/* Compare if a and b are equal or not */
uint8_t compare_equal(uint64_t *seq, int32_t a, int32_t b) {
  int32_t segments = (epw*bpe-1)/32;
  uint32_t d = epw*bpe-segments*32;
  uint32_t mask = ((uint64_t)1 << d)-1;

  for(int32_t i=0; i < segments; i++) {
    uint32_t v_a = get_32b(seq, a, 32*i);
    uint32_t v_b = get_32b(seq, b, 32*i);
    if(v_a != v_b)
      return 0;
  }

  uint32_t v_a = get_32b(seq, a, 32*segments);
  uint32_t v_b = get_32b(seq, b, 32*segments);
  v_a &= mask;
  v_b &= mask;

  if(v_a != v_b)
    return 0;
	
  return 1;
}

/* Compare if a is less than b */
uint8_t compare_less(uint64_t *seq, int32_t a, int32_t b) {
  int32_t segments = (epw*bpe-1)/32;
  uint32_t d = epw*bpe-segments*32;
  uint32_t mask = ((uint64_t)1 << d)-1;

  uint32_t v_a = get_32b(seq, a, 32*segments);
  uint32_t v_b = get_32b(seq, b, 32*segments);
  v_a &= mask;
  v_b &= mask;

  if(v_a < v_b) return 1;
  else if(v_a > v_b) return 0;
  
  for(int32_t i=segments-1; i >= 0; i--) {
    v_a = get_32b(seq, a, 32*i);
    v_b = get_32b(seq, b, 32*i);

    if(v_a < v_b) return 1;
    else if(v_a > v_b) return 0;
  }
	
  return 0;
}


/*******************
 *** Quick Sort  ***
 *******************/

/* Let T_c be the concatenation of two rotated sequences. Let i be the index of
   a metasymbol in T_c. This function returns the real index of such metasymbol
   in the input sequence
 */
template <class E>
E get_idx(E *rotation, E i, E n, E iter) {
  if((i >= n_m && i < n_c-1)) {
    return rotation[i-n_m];
  }
  else if(i < n_m) {
    return rotation[i]-iter+1;
  }
}

/* Function to compute the median of three values */
template <class S, class E>
E median(S* seq, E a, E b, E c) {
  if(compare_less(seq, a, b)) {
    if(compare_less(seq, b, c)) return b;
    else if(compare_less(seq, a, c)) return c;
    else return a;
  }
  else {
    if(compare_less(seq, a, c)) return a;
    else if(compare_less(seq, b, c)) return c;
    else return b;
  }
}

/* Sequential insertion sort

   The input has the same definition than the function split_par()
 */
template <class S, class E, class B>
void insertion_sort(S* seq, E* rotation, E *ref, E n, E iter, B *bits) {
  for (int32_t i=0; i < n; i++) {
    E r = ref[i];

    E v = get_idx(rotation, r, n, iter);
    E *A = ref + i - 1;
    E id_A = get_idx(rotation, *A, n, iter);

    while (A >= ref && compare_less(seq, v, id_A)) {
      *(A+1) = *A;
      A--;
      id_A = get_idx(rotation, *A, n, iter);
    }
    *(A+1) = r;
  }

  for (int32_t i=0; i < n-1; i++) {
    int32_t a = get_idx(rotation, ref[i], n, iter);
    int32_t b = get_idx(rotation, ref[i+1], n, iter);

    if(!compare_equal(seq, a, b)) {
      bits[i] = 1;
    }
  }
}

/* Sequential split function.

   The input has the same definition than the function split_par()
 */
template <class S, class E, class B>
std::pair<E*,E*> split_seq(S *seq, E *rotation, E *ref, E n, E iter, B *bits) {
  E id_a = get_idx(rotation, ref[n/4], n, iter);
  E id_b = get_idx(rotation, ref[n/2], n, iter);
  E id_c = get_idx(rotation, ref[(3*n)/4], n, iter);
  E p = median(seq, id_a, id_b, id_c);

  E *L = ref;   // below L are less than pivot
  E *M = ref;   // between L and M are equal to pivot
  E *R = ref+n-1; // above R are greater than pivot
  E id_M = 0, id_R = 0;
  while (1) {
    id_M = get_idx(rotation, *M, n, iter);
    while (!compare_less(seq, p, id_M)) {
      if (compare_less(seq, id_M, p)) std::swap(*M,*(L++));
      if (M >= R) break; 
      M++;

      id_M = get_idx(rotation, *M, n, iter);
    }

    id_R = get_idx(rotation, *R, n, iter);
    
    while (compare_less(seq, p, id_R)) {
      R--;
      id_R = get_idx(rotation, *R, n, iter);
    }
    
    if (M >= R) break; 
    std::swap(*M,*R--); 
    id_M = get_idx(rotation, *M, n, iter);
    if (compare_less(seq, id_M, p)) std::swap(*M,*(L++));
    M++;
  }
  if(M == R)
    M++;

  E b = L - ref;
  if(b)
    bits[b-1] = 1;
  b = R - ref;
  bits[b] = 1;

  return std::pair<E*,E*>(L,M);
}

/* Sequential split function. This function is used as part of the function
   split_par()

   Input
   - seq: Input sequence
   - rotation: Array with the starting position of the metasymbols in seq
   considering a rotation of 'iter' symbols to the right
   - ref: Position of the metasymbols in the rotated sequence. These position
   will be sorted according to the values of the metasymbols (it is the output)
   - n: Size of the sequence to be sorted (it is two times the number of
   metasymbols, since we are sorting the concatenation of two rotated sequences)
   - p: Pivot value
   - iter: The current iteration. Also, it can be interpreted as the value of the
   current rotation
   - l: In the output, all the values less than the pivot are below 'l'
   - m: In the output, all the values equal to the pivot are between 'l' and 'm'
   - r: In the output, all the values greater than the pivot are above 'r'
 */
template <class S, class E>
void _split(S *seq, E *rotation, E *ref, E n, E p, E iter, E *l, E *m, E *r) {
  E *L = ref;   // below L are less than pivot
  E *M = ref;   // between L and M are equal to pivot
  E *R = ref+n-1; // above R are greater than pivot
  E id_M = 0, id_R = 0;

  while (1) {
    id_M = get_idx(rotation, *M, n, iter);
    while (!compare_less(seq, p, id_M)) {
      if (compare_less(seq, id_M, p)) std::swap(*M,*(L++));
      if (M >= R) break; 
      M++;

      id_M = get_idx(rotation, *M, n, iter);
    }

    id_R = get_idx(rotation, *R, n, iter);
    
    while (compare_less(seq, p, id_R)) {
      R--;
      id_R = get_idx(rotation, *R, n, iter);
    }
    
    if (M >= R) break; 
    std::swap(*M,*R--); 
    id_M = get_idx(rotation, *M, n, iter);
    if (compare_less(seq, id_M, p)) std::swap(*M,*(L++));
    M++;
  }
  if(M == R)
    M++;

  *l = L - ref;
  *m = M - L;
  *r = ref + n - M;
}

/* Parallel split function
   
   The input has the same definition than the function quick_sort_par()
*/
template <class S, class E, class B>
std::pair<E*,E*> split_par(S *seq, E *rotation, E *ref, E n, E iter, B *bits, E *cpy) {
  E id_a = get_idx(rotation, ref[n/4], n, iter);
  E id_b = get_idx(rotation, ref[n/2], n, iter);
  E id_c = get_idx(rotation, ref[(3*n)/4], n, iter);
  E p = median(seq, id_a, id_b, id_c);

  int32_t chunk = n/num_threads;
  uint32_t *L = new uint32_t[num_threads]();
  uint32_t *M = new uint32_t[num_threads]();
  uint32_t *R = new uint32_t[num_threads]();

  cilk_for(int32_t i = 0; i < num_threads; i++) {
    int32_t os = i*chunk;
    int32_t nn = (i == num_threads-1)? n - i*chunk : chunk;
    if(nn > 0)
      _split(seq, rotation, ref+os, (uint32_t)nn, p, iter, &L[i], &M[i], &R[i]);
  }

  int32_t tot_l = 0, tot_m = 0, tot_r = 0;
  int32_t tmp = 0, acc = 0;
  for(int32_t i = 0; i < num_threads; i++) {
    tmp = L[i];
    L[i] = acc;
    acc += tmp;
  }
  tot_l = acc;
  for(int32_t i = 0; i < num_threads; i++) {
    tmp = M[i];
    M[i] = acc;
    acc += tmp;
  }
  tot_m = acc;
  for(int32_t i = 0; i < num_threads; i++) {
    tmp = R[i];
    R[i] = acc;
    acc += tmp;
  }
  tot_r = acc;

  cilk_for(int32_t i = 0; i < num_threads; i++) {
    // Lower
    int32_t os_c = L[i];
    int32_t nn = (i == num_threads-1)? tot_l-L[i] : L[i+1]-L[i];
    int32_t ll = i*chunk;
    int32_t ul = ll+nn;
    
    for(int32_t j=ll, k=0; j < ul; j++, k++)
      cpy[os_c+k] = ref[j];

    // Equal
    os_c = M[i];
    ll += nn;
    nn = (i == num_threads-1)? tot_m-M[i] : M[i+1]-M[i];
    ul = ll+nn;
    for(int32_t j=ll, k=0; j < ul; j++, k++)
      cpy[os_c+k] = ref[j];
    
    // Greater
    os_c = R[i];
    ll += nn;
    nn = (i == num_threads-1)? tot_r-R[i] : R[i+1]-R[i];
    ul = ll+nn;
    for(int32_t j=ll, k=0; j < ul; j++, k++)
      cpy[os_c+k] = ref[j];
  }

  cilk_for(int32_t i = 0; i < n; i++)
    ref[i] = cpy[i];

  E* LL = ref;
  E* MM = ref + tot_l;
  E* RR = ref + tot_m;

  E b = MM - ref;
  if(b)
    bits[b-1] = 1;
  b = RR - ref - 1;
  bits[b] = 1;
  
  delete[] L;
  delete[] M;
  delete[] R;
  return std::pair<E*,E*>(MM,RR);
}

/* Sequential quick sort 

   The input has the same definition than the function quick_sort_par()
*/
template <class S, class E, class B>
void quick_sort_seq(S* seq, E* rotation, E *ref, E n, E iter, B *bits) {
  if(!n) return;

  if (n < ISORT) insertion_sort(seq, rotation, ref, n, iter, bits);
  else {
    std::pair<E*,E*> X = split_seq(seq, rotation, ref, n, iter, bits);
    E nn = n - (ref+n-X.second);

    quick_sort_seq(seq, rotation, ref, E(X.first - ref), iter, bits);
    // Exclude all elts that equal pivot
    quick_sort_seq(seq, rotation, X.second, E(ref+n-X.second), iter, bits+nn); 
  }
}

/* Parallel quick sort based on the median of three elements as pivot and uses
   insertion sort for small inputs 
   
   Input
   - seq: Input sequence
   - rotation: Array with the starting position of the metasymbols in seq
   considering a rotation of 'iter' symbols to the right
   - ref: Position of the metasymbols in the rotated sequence. These position
   will be sorted according to the values of the metasymbols (it is the output)
   - n: Size of the sequence to be sorted (it is two times the number of
   metasymbols, since we are sorting the concatenation of two rotated sequences)
   - iter: The current iteration. Also, it can be interpreted as the value of the
   current rotation
   - spawn_th: Threshold in the number of parallel recursive calls
   - bits: Array of char to mark the ranges of equal symbols
   - cpy: Auxiliary array to copy elements of the ref array during the splitting
*/
template <class S, class E, class B>
void quick_sort_par(S* seq, E* rotation, E *ref, E n, E iter, E spawn_th, B *bits, E *cpy) {
  if (n < spawn_th) quick_sort_seq(seq, rotation, ref, n, iter, bits);
  else {
    std::pair<E*,E*> X = split_par(seq, rotation, ref, n, iter, bits, cpy);
    E nn = n - (ref+n-X.second);

    cilk_spawn quick_sort_par(seq, rotation, ref, E(X.first - ref), iter, spawn_th, bits, cpy);
    quick_sort_par(seq, rotation, X.second, E(ref+n-X.second), iter, spawn_th, bits+nn, cpy+nn);
    cilk_sync;
  }
}

/* Given the current rotation of the sequence 'seq', compute an equivalent sequence
   over a contiguous alphabet
   
   Input
   - seq: Input sequence
   - rotation: Array with the starting position of the metasymbols in seq
   considering a rotation of 'iter' symbols to the right
   - iter: The current iteration. Also, it can be interpreted as the value of the
   current rotation
   - delta: The width (in bits) of the metasymbols
*/
uint32_t* contig_alphabet(uint64_t *seq, uint32_t *rotation, uint32_t iter,
			 int32_t delta) {
  uint32_t *ref = new uint32_t[n_c-1];
  uint32_t* out = new uint32_t[n_c];
  uint8_t *bits = new uint8_t[n_c-1]();
  
  cilk_for(int i=0; i < n_c-1; i++)
    ref[i] = i;
  
  quick_sort_par(seq, rotation, ref, n_c-1, iter, spawn_th, bits, out);

  out[n_c-1] = 0;

  int32_t chunk = (int32_t)ceil((double)(n_c-1)/num_threads);
  int32_t *partial_scan = new int32_t[num_threads];

  for(uint32_t i = 0; i < num_threads; i++)
    partial_scan[i] = 0;
  
  cilk_for(uint32_t i = 0; i < num_threads; i++) {
    uint32_t ll = i*chunk;
    uint32_t ul = ((ll + chunk) > n_c-1)? n_c-1: ll + chunk;

    for(int32_t j = ll; j < ul; j++)
      partial_scan[i] += bits[j];
  }

  for(uint i = 1; i < num_threads; i++)
    partial_scan[i] += partial_scan[i-1];

  cilk_for(uint32_t i = 0; i < num_threads; i++) {
    uint32_t ll = i*chunk;
    uint32_t ul = ((ll + chunk) > n_c-1)? n_c-1: ll + chunk;

    int32_t acc = 0;
    int32_t o = (i > 0)? partial_scan[i-1]: 0;
    for(int32_t j = ll; j < ul; j++) {
      out[ref[j]] = acc + o;
      acc += bits[j];
    }
  }

  delete[] bits;
  delete[] ref;
  delete[] partial_scan;

  return out;
}

template<typename T_in>
void prepare_bwt(std::vector<T_in> &B, uint32_t r, T_in *map) {
  B[r] = B[0]; // Necessary to avoid problems with the rotations  
  int32_t n = B.size()-r;
  for(int i=0; i < n; i++)
    B[i] = map[B[i+r]];
}


#endif // UTILS_HPP
