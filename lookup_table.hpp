/******************************************************************************
 * lookup_table.hpp
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

#ifndef LOOKUPT_HPP
#define LOOKUPT_HPP

#include "utils.hpp"

typedef uint16_t lt_w; // lookup table word

struct lookup_table {
  uint8_t** T;
  lt_w mask; // mask used to delete unused bits. Notice that not all the bits of
             // lt_w are necessarily used to represent entries, since nub may be
             // different to sizeof(lt_w)*8 
  uint32_t sigma;
  uint32_t npw; // Number of elements per entry
  uint32_t bpe; // bits per element
  uint32_t nub; // Number of used bits
  lt_w width;
  
  // Constructor
  lookup_table(uint32_t _sigma): sigma(_sigma) {
    width = sizeof(lt_w)*8;
    bpe = bits::hi(sigma)+1; // bits per element
    npw = width/bpe; // Number of elements per entry
    nub = npw*bpe;
    uint32_t ne = pow(2,nub); // Number of entries
    mask = ((lt_w)1 << nub)-1;

    T = (uint8_t **)calloc(ne, sizeof(uint8_t *));
    for(uint32_t i=0; i < ne; i++)
      T[i] = (uint8_t *)calloc(sigma, sizeof(uint8_t));
    
    for(uint32_t i=0; i < ne; i++)
      for(uint32_t j=0; j < npw; j++) {
    	lt_w desp = bpe*j;
    	lt_w mask =((lt_w)1 << bpe)-1;
    	mask <<= desp;
    	uint8_t sym = (i & mask) >> desp;
    	if(sym < sigma)
    	  T[i][sym]++;
      }    
  }

  // Count the occurrence of each symbol in the word w
  uint8_t* popcount(lt_w w) {        
    return T[w &= mask];
  }

  // Count the occurrence of the symbol s in the word w up to the k-th symbol (exclusive)
  uint8_t popcount_i(lt_w w, uint8_t s, uint32_t k) {
    if(k > npw) {
      std::cerr << "Error: the index k must be less than " << npw;
      std::cerr << " (current is " << k << ")" << endl;
      exit(EXIT_FAILURE);
    }

    uint8_t out=0;

    w &= mask;

    for(uint32_t j=0; j < k; j++) {
      lt_w desp = bpe*j;
      lt_w mask_l =((lt_w)1 << bpe)-1; // local mask
      mask_l <<= desp;
      
      uint8_t ss = (w & mask_l);
      uint8_t sym = ((w & mask_l) >> desp);

      if(sym == s)
	out++;
    }
    
    return out;
  }
};

#endif // LOOKUPT_HPP
