/******************************************************************************
 * defs.hpp
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

#ifndef DEFS_HPP
#define DEFS_HPP

#include <stdint.h>

uint64_t mask_get_sym; // Mask to recover a symbol from a 64-bit word
uint64_t mask_get_32; // Mask to recover 32 bits from a 64-bit word
int32_t n_m; // Number of meta-symbols in the packed sequence
uint32_t n_c; // Number of meta-symbols in the concatenation of two packed sequences
uint32_t bpe; // Bits per element
uint32_t epw; // Elements per packed word

#define ISORT 25
uint32_t spawn_th; // Threshold for quicksort

#endif
