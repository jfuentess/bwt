#ifndef A_RADIX_INCLUDED
#define A_RADIX_INCLUDED

#include <iostream>
#include <math.h>
#include "parallel.h"
#include "sequence.h"
#include "utils.h"
#include "transpose.h"

void integerSort(uintT *A, long n);
void integerSort(unsigned long *A, long n);
void integerSort(long *A, long n);

static void integerSort(uintT *A, long n, char* s);

void integerSort2(pair<ulong,uintT> *A, long n);

void integerSort(pair<ulong,uintT> *A, long n);

template <class T>
void integerSort(pair<uintT,T> *A, long n);

void integerSort(pair<uintT,uintT> *A, long n);

template <class T>
void integerSort(pair<uintT,T> *A, long n, char* s);

template <class T>
void integerSort2(pair<uintT,T> *A, long n);

void integerSort2(pair<uintT,uintT> *A, long n);


#endif

