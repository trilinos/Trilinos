// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ZOLTAN_SORT_H
#define ZOLTAN_SORT_H

#include "zz_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Sorting */
void Zoltan_quicksort_pointer_dec_float_int (int*, float*, int*, int, int);
void Zoltan_quicksort_pointer_dec_float     (int*, float*, int,  int);
void Zoltan_quicksort_pointer_dec_double    (int*, double*, int,  int);

void Zoltan_quicksort_pointer_inc_float     (int*, float*, int,  int);
void Zoltan_quicksort_pointer_inc_double    (int*, double*, int,  int);
void Zoltan_quicksort_pointer_inc_int_int   (int*, int*,   int*, int, int);
void Zoltan_quicksort_pointer_inc_gno_int   (int*, ZOLTAN_GNO_TYPE *, int *, int, int);
void Zoltan_quicksort_pointer_inc_long_int  (int*, long*, int *, int, int);
void Zoltan_quicksort_list_inc_short        (short*, int*,   int,  int);
void Zoltan_quicksort_list_inc_one_int      (int*, int, int);
void Zoltan_quicksort_list_inc_int          (int*, int*,   int,  int);
void Zoltan_quicksort_list_inc_gno          (ZOLTAN_GNO_TYPE *, int*,   int,  int);
void Zoltan_quicksort_list_inc_long         (long*, int*,   int,  int);

void Zoltan_quicksort_pointer_inc_long_long_int   (int*, long long*, int *, int, int);
void Zoltan_quicksort_list_inc_long_long          (int64_t*, int*,   int,  int);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_SORT_H_ */
