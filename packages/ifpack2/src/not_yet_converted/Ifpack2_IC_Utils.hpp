// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_IC_UTILS_HPP
#define IFPACK2_IC_UTILS_HPP

typedef struct {
    double *val;  /* also known as A  */
    int    *col;  /* also known as JA; first column is column 0 */
    int    *ptr;  /* also known as IA; with ptr[0] = 0 */
} Ifpack2_AIJMatrix;

extern "C" {
void quicksort (int *const pbase, double *const daux, int total_elems);
}

void Ifpack2_AIJMatrix_dealloc(Ifpack2_AIJMatrix *a);

void crout_ict(
    int n,
#ifdef TIFPACK
    void * A,
    int maxentries,
    int (*getcol)( void * A, int col, int ** nentries, double * val, int * ind),
    int (*getdiag)( void *A, double * diag),
#else
    const Ifpack2_AIJMatrix *AL,
    const double *Adiag,
#endif
    double droptol,
    int lfil,
    Ifpack2_AIJMatrix *L,
    double **pdiag);


#endif
