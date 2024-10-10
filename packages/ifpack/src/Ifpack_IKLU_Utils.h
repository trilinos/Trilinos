/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_IKLU_UTILS_H
#define IFPACK_IKLU_UTILS_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

/* The code found in this file is adapted from CSparse Version 2.0.0
   written by Tim Davis, UFL
*/

/* --- primary CSparse routines and data structures ------------------------- */

typedef struct row_matrix    /* matrix in compressed-row or triplet form */
{
    int nzmax ;	    /* maximum number of entries */
    int m ;	    /* number of rows */
    int n ;	    /* number of columns */
    int *p ;	    /* row pointers (size m+1) or col indices (size nzmax) */
    int *j ;	    /* col indices, size nzmax */
    double *x ;	    /* numerical values, size nzmax */
    int nz ;	    /* # of entries in triplet matrix, -1 for compressed-row */
} csr ;

csr *csr_add (const csr *A, const csr *B, double alpha, double beta) ;
csr *csr_multiply (const csr *A, const csr *B) ;
double csr_norm (const csr *A) ;
int csr_print (const csr *A, int brief) ;
csr *csr_transpose (const csr *A, int values) ;

/* utilities */
void *csr_realloc (void *p, int n, size_t size, int *ok) ;

/* csr utilities */
csr *csr_spalloc (int m, int n, int nzmax, int values, int triplet) ;
csr *csr_spfree (csr *A) ;
int csr_sprealloc (csr *A, int nzmax) ;



/* --- secondary CSparse routines and data structures ----------------------- */
typedef struct cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    int *pinv ;	    /* inverse row perm. for QR, fill red. perm for Chol */
    int *q ;	    /* fill-reducing column permutation for LU and QR */
    int *parent ;   /* elimination tree for Cholesky and QR */
    int *cp ;	    /* column pointers for Cholesky, row counts for QR */
    int *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    int m2 ;	    /* # of rows for QR, after adding fictitious rows */
    double lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;    /* # entries in U for LU; in R for QR */
} css ;

typedef struct csr_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    csr *L ;	    /* L for LU and Cholesky, V for QR */
    csr *U ;	    /* U for LU, R for QR, not used for Cholesky */
    int *pinv ;	    /* partial pivoting for LU */
    int *perm ;	    /* partial pivoting for LU */
    double *B ;	    /* beta [0..n-1] for QR */
} csrn ;

typedef struct csr_dmperm_results    /* csr_dmperm or csr_scc output */
{
    int *p ;	    /* size m, row permutation */     /* may be back wards */
    int *q ;	    /* size n, column permutation */
    int *r ;	    /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    int *s ;	    /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    int nb ;	    /* # of blocks in fine dmperm decomposition */
    int rr [5] ;    /* coarse row decomposition */
    int cc [5] ;    /* coarse column decomposition */
} csrd ;

int *csr_amd (int order, const csr *A) ;

int csr_droptol (csr *A, double tol);
int csr_dropzeros (csr *A);
int csr_lsolve (const csr *L, double *x);
csrn *csr_lu (const csr *A, const css *S, double tol);
csr *csr_permute (const csr *A, const int *pinv, const int *q, int values);
css *csr_sqr (int order, const csr *A);
int csr_usolve (const csr *U, double *x);

/* utilities */
css *csr_sfree (css *S) ;
csrd *csr_dfree (csrd *D);
csrn *csr_nfree (csrn *N);



/* --- tertiary CSparse routines -------------------------------------------- */
double csr_cumsum (int *p, int *c, int n) ;
int csr_dfs (int j, csr *G, int top, int *xi, int *pstack, const int *pinv);
int csr_reach (csr *G, const csr *B, int k, int *xi, const int *pinv);
int csr_scatter (const csr *A, int j, double beta, int *w, double *x, int mark,
    csr *C, int nz) ;
csrd *csr_scc (csr *A);
int csr_spsolve (csr *G, const csr *B, int k, int *xi,
                 double *x, const int *pinv, int up);
int csr_tdfs (int j, int k, int *head, const int *next, int *post,
    int *stack) ;
/* utilities */
csrd *csr_dalloc (int m, int n);
/* csr utilities */
csrd *csr_ddone (csrd *D, csr *C, void *w, int ok) ;
csr *csr_done (csr *C, void *w, void *x, int ok) ;
int *csr_idone (int *p, csr *C, void *w, int ok) ;
csrn *csr_ndone (csrn *N, csr *C, void *w, void *x, int ok) ;

int csr_fkeep (csr *A, int (*fkeep) (int, int, double, void *), void *other);

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))

#endif // IFPACK_IKLU_UTILS_H
