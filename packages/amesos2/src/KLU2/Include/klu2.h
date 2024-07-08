/* ========================================================================== */
/* === klu include file ===================================================== */
/* ========================================================================== */
// @HEADER
// *****************************************************************************
//                   KLU2: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the KLU2 contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER


/* Include file for user programs that call klu_* routines */

#ifndef _TKLU_H
#define _TKLU_H


#include "trilinos_amd.h"
#include "trilinos_colamd.h"
#include "trilinos_btf_decl.h"

/* -------------------------------------------------------------------------- */
/* Symbolic object - contains the pre-ordering computed by klu_analyze */
/* -------------------------------------------------------------------------- */

/* TODO : Entry is not needed in symbolic, numeric and common. Remove TODO */
template <typename Entry, typename Int> struct klu_symbolic
{
    /* A (P,Q) is in upper block triangular form.  The kth block goes from
     * row/col index R [k] to R [k+1]-1.  The estimated number of nonzeros
     * in the L factor of the kth block is Lnz [k]. 
     */

    /* only computed if the AMD ordering is chosen: */
    double symmetry ;   /* symmetry of largest block */
    double est_flops ;  /* est. factorization flop count */
    double lnz, unz ;   /* estimated nz in L and U, including diagonals */
    double *Lnz ;       /* size n, but only Lnz [0..nblocks-1] is used */

    /* computed for all orderings: */
    Int
        n,              /* input matrix A is n-by-n */
        nz,             /* # entries in input matrix */
        *P,             /* size n */
        *Q,             /* size n */
        *R,             /* size n+1, but only R [0..nblocks] is used */
        nzoff,          /* nz in off-diagonal blocks */
        nblocks,        /* number of blocks */
        maxblock,       /* size of largest block */
        ordering,       /* ordering used (AMD, COLAMD, or GIVEN) */
        do_btf ;        /* whether or not BTF preordering was requested */

    /* only computed if BTF preordering requested */
    Int structural_rank ;   /* 0 to n-1 if the matrix is structurally rank
                        * deficient.  -1 if not computed.  n if the matrix has
                        * full structural rank */

} ;

/* -------------------------------------------------------------------------- */
/* Numeric object - contains the factors computed by klu_factor */
/* -------------------------------------------------------------------------- */
template <typename Entry, typename Int> struct klu_numeric 
{
    /* LU factors of each block, the pivot row permutation, and the
     * entries in the off-diagonal blocks */

    Int n ;             /* A is n-by-n */
    Int nblocks ;       /* number of diagonal blocks */
    Int lnz ;           /* actual nz in L, including diagonal */
    Int unz ;           /* actual nz in U, including diagonal */
    Int max_lnz_block ; /* max actual nz in L in any one block, incl. diag */
    Int max_unz_block ; /* max actual nz in U in any one block, incl. diag */
    Int *Pnum ;         /* size n. final pivot permutation */
    Int *Pinv ;         /* size n. inverse of final pivot permutation */

    /* LU factors of each block */
    Int *Lip ;          /* size n. pointers into LUbx[block] for L */
    Int *Uip ;          /* size n. pointers into LUbx[block] for U */
    Int *Llen ;         /* size n. Llen [k] = # of entries in kth column of L */
    Int *Ulen ;         /* size n. Ulen [k] = # of entries in kth column of U */
    void **LUbx ;       /* L and U indices and entries (excl. diagonal of U) */
    size_t *LUsize ;    /* size of each LUbx [block], in sizeof (Unit) */
    void *Udiag ;       /* diagonal of U */

    /* scale factors; can be NULL if no scaling */
    double *Rs ;        /* size n. Rs [i] is scale factor for row i */

    /* permanent workspace for factorization and solve */
    size_t worksize ;   /* size (in bytes) of Work */
    void *Work ;        /* workspace */
    void *Xwork ;       /* alias into Numeric->Work */
    Int *Iwork ;        /* alias into Numeric->Work */

    /* off-diagonal entries in a conventional compressed-column sparse matrix */
    Int *Offp ;         /* size n+1, column pointers */
    Int *Offi ;         /* size nzoff, row indices */
    void *Offx ;        /* size nzoff, numerical values */
    Int nzoff ;

} ;

/* -------------------------------------------------------------------------- */
/* KLU control parameters and statistics */
/* -------------------------------------------------------------------------- */

/* Common->status values */
#define KLU_OK 0
#define KLU_SINGULAR (1)            /* status > 0 is a warning, not an error */
#define KLU_OUT_OF_MEMORY (-2)
#define KLU_INVALID (-3)
#define KLU_TOO_LARGE (-4)          /* integer overflow has occured */

template <typename Entry, typename Int>  struct klu_common
{

    /* --------------------------------------------------------------------- */
    /* parameters */
    /* --------------------------------------------------------------------- */

    double tol ;            /* pivot tolerance for diagonal preference */
    double memgrow ;        /* realloc memory growth size for LU factors */
    double initmem_amd ;    /* init. memory size with AMD: c*nnz(L) + n */
    double initmem ;        /* init. memory size: c*nnz(A) + n */
    double maxwork ;        /* maxwork for BTF, <= 0 if no limit */

    Int btf ;               /* use BTF pre-ordering, or not */
    Int ordering ;          /* 0: AMD, 1: COLAMD, 2: user P and Q,
                             * 3: user function */
    Int scale ;             /* row scaling: -1: none (and no error check),
                             * 0: none, 1: sum, 2: max */

    /* memory management routines */
    void *(*malloc_memory) (size_t) ;           /* pointer to malloc */
    void *(*realloc_memory) (void *, size_t) ;  /* pointer to realloc */
    void (*free_memory) (void *) ;              /* pointer to free */
    void *(*calloc_memory) (size_t, size_t) ;   /* pointer to calloc */

    /* pointer to user ordering function */
    int (*user_order) (int, int *, int *, int *, struct klu_common<Entry, Int> *) ;

    /* pointer to user data, passed unchanged as the last parameter to the
     * user ordering function (optional, the user function need not use this
     * information). */
    void *user_data ;

    Int halt_if_singular ;      /* how to handle a singular matrix:
        * FALSE: keep going.  Return a Numeric object with a zero U(k,k).  A
        *   divide-by-zero may occur when computing L(:,k).  The Numeric object
        *   can be passed to klu_solve (a divide-by-zero will occur).  It can
        *   also be safely passed to klu_refactor.
        * TRUE: stop quickly.  klu_factor will free the partially-constructed
        *   Numeric object.  klu_refactor will not free it, but will leave the
        *   numerical values only partially defined.  This is the default. */

    /* ---------------------------------------------------------------------- */
    /* statistics */
    /* ---------------------------------------------------------------------- */

    Int status ;                /* KLU_OK if OK, < 0 if error */
    Int nrealloc ;              /* # of reallocations of L and U */

    Int structural_rank ;       /* 0 to n-1 if the matrix is structurally rank
        * deficient (as determined by maxtrans).  -1 if not computed.  n if the
        * matrix has full structural rank.  This is computed by klu_analyze
        * if a BTF preordering is requested. */

    Int numerical_rank ;        /* First k for which a zero U(k,k) was found,
        * if the matrix was singular (in the range 0 to n-1).  n if the matrix
        * has full rank. This is not a true rank-estimation.  It just reports
        * where the first zero pivot was found.  -1 if not computed.
        * Computed by klu_factor and klu_refactor. */

    Int singular_col ;          /* n if the matrix is not singular.  If in the
        * range 0 to n-1, this is the column index of the original matrix A that
        * corresponds to the column of U that contains a zero diagonal entry.
        * -1 if not computed.  Computed by klu_factor and klu_refactor. */

    Int noffdiag ;      /* # of off-diagonal pivots, -1 if not computed */

    double flops ;      /* actual factorization flop count, from klu_flops */
    double rcond ;      /* crude reciprocal condition est., from klu_rcond */
    double condest ;    /* accurate condition est., from klu_condest */
    double rgrowth ;    /* reciprocal pivot rgrowth, from klu_rgrowth */
    double work ;       /* actual work done in BTF, in klu_analyze */

    size_t memusage ;   /* current memory usage, in bytes */
    size_t mempeak ;    /* peak memory usage, in bytes */

} ;

/* -------------------------------------------------------------------------- */
/* klu_defaults: sets default control parameters */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_defaults
(
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_analyze:  orders and analyzes a matrix */
/* -------------------------------------------------------------------------- */

/* Order the matrix with BTF (or not), then order each block with AMD, COLAMD,
 * a natural ordering, or with a user-provided ordering function */

template <typename Entry, typename Int>
klu_symbolic<Entry, Int> *klu_analyze
(
    /* inputs, not modified */
    Int n,              /* A is n-by-n */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_analyze_given: analyzes a matrix using given P and Q */
/* -------------------------------------------------------------------------- */

/* Order the matrix with BTF (or not), then use natural or given ordering
 * P and Q on the blocks.  P and Q are interpretted as identity
 * if NULL. */

template <typename Entry, typename Int>
klu_symbolic<Entry, Int> *klu_analyze_given
(
    /* inputs, not modified */
    Int n,              /* A is n-by-n */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    Int P [ ],          /* size n, user's row permutation (may be NULL) */
    Int Q [ ],          /* size n, user's column permutation (may be NULL) */
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_factor:  factors a matrix using the klu_analyze results */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
klu_numeric<Entry, Int> *klu_factor /* returns KLU_OK if OK, < 0 if error */
(
    /* inputs, not modified */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    Entry Ax [ ],      /* size nz, numerical values */
    klu_symbolic<Entry, Int> *Symbolic,
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_solve: solves Ax=b using the Symbolic and Numeric objects */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_solve
(
    /* inputs, not modified */
    klu_symbolic<Entry, Int> *Symbolic,
    klu_numeric<Entry, Int> *Numeric,
    Int ldim,               /* leading dimension of B */
    Int nrhs,               /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    Entry B [ ],           /* size ldim*nrhs */
    klu_common<Entry, Int> *Common
) ;


/* -------------------------------------------------------------------------- */
/* klu_tsolve: solves A'x=b using the Symbolic and Numeric objects */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_tsolve
(
    /* inputs, not modified */
    klu_symbolic<Entry, Int> *Symbolic,
    klu_numeric<Entry, Int> *Numeric,
    Int ldim,               /* leading dimension of B */
    Int nrhs,               /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    Entry B [ ],           /* size ldim*nrhs */
    Entry Xout [ ],           /* size ldim*1 */
    klu_common<Entry, Int> *Common
) ;


/* -------------------------------------------------------------------------- */
/* klu_solve2: solves Ax=b using the Symbolic and Numeric objects */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_solve2
(
    /* inputs, not modified */
    klu_symbolic<Entry, Int> *Symbolic,
    klu_numeric<Entry, Int> *Numeric,
    Int ldim,               /* leading dimension of B */
    Int nrhs,               /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    Entry B [ ],           /* size ldim*nrhs */
    Entry Xout [ ],           /* size ldim*1 */
    klu_common<Entry, Int> *Common
) ;


/* -------------------------------------------------------------------------- */
/* klu_tsolve2: solves A'x=b using the Symbolic and Numeric objects */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_tsolve2
(
    /* inputs, not modified */
    klu_symbolic<Entry, Int> *Symbolic,
    klu_numeric<Entry, Int> *Numeric,
    Int ldim,               /* leading dimension of B */
    Int nrhs,               /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    Entry B [ ],           /* size ldim*nrhs */
    klu_common<Entry, Int> *Common
) ;


/* -------------------------------------------------------------------------- */
/* klu_refactor: refactorizes matrix with same ordering as klu_factor */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_refactor         /* return TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    Entry Ax [ ],      /* size nz, numerical values */
    klu_symbolic<Entry, Int> *Symbolic,
    /* input, and numerical values modified on output */
    klu_numeric<Entry, Int> *Numeric,
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_free_symbolic: destroys the Symbolic object */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_free_symbolic
(
    klu_symbolic<Entry, Int> **Symbolic,
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_free_numeric: destroys the Numeric object */
/* -------------------------------------------------------------------------- */

/* Note that klu_free_numeric and klu_z_free_numeric are identical; each can
 * free both kinds of Numeric objects (real and complex) */

template <typename Entry, typename Int>
Int klu_free_numeric
(
    klu_numeric<Entry, Int> **Numeric,
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_sort: sorts the columns of the LU factorization */
/* -------------------------------------------------------------------------- */

/* this is not needed except for the MATLAB interface */

template <typename Entry, typename Int>
Int klu_sort
(
    /* inputs, not modified */
    klu_symbolic<Entry, Int> *Symbolic,
    /* input/output */
    klu_numeric<Entry, Int> *Numeric,
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_flops: determines # of flops performed in numeric factorzation */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_flops
(
    /* inputs, not modified */
    klu_symbolic<Entry, Int> *Symbolic,
    klu_numeric<Entry, Int> *Numeric,
    /* input/output */
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_rgrowth : compute the reciprocal pivot growth */
/* -------------------------------------------------------------------------- */

/* Pivot growth is computed after the input matrix is permuted, scaled, and
 * off-diagonal entries pruned.  This is because the LU factorization of each
 * block takes as input the scaled diagonal blocks of the BTF form.  The
 * reciprocal pivot growth in column j of an LU factorization of a matrix C
 * is the largest entry in C divided by the largest entry in U; then the overall
 * reciprocal pivot growth is the smallest such value for all columns j.  Note
 * that the off-diagonal entries are not scaled, since they do not take part in
 * the LU factorization of the diagonal blocks.
 *
 * In MATLAB notation:
 *
 * rgrowth = min (max (abs ((R \ A(p,q)) - F)) ./ max (abs (U))) */

template <typename Entry, typename Int>
Int klu_rgrowth
(
    Int Ap [ ],
    Int Ai [ ],
    Entry Ax [ ],
    klu_symbolic<Entry, Int> *Symbolic,
    klu_numeric<Entry, Int> *Numeric,
    klu_common<Entry, Int> *Common  /* Common->rgrowth = reciprocal pivot growth */
) ;

/* -------------------------------------------------------------------------- */
/* klu_condest */
/* -------------------------------------------------------------------------- */

/* Computes a reasonably accurate estimate of the 1-norm condition number, using
 * Hager's method, as modified by Higham and Tisseur (same method as used in
 * MATLAB's condest */

template <typename Entry, typename Int>
Int klu_condest
(
    Int Ap [ ],             /* size n+1, column pointers, not modified */
    Entry Ax [ ],           /* size nz = Ap[n], numerical values, not modified*/
    klu_symbolic<Entry, Int> *Symbolic, /* symbolic analysis, not modified */
    klu_numeric<Entry, Int> *Numeric,   /* numeric factorization, not modified */
    klu_common<Entry, Int> *Common      /* result returned in Common->condest */
) ;

/* -------------------------------------------------------------------------- */
/* klu_rcond: compute min(abs(diag(U))) / max(abs(diag(U))) */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_rcond
(
    klu_symbolic<Entry, Int> *Symbolic,         /* input, not modified */
    klu_numeric<Entry, Int> *Numeric,           /* input, not modified */
    klu_common<Entry, Int> *Common              /* result in Common->rcond */
) ;

/* -------------------------------------------------------------------------- */
/* klu_scale */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_scale           /* return TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    Int scale,          /* <0: none, no error check; 0: none, 1: sum, 2: max */
    Int n,
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    Entry Ax [ ],
    /* outputs, not defined on input */
    double Rs [ ],
    /* workspace, not defined on input or output */
    Int W [ ],          /* size n, can be NULL */
    klu_common<Entry, Int> *Common
) ;

/* -------------------------------------------------------------------------- */
/* klu_extract  */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
Int klu_extract     /* returns TRUE if successful, FALSE otherwise */
(
    /* inputs: */
    klu_numeric<Entry, Int> *Numeric,
    klu_symbolic<Entry, Int> *Symbolic,

    /* outputs, either allocated on input, or ignored otherwise */

    /* L */
    Int *Lp,        /* size n+1 */
    Int *Li,        /* size Numeric->lnz */
    Entry *Lx,     /* size Numeric->lnz */

    /* U */
    Int *Up,        /* size n+1 */
    Int *Ui,        /* size Numeric->unz */
    Entry *Ux,     /* size Numeric->unz */

    /* F */
    Int *Fp,        /* size n+1 */
    Int *Fi,        /* size Numeric->nzoff */
    Entry *Fx,     /* size Numeric->nzoff */

    /* P, row permutation */
    Int *P,         /* size n */

    /* Q, column permutation */
    Int *Q,         /* size n */

    /* Rs, scale factors */
    Entry *Rs,     /* size n */

    /* R, block boundaries */
    Int *R,         /* size Symbolic->nblocks+1 (nblocks is at most n) */

    klu_common<Entry, Int> *Common
) ;


/* -------------------------------------------------------------------------- */
/* KLU memory management routines */
/* -------------------------------------------------------------------------- */

template <typename Entry, typename Int>
void *klu_malloc        /* returns pointer to the newly malloc'd block */
(
    /* ---- input ---- */
    size_t n,           /* number of items */
    size_t size,        /* size of each item */
    /* --------------- */
    klu_common<Entry, Int> *Common
) ;

template <typename Entry, typename Int>
void *klu_free          /* always returns NULL */
(
    /* ---- in/out --- */
    void *p,            /* block of memory to free */
    size_t n,           /* number of items */
    size_t size,        /* size of each item */
    /* --------------- */
    klu_common<Entry, Int> *Common
) ;

template <typename Entry, typename Int>
void *klu_realloc       /* returns pointer to reallocated block */
(
    /* ---- input ---- */
    size_t nnew,        /* requested # of items in reallocated block */
    size_t nold,        /* current size of block, in # of items */
    size_t size,        /* size of each item */
    /* ---- in/out --- */
    void *p,            /* block of memory to realloc */
    /* --------------- */
    klu_common<Entry, Int> *Common
) ;

/* ========================================================================== */
/* === KLU version ========================================================== */
/* ========================================================================== */

/* All versions of KLU include these definitions.
 * As an example, to test if the version you are using is 1.2 or later:
 *
 *      if (KLU_VERSION >= KLU_VERSION_CODE (1,2)) ...
 *
 * This also works during compile-time:
 *
 *      #if (KLU >= KLU_VERSION_CODE (1,2))
 *          printf ("This is version 1.2 or later\n") ;
 *      #else
 *          printf ("This is an early version\n") ;
 *      #endif
 */

#define KLU_DATE "Mar 24, 2009"
#define KLU_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define KLU_MAIN_VERSION 1
#define KLU_SUB_VERSION 1
#define KLU_SUBSUB_VERSION 0
#define KLU_VERSION KLU_VERSION_CODE(KLU_MAIN_VERSION,KLU_SUB_VERSION)

#endif
