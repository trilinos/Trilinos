/* ========================================================================== */
/* === klu include file ===================================================== */
/* ========================================================================== */

/* Include file for user programs that call klu_* routines */

#ifndef TRILINOS_KLU_DECL_H
#define TRILINOS_KLU_DECL_H

/* make it easy for C++ programs to include KLU */
#ifdef __cplusplus
extern "C" {
#endif

#include "trilinos_amd.h"
#include "trilinos_colamd.h"
#include "trilinos_btf_decl.h"

/* -------------------------------------------------------------------------- */
/* Symbolic object - contains the pre-ordering computed by klu_analyze */
/* -------------------------------------------------------------------------- */

typedef struct
{
    /* A (P,Q) is in upper block triangular form.  The kth block goes from
     * row/col index R [k] to R [k+1]-1.  The estimated number of nonzeros
     * in the L factor of the kth block is Lnz [k]. 
     */

    /* only computed if the AMD ordering is chosen: */
    double symmetry ;	/* symmetry of largest block */
    double est_flops ;	/* est. factorization flop count */
    double lnz, unz ;	/* estimated nz in L and U, including diagonals */
    double *Lnz ;	/* size n, but only Lnz [0..nblocks-1] is used */

    /* computed for all orderings: */
    int
	n,		/* input matrix A is n-by-n */
	nz,		/* # entries in input matrix */
	*P, 		/* size n */
	*Q,		/* size n */
	*R,		/* size n+1, but only R [0..nblocks] is used */
	nzoff,		/* nz in off-diagonal blocks */
	nblocks,	/* number of blocks */
	maxblock,	/* size of largest block */
	ordering,	/* ordering used (AMD, TRILINOS_COLAMD, or GIVEN) */
	do_btf ;	/* whether or not TRILINOS_BTF preordering was requested */

    /* only computed if TRILINOS_BTF preordering requested */
    int structural_rank ;   /* 0 to n-1 if the matrix is structurally rank
			* deficient.  -1 if not computed.  n if the matrix has
			* full structural rank */

} trilinos_klu_symbolic ;

typedef struct		/* 64-bit version (otherwise same as above) */
{
    double symmetry, est_flops, lnz, unz ;
    double *Lnz ;
    UF_long n, nz, *P, *Q, *R, nzoff, nblocks, maxblock, ordering, do_btf,
	structural_rank ;

} trilinos_klu_l_symbolic ;

/* -------------------------------------------------------------------------- */
/* Numeric object - contains the factors computed by klu_factor */
/* -------------------------------------------------------------------------- */

typedef struct
{
    /* LU factors of each block, the pivot row permutation, and the
     * entries in the off-diagonal blocks */

    int n ;		/* A is n-by-n */
    int nblocks ;	/* number of diagonal blocks */
    int lnz ;		/* actual nz in L, including diagonal */
    int unz ;		/* actual nz in U, including diagonal */
    int max_lnz_block ;	/* max actual nz in L in any one block, incl. diag */
    int max_unz_block ;	/* max actual nz in U in any one block, incl. diag */
    int *Pnum ;		/* size n. final pivot permutation */
    int *Pinv ;		/* size n. inverse of final pivot permutation */

    /* LU factors of each block */
    int *Lip ;		/* size n. pointers into LUbx[block] for L */
    int *Uip ;		/* size n. pointers into LUbx[block] for U */
    int *Llen ;		/* size n. Llen [k] = # of entries in kth column of L */
    int *Ulen ;		/* size n. Ulen [k] = # of entries in kth column of U */
    void **LUbx ;	/* L and U indices and entries (excl. diagonal of U) */
    size_t *LUsize ;	/* size of each LUbx [block], in sizeof (Unit) */
    void *Udiag ;	/* diagonal of U */

    /* scale factors; can be NULL if no scaling */
    double *Rs ;	/* size n. Rs [i] is scale factor for row i */

    /* permanent workspace for factorization and solve */
    size_t worksize ;	/* size (in bytes) of Work */
    void *Work ;	/* workspace */
    void *Xwork ;	/* alias into Numeric->Work */
    int *Iwork ;	/* alias into Numeric->Work */

    /* off-diagonal entries in a conventional compressed-column sparse matrix */
    int *Offp ;		/* size n+1, column pointers */
    int *Offi ;		/* size nzoff, row indices */
    void *Offx ;	/* size nzoff, numerical values */
    int nzoff ;

} trilinos_klu_numeric ;

typedef struct		/* 64-bit version (otherwise same as above) */
{
    UF_long n, nblocks, lnz, unz, max_lnz_block, max_unz_block, *Pnum, *Pinv,
	*Lip, *Uip, *Llen, *Ulen ;
    void **LUbx ;
    size_t *LUsize ;
    void *Udiag ;
    double *Rs ;
    size_t worksize ;
    void *Work, *Xwork ;
    UF_long *Iwork ;
    UF_long *Offp, *Offi ;
    void *Offx ;
    UF_long nzoff ;

} trilinos_klu_l_numeric ;

/* -------------------------------------------------------------------------- */
/* KLU control parameters and statistics */
/* -------------------------------------------------------------------------- */

/* Common->status values */
#define TRILINOS_KLU_OK 0
#define TRILINOS_KLU_SINGULAR (1)	    /* status > 0 is a warning, not an error */
#define TRILINOS_KLU_OUT_OF_MEMORY (-2)
#define TRILINOS_KLU_INVALID (-3)
#define TRILINOS_KLU_TOO_LARGE (-4)	    /* integer overflow has occured */

typedef struct trilinos_klu_common_struct
{

    /* ---------------------------------------------------------------------- */
    /* parameters */
    /* ---------------------------------------------------------------------- */

    double tol ;	    /* pivot tolerance for diagonal preference */
    double memgrow ;	    /* realloc memory growth size for LU factors */
    double initmem_amd ;    /* init. memory size with AMD: c*nnz(L) + n */
    double initmem ;	    /* init. memory size: c*nnz(A) + n */
    double maxwork ;	    /* maxwork for TRILINOS_BTF, <= 0 if no limit */

    int btf ;		    /* use TRILINOS_BTF pre-ordering, or not */
    int ordering ;	    /* 0: AMD, 1: TRILINOS_COLAMD, 2: user P and Q,
			     * 3: user function */
    int scale ;		    /* row scaling: -1: none (and no error check),
			     * 0: none, 1: sum, 2: max */

    /* memory management routines */
    void *(*malloc_memory) (size_t) ;		/* pointer to malloc */
    void *(*realloc_memory) (void *, size_t) ;  /* pointer to realloc */
    void (*free_memory) (void *) ;		/* pointer to free */
    void *(*calloc_memory) (size_t, size_t) ;	/* pointer to calloc */

    /* pointer to user ordering function */
    int (*user_order) (int, int *, int *, int *, struct trilinos_klu_common_struct *) ;

    /* pointer to user data, passed unchanged as the last parameter to the
     * user ordering function (optional, the user function need not use this
     * information). */
    void *user_data ;

    int halt_if_singular ;	/* how to handle a singular matrix:
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

    int status ;	        /* TRILINOS_KLU_OK if OK, < 0 if error */
    int nrealloc ;		/* # of reallocations of L and U */

    int structural_rank ;	/* 0 to n-1 if the matrix is structurally rank
	* deficient (as determined by maxtrans).  -1 if not computed.  n if the
	* matrix has full structural rank.  This is computed by klu_analyze
	* if a TRILINOS_BTF preordering is requested. */

    int numerical_rank ;	/* First k for which a zero U(k,k) was found,
	* if the matrix was singular (in the range 0 to n-1).  n if the matrix
	* has full rank. This is not a true rank-estimation.  It just reports
	* where the first zero pivot was found.  -1 if not computed.
	* Computed by klu_factor and klu_refactor. */

    int singular_col ;		/* n if the matrix is not singular.  If in the
	* range 0 to n-1, this is the column index of the original matrix A that
	* corresponds to the column of U that contains a zero diagonal entry.
	* -1 if not computed.  Computed by klu_factor and klu_refactor. */

    int noffdiag ;	/* # of off-diagonal pivots, -1 if not computed */

    double flops ;	/* actual factorization flop count, from klu_flops */
    double rcond ;	/* crude reciprocal condition est., from klu_rcond */
    double condest ;	/* accurate condition est., from klu_condest */
    double rgrowth ;	/* reciprocal pivot rgrowth, from klu_rgrowth */
    double work ;	/* actual work done in TRILINOS_BTF, in klu_analyze */

    size_t memusage ;	/* current memory usage, in bytes */
    size_t mempeak ;	/* peak memory usage, in bytes */

} trilinos_klu_common ;

typedef struct trilinos_klu_l_common_struct /* 64-bit version (otherwise same as above)*/
{

    double tol, memgrow, initmem_amd, initmem, maxwork ;
    UF_long btf, ordering, scale ;
    void *(*malloc_memory) (size_t) ;
    void *(*realloc_memory) (void *, size_t) ;
    void (*free_memory) (void *) ;
    void *(*calloc_memory) (size_t, size_t) ;
    UF_long (*user_order) (UF_long, UF_long *, UF_long *, UF_long *,
	struct trilinos_klu_l_common_struct *) ;
    void *user_data ;
    UF_long halt_if_singular ;
    UF_long status, nrealloc, structural_rank, numerical_rank, singular_col,
	noffdiag ;
    double flops, rcond, condest, rgrowth, work ;
    size_t memusage, mempeak ;

} trilinos_klu_l_common ;

/* -------------------------------------------------------------------------- */
/* klu_defaults: sets default control parameters */
/* -------------------------------------------------------------------------- */

int trilinos_klu_defaults
(
    trilinos_klu_common *Common
) ;

UF_long trilinos_klu_l_defaults (trilinos_klu_l_common *Common) ;

/* -------------------------------------------------------------------------- */
/* klu_analyze:  orders and analyzes a matrix */
/* -------------------------------------------------------------------------- */

/* Order the matrix with TRILINOS_BTF (or not), then order each block with AMD, TRILINOS_COLAMD,
 * a natural ordering, or with a user-provided ordering function */

trilinos_klu_symbolic *trilinos_klu_analyze
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    trilinos_klu_common *Common
) ;

trilinos_klu_l_symbolic *trilinos_klu_l_analyze (UF_long, UF_long *, UF_long *,
    trilinos_klu_l_common *Common) ;


/* -------------------------------------------------------------------------- */
/* klu_analyze_given: analyzes a matrix using given P and Q */
/* -------------------------------------------------------------------------- */

/* Order the matrix with TRILINOS_BTF (or not), then use natural or given ordering
 * P and Q on the blocks.  P and Q are interpretted as identity
 * if NULL. */

trilinos_klu_symbolic *trilinos_klu_analyze_given
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    int P [ ],		/* size n, user's row permutation (may be NULL) */
    int Q [ ],		/* size n, user's column permutation (may be NULL) */
    trilinos_klu_common *Common
) ;

trilinos_klu_l_symbolic *trilinos_klu_l_analyze_given (UF_long, UF_long *, UF_long *, UF_long *,
    UF_long *, trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_factor:  factors a matrix using the klu_analyze results */
/* -------------------------------------------------------------------------- */

trilinos_klu_numeric *trilinos_klu_factor	/* returns TRILINOS_KLU_OK if OK, < 0 if error */
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],	/* size nz, numerical values */
    trilinos_klu_symbolic *Symbolic,
    trilinos_klu_common *Common
) ;

trilinos_klu_numeric *trilinos_klu_z_factor      /* returns TRILINOS_KLU_OK if OK, < 0 if error */
(
     /* inputs, not modified */
     int Ap [ ],        /* size n+1, column pointers */
     int Ai [ ],        /* size nz, row indices */
     double Ax [ ],	/* size 2*nz, numerical values (real,imag pairs) */
     trilinos_klu_symbolic *Symbolic,
     trilinos_klu_common *Common
) ;

/* long / real version */
trilinos_klu_l_numeric *trilinos_klu_l_factor (UF_long *, UF_long *, double *, trilinos_klu_l_symbolic *,
    trilinos_klu_l_common *) ;

/* long / complex version */
trilinos_klu_l_numeric *trilinos_klu_zl_factor (UF_long *, UF_long *, double *, trilinos_klu_l_symbolic *,
    trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_solve: solves Ax=b using the Symbolic and Numeric objects */
/* -------------------------------------------------------------------------- */

int trilinos_klu_solve
(
    /* inputs, not modified */
    trilinos_klu_symbolic *Symbolic,
    trilinos_klu_numeric *Numeric,
    int ldim,		    /* leading dimension of B */
    int nrhs,		    /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ],	    /* size ldim*nrhs */
    trilinos_klu_common *Common
) ;

int trilinos_klu_z_solve
(
     /* inputs, not modified */
     trilinos_klu_symbolic *Symbolic,
     trilinos_klu_numeric *Numeric,
     int ldim,               /* leading dimension of B */
     int nrhs,               /* number of right-hand-sides */

     /* right-hand-side on input, overwritten with solution to Ax=b on output */
     double B [ ],	    /* size 2*ldim*nrhs */
     trilinos_klu_common *Common
) ;

UF_long trilinos_klu_l_solve (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, UF_long, UF_long,
    double *, trilinos_klu_l_common *) ;

UF_long trilinos_klu_zl_solve (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, UF_long, UF_long,
    double *, trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_tsolve: solves A'x=b using the Symbolic and Numeric objects */
/* -------------------------------------------------------------------------- */

int trilinos_klu_tsolve
(
    /* inputs, not modified */
    trilinos_klu_symbolic *Symbolic,
    trilinos_klu_numeric *Numeric,
    int ldim,		    /* leading dimension of B */
    int nrhs,		    /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ],	    /* size ldim*nrhs */
    trilinos_klu_common *Common
) ;

int trilinos_klu_z_tsolve
(
    /* inputs, not modified */
    trilinos_klu_symbolic *Symbolic,
    trilinos_klu_numeric *Numeric,
    int ldim,		    /* leading dimension of B */
    int nrhs,		    /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ],	    /* size 2*ldim*nrhs */
    int conj_solve,	    /* TRUE: conjugate solve, FALSE: solve A.'x=b */
    trilinos_klu_common *Common
     
) ;

UF_long trilinos_klu_l_tsolve (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, UF_long, UF_long,
    double *, trilinos_klu_l_common *) ;

UF_long trilinos_klu_zl_tsolve (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, UF_long, UF_long,
    double *, UF_long, trilinos_klu_l_common * ) ;


/* -------------------------------------------------------------------------- */
/* klu_refactor: refactorizes matrix with same ordering as klu_factor */
/* -------------------------------------------------------------------------- */

int trilinos_klu_refactor	    /* return TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],	/* size nz, numerical values */
    trilinos_klu_symbolic *Symbolic,
    /* input, and numerical values modified on output */
    trilinos_klu_numeric *Numeric,
    trilinos_klu_common *Common
) ;

int trilinos_klu_z_refactor	    /* return TRUE if successful, FALSE otherwise */
(
     /* inputs, not modified */
     int Ap [ ],	/* size n+1, column pointers */
     int Ai [ ],	/* size nz, row indices */
     double Ax [ ],	/* size 2*nz, numerical values */
     trilinos_klu_symbolic *Symbolic,
     /* input, and numerical values modified on output */
     trilinos_klu_numeric *Numeric,
     trilinos_klu_common *Common
) ;

UF_long trilinos_klu_l_refactor (UF_long *, UF_long *, double *, trilinos_klu_l_symbolic *,
    trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;

UF_long trilinos_klu_zl_refactor (UF_long *, UF_long *, double *, trilinos_klu_l_symbolic *,
    trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_free_symbolic: destroys the Symbolic object */
/* -------------------------------------------------------------------------- */

int trilinos_klu_free_symbolic
(
    trilinos_klu_symbolic **Symbolic,
    trilinos_klu_common *Common
) ;

UF_long trilinos_klu_l_free_symbolic (trilinos_klu_l_symbolic **, trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_free_numeric: destroys the Numeric object */
/* -------------------------------------------------------------------------- */

/* Note that klu_free_numeric and klu_z_free_numeric are identical; each can
 * free both kinds of Numeric objects (real and complex) */

int trilinos_klu_free_numeric
(
    trilinos_klu_numeric **Numeric,
    trilinos_klu_common *Common
) ;

int trilinos_klu_z_free_numeric
(
     trilinos_klu_numeric **Numeric,
     trilinos_klu_common *Common
) ;

UF_long trilinos_klu_l_free_numeric (trilinos_klu_l_numeric **, trilinos_klu_l_common *) ;
UF_long trilinos_klu_zl_free_numeric (trilinos_klu_l_numeric **, trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_sort: sorts the columns of the LU factorization */
/* -------------------------------------------------------------------------- */

/* this is not needed except for the MATLAB interface */

int trilinos_klu_sort
(
    /* inputs, not modified */
    trilinos_klu_symbolic *Symbolic,
    /* input/output */
    trilinos_klu_numeric *Numeric,
    trilinos_klu_common *Common
) ;

int trilinos_klu_z_sort
(
    /* inputs, not modified */
    trilinos_klu_symbolic *Symbolic,
    /* input/output */
    trilinos_klu_numeric *Numeric,
    trilinos_klu_common *Common
) ;

UF_long trilinos_klu_l_sort (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;
UF_long trilinos_klu_zl_sort (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_flops: determines # of flops performed in numeric factorzation */
/* -------------------------------------------------------------------------- */

int trilinos_klu_flops
(
    /* inputs, not modified */
    trilinos_klu_symbolic *Symbolic,
    trilinos_klu_numeric *Numeric,
    /* input/output */
    trilinos_klu_common *Common
) ;

int trilinos_klu_z_flops
(
    /* inputs, not modified */
    trilinos_klu_symbolic *Symbolic,
    trilinos_klu_numeric *Numeric,
    /* input/output */
    trilinos_klu_common *Common
) ;

UF_long trilinos_klu_l_flops (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;
UF_long trilinos_klu_zl_flops (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;



/* -------------------------------------------------------------------------- */
/* klu_rgrowth : compute the reciprocal pivot growth */
/* -------------------------------------------------------------------------- */

/* Pivot growth is computed after the input matrix is permuted, scaled, and
 * off-diagonal entries pruned.  This is because the LU factorization of each
 * block takes as input the scaled diagonal blocks of the TRILINOS_BTF form.  The
 * reciprocal pivot growth in column j of an LU factorization of a matrix C
 * is the largest entry in C divided by the largest entry in U; then the overall
 * reciprocal pivot growth is the smallest such value for all columns j.  Note
 * that the off-diagonal entries are not scaled, since they do not take part in
 * the LU factorization of the diagonal blocks.
 *
 * In MATLAB notation:
 *
 * rgrowth = min (max (abs ((R \ A(p,q)) - F)) ./ max (abs (U))) */

int trilinos_klu_rgrowth
(
    int Ap [ ],
    int Ai [ ],
    double Ax [ ],
    trilinos_klu_symbolic *Symbolic,
    trilinos_klu_numeric *Numeric,
    trilinos_klu_common *Common		/* Common->rgrowth = reciprocal pivot growth */
) ;

int trilinos_klu_z_rgrowth
(
    int Ap [ ],
    int Ai [ ],
    double Ax [ ],
    trilinos_klu_symbolic *Symbolic,
    trilinos_klu_numeric *Numeric,
    trilinos_klu_common *Common		/* Common->rgrowth = reciprocal pivot growth */
) ;

UF_long trilinos_klu_l_rgrowth (UF_long *, UF_long *, double *, trilinos_klu_l_symbolic *,
    trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;

UF_long trilinos_klu_zl_rgrowth (UF_long *, UF_long *, double *, trilinos_klu_l_symbolic *,
    trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_condest */
/* -------------------------------------------------------------------------- */

/* Computes a reasonably accurate estimate of the 1-norm condition number, using
 * Hager's method, as modified by Higham and Tisseur (same method as used in
 * MATLAB's condest */

int trilinos_klu_condest
(
    int Ap [ ],		    /* size n+1, column pointers, not modified */
    double Ax [ ],	    /* size nz = Ap[n], numerical values, not modified*/
    trilinos_klu_symbolic *Symbolic, /* symbolic analysis, not modified */
    trilinos_klu_numeric *Numeric,   /* numeric factorization, not modified */
    trilinos_klu_common *Common	    /* result returned in Common->condest */
) ;

int trilinos_klu_z_condest
(
    int Ap [ ],
    double Ax [ ],	    /* size 2*nz */
    trilinos_klu_symbolic *Symbolic,
    trilinos_klu_numeric *Numeric,
    trilinos_klu_common *Common	    /* result returned in Common->condest */
) ;

UF_long trilinos_klu_l_condest (UF_long *, double *, trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *,
    trilinos_klu_l_common *) ;

UF_long trilinos_klu_zl_condest (UF_long *, double *, trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *,
    trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_rcond: compute min(abs(diag(U))) / max(abs(diag(U))) */
/* -------------------------------------------------------------------------- */

int trilinos_klu_rcond
(
    trilinos_klu_symbolic *Symbolic,	    /* input, not modified */
    trilinos_klu_numeric *Numeric,	    /* input, not modified */
    trilinos_klu_common *Common		    /* result in Common->rcond */
) ;

int trilinos_klu_z_rcond
(
    trilinos_klu_symbolic *Symbolic,	    /* input, not modified */
    trilinos_klu_numeric *Numeric,	    /* input, not modified */
    trilinos_klu_common *Common		    /* result in Common->rcond */
) ;

UF_long trilinos_klu_l_rcond (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;

UF_long trilinos_klu_zl_rcond (trilinos_klu_l_symbolic *, trilinos_klu_l_numeric *, trilinos_klu_l_common *) ;



/* -------------------------------------------------------------------------- */
/* klu_scale */
/* -------------------------------------------------------------------------- */

int trilinos_klu_scale		/* return TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    int scale,		/* <0: none, no error check; 0: none, 1: sum, 2: max */
    int n,
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    /* outputs, not defined on input */
    double Rs [ ],
    /* workspace, not defined on input or output */
    int W [ ],		/* size n, can be NULL */
    trilinos_klu_common *Common
) ;

int trilinos_klu_z_scale		/* return TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    int scale,		/* <0: none, no error check; 0: none, 1: sum, 2: max */
    int n,
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    /* outputs, not defined on input */
    double Rs [ ],
    /* workspace, not defined on input or output */
    int W [ ],		/* size n, can be NULL */
    trilinos_klu_common *Common
) ;

UF_long trilinos_klu_l_scale (UF_long, UF_long, UF_long *, UF_long *, double *,
    double *, UF_long *, trilinos_klu_l_common *) ;

UF_long trilinos_klu_zl_scale (UF_long, UF_long, UF_long *, UF_long *, double *,
    double *, UF_long *, trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* klu_extract  */
/* -------------------------------------------------------------------------- */

int trilinos_klu_extract	    /* returns TRUE if successful, FALSE otherwise */
(
    /* inputs: */
    trilinos_klu_numeric *Numeric,
    trilinos_klu_symbolic *Symbolic,

    /* outputs, either allocated on input, or ignored otherwise */

    /* L */
    int *Lp,	    /* size n+1 */
    int *Li,	    /* size Numeric->lnz */
    double *Lx,	    /* size Numeric->lnz */

    /* U */
    int *Up,	    /* size n+1 */
    int *Ui,	    /* size Numeric->unz */
    double *Ux,	    /* size Numeric->unz */

    /* F */
    int *Fp,	    /* size n+1 */
    int *Fi,	    /* size Numeric->nzoff */
    double *Fx,	    /* size Numeric->nzoff */

    /* P, row permutation */
    int *P,	    /* size n */

    /* Q, column permutation */
    int *Q,	    /* size n */

    /* Rs, scale factors */
    double *Rs,	    /* size n */

    /* R, block boundaries */
    int *R,	    /* size Symbolic->nblocks+1 (nblocks is at most n) */

    trilinos_klu_common *Common
) ;


int trilinos_klu_z_extract	    /* returns TRUE if successful, FALSE otherwise */
(
    /* inputs: */
    trilinos_klu_numeric *Numeric,
    trilinos_klu_symbolic *Symbolic,

    /* outputs, all of which must be allocated on input */

    /* L */
    int *Lp,	    /* size n+1 */
    int *Li,	    /* size nnz(L) */
    double *Lx,	    /* size nnz(L) */
    double *Lz,	    /* size nnz(L) for the complex case, ignored if real */

    /* U */
    int *Up,	    /* size n+1 */
    int *Ui,	    /* size nnz(U) */
    double *Ux,	    /* size nnz(U) */
    double *Uz,	    /* size nnz(U) for the complex case, ignored if real */

    /* F */
    int *Fp,	    /* size n+1 */
    int *Fi,	    /* size nnz(F) */
    double *Fx,	    /* size nnz(F) */
    double *Fz,	    /* size nnz(F) for the complex case, ignored if real */

    /* P, row permutation */
    int *P,	    /* size n */

    /* Q, column permutation */
    int *Q,	    /* size n */

    /* Rs, scale factors */
    double *Rs,	    /* size n */

    /* R, block boundaries */
    int *R,	    /* size Symbolic->nblocks+1 (nblocks is at most n) */

    trilinos_klu_common *Common
) ;

UF_long trilinos_klu_l_extract (trilinos_klu_l_numeric *, trilinos_klu_l_symbolic *,
    UF_long *, UF_long *, double *,
    UF_long *, UF_long *, double *,
    UF_long *, UF_long *, double *,
    UF_long *, UF_long *, double *, UF_long *, trilinos_klu_l_common *) ;

UF_long trilinos_klu_zl_extract (trilinos_klu_l_numeric *, trilinos_klu_l_symbolic *,
    UF_long *, UF_long *, double *, double *,
    UF_long *, UF_long *, double *, double *,
    UF_long *, UF_long *, double *, double *,
    UF_long *, UF_long *, double *, UF_long *, trilinos_klu_l_common *) ;


/* -------------------------------------------------------------------------- */
/* KLU memory management routines */
/* -------------------------------------------------------------------------- */

void *trilinos_klu_malloc	/* returns pointer to the newly malloc'd block */
(
    /* ---- input ---- */
    size_t n,		/* number of items */
    size_t size,	/* size of each item */
    /* --------------- */
    trilinos_klu_common *Common
) ;

void *trilinos_klu_free		/* always returns NULL */
(
    /* ---- in/out --- */
    void *p,		/* block of memory to free */
    size_t n,		/* number of items */
    size_t size,	/* size of each item */
    /* --------------- */
    trilinos_klu_common *Common
) ;

void *trilinos_klu_realloc	/* returns pointer to reallocated block */
(
    /* ---- input ---- */
    size_t nnew,	/* requested # of items in reallocated block */
    size_t nold,	/* current size of block, in # of items */
    size_t size,	/* size of each item */
    /* ---- in/out --- */
    void *p,		/* block of memory to realloc */
    /* --------------- */
    trilinos_klu_common *Common
) ;

void *trilinos_klu_l_malloc (size_t, size_t, trilinos_klu_l_common *) ;
void *trilinos_klu_l_free (void *, size_t, size_t, trilinos_klu_l_common *) ;
void *trilinos_klu_l_realloc (size_t, size_t, size_t, void *, trilinos_klu_l_common *) ;


/* ========================================================================== */
/* === KLU version ========================================================== */
/* ========================================================================== */

/* All versions of KLU include these definitions.
 * As an example, to test if the version you are using is 1.2 or later:
 *
 *	if (TRILINOS_KLU_VERSION >= TRILINOS_KLU_VERSION_CODE (1,2)) ...
 *
 * This also works during compile-time:
 *
 *	#if (KLU >= TRILINOS_KLU_VERSION_CODE (1,2))
 *	    printf ("This is version 1.2 or later\n") ;
 *	#else
 *	    printf ("This is an early version\n") ;
 *	#endif
 */

#define TRILINOS_KLU_DATE "May 31, 2007"
#define TRILINOS_KLU_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define TRILINOS_KLU_MAIN_VERSION 1
#define TRILINOS_KLU_SUB_VERSION 0
#define TRILINOS_KLU_SUBSUB_VERSION 0
#define TRILINOS_KLU_VERSION TRILINOS_KLU_VERSION_CODE(TRILINOS_KLU_MAIN_VERSION,TRILINOS_KLU_SUB_VERSION)

#ifdef __cplusplus
}
#endif
#endif
