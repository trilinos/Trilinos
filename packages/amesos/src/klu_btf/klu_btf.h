/* ========================================================================== */
/* === klu_btf include file ================================================= */
/* ========================================================================== */

/* Include file for user programs that call klu_btf_* routines */

#ifndef _KLU_BTF_H
#define _KLU_BTF_H

#include "klu.h"
#include "amd.h"
#include "colamd.h"

/* -------------------------------------------------------------------------- */
/* Symbolic object - contains the pre-ordering computed by klu_btf_analyze */
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
	ordering,	/* ordering used (AMD, COLAMD, or GIVEN) */
	do_btf ;	/* whether or not BTF preordering was requested */

} klu_symbolic ;

/* -------------------------------------------------------------------------- */
/* Numeric object - contains the factors computed by klu_btf_factor */
/* -------------------------------------------------------------------------- */

typedef struct
{
    /* LU factors of each block, the pivot row permutation, and the
     * entries in the off-diagonal blocks */

    double umin ;	/* min abs diagonal entry in U */
    double umax ;	/* max abs diagonal entry in U */

    int nblocks ;
    int lnz, unz ;	/* actual nz in L and U, including diagonals */
    int *Pnum ;		/* final pivot permutation */
    int *Pinv ;		/* inverse of final pivot permutation */
    int noffdiag ;	/* # of off-diagonal pivots */

    /* LU factors of each block */
    int
	**Lbp,		/* Lbp [k] is pointer to Lp array for kth block */
	**Lbi,
	**Ubp,
	**Ubi ;
    double
	**Lbx,
	**Ubx ;

    double *Singleton ;	/* singleton values */

    double *Rs ;	/* row scaling factors */
    int scale ;		/* 0: none, 1: sum, 2: max */

    /* permanent workspace for factorization and solve */
    size_t worksize ;		/* size (in bytes) of Work */
    double *Work ;		/* of size MAX (4n doubles,
				   n doubles + 5*maxblock int's */
    double *Xwork ;		/* aliases into Numeric->Work */
    int *Iwork ;

    /* off-diagonal entries */
    int *Offp, *Offi ;
    double *Offx ;

    /* statistics determined in klu_btf_factor: */
    /* double flops ;	TODO: actual factorization flop count */
    int nlrealloc ;	/* # of reallocations of L */
    int nurealloc ;	/* # of reallocations of U */

} klu_numeric ;

/* -------------------------------------------------------------------------- */
/* KLU control parameters */
/* -------------------------------------------------------------------------- */

typedef struct
{
    double tol ;		/* pivot tolerance for diagonal preference */
    double growth ;		/* realloc growth size */
    double initmem_amd ;	/* init. memory size with AMD: c*nnz(L) + n */
    double initmem ;		/* init. memory size: c*nnz(A) + n */
    int btf ;			/* use BTF pre-ordering, or not */
    int ordering ;		/* 0: AMD, 1: COLAMD, 2: user P and Q */
    int scale ;			/* row scaling: 0: none, 1: sum, 2: max */

} klu_control ;

/* -------------------------------------------------------------------------- */
/* klu_btf_analyze:  pre-orderings and analyzes a matrix with BTF and AMD */
/* -------------------------------------------------------------------------- */

/* Order the matrix with BTF (or not), then AMD, COLAMD, or natural ordering
 * on the blocks. */

klu_symbolic *klu_btf_analyze
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    klu_control *control    /* optional; may be NULL */
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_analyze_given: analyzes a matrix using given P and Q */
/* -------------------------------------------------------------------------- */

/* Order the matrix with BTF (or not), then use natural or given ordering
 * Puser and Quser on the blocks.  Puser and Quser are interpretted as identity
 * if NULL. */

klu_symbolic *klu_btf_analyze_given
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    int Puser [ ],	/* size n, user's row permutation (may be NULL) */
    int Quser [ ],	/* size n, user's column permutation (may be NULL) */
    klu_control *control    /* optional; may be NULL */
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_factor:  factors a matrix using the klu_btf_analyze results */
/* -------------------------------------------------------------------------- */

klu_numeric *klu_btf_factor
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    klu_symbolic *Symbolic,
    klu_control *control    /* optional; may be NULL */
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_free_symbolic: destroys the Symbolic object */
/* -------------------------------------------------------------------------- */

void klu_btf_free_symbolic
(
    klu_symbolic **Symbolic
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_free_numeric: destroys the Numeric object */
/* -------------------------------------------------------------------------- */

void klu_btf_free_numeric
(
    klu_numeric **Numeric
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_solve: solves Ax=b using the Symbolic and Numeric */
/* -------------------------------------------------------------------------- */

void klu_btf_solve
(
    /* inputs, not modified */
    klu_symbolic *Symbolic,
    klu_numeric *Numeric,
    int ldim,		    /* leading dimension of B */
    int nrhs,		    /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ]
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_tsolve: solves A'x=b using the Symbolic and Numeric */
/* -------------------------------------------------------------------------- */

void klu_btf_tsolve
(
    /* inputs, not modified */
    klu_symbolic *Symbolic,
    klu_numeric *Numeric,
    int ldim,		    /* leading dimension of B */
    int nrhs,		    /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ]
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_scale: computes scale factors (user doesn't need to call directly) */
/* -------------------------------------------------------------------------- */

int klu_btf_scale
(
    /* inputs, not modified */
    int scale,		/* row scaling method: 0: none, 1: sum, 2: max */
    int n,
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    /* outputs, not defined on input */
    double Rs [ ],
    /* workspace, not defined on input or output */
    int W [ ]		/* size n, if present (can be a NULL pointer) */
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_refactor: refactorizes matrix with same ordering as klu_btf_factor */
/* -------------------------------------------------------------------------- */

int klu_btf_refactor	/* returns KLU_OK if OK, < 0 if error */
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    klu_symbolic *Symbolic,
    klu_control *control,	/* optional; may be NULL */
    /* input, and numerical values modified on output */
    klu_numeric *Numeric
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_defaults: sets default control parameters */
/* -------------------------------------------------------------------------- */

void klu_btf_defaults
(
    klu_control *control	/* optional; may be NULL */
) ;

#endif
