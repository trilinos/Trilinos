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
    int
	n,
	*P, 		/* size n */
	*Q,		/* size n */
	*R,		/* size n+1, but only R [0..nblocks] is used */
	nzoff,		/* nz in off-diagonal blocks */
	nblocks,	/* number of blocks */
	maxnz,		/* max nz in any block */
	maxblock,	/* size of largest block */
	ordering ;	/* ordering used (AMD or COLAMD) */

    /* this info is stored as double, to avoid integer overflow: */
    double lnz, unz ;	/* estimated nz in L and U, including diagonals */
    double *Lnz ;	/* size n, but only Lnz [0..nblocks-1] is used */

} klu_symbolic ;

/* -------------------------------------------------------------------------- */
/* Numeric object - contains the factors computed by klu_btf_factor */
/* -------------------------------------------------------------------------- */

typedef struct
{
    /* LU factors of each block, the pivot row permutation, and the
     * entries in the off-diagonal blocks */

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
    double *X ;		/* size n workspace, all zero */

    /* off-diagonal entries */
    int *Offp, *Offi ;
    double *Offx ;

} klu_numeric ;

/* -------------------------------------------------------------------------- */
/* klu_btf_analyze:  pre-orderings and analyzes a matrix with BTF and AMD */
/* -------------------------------------------------------------------------- */

klu_symbolic *klu_btf_analyze
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Control [ ],
    double Info [ ]
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
    double Control [ ],
    double Info [ ]
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

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ],

    /* workspace of size n, undefined on input and output */
    double W [ ]
) ;

/* -------------------------------------------------------------------------- */
/* klu_btf_scale: computes scale factors (user doesn't need to call directly) */
/* -------------------------------------------------------------------------- */

int klu_btf_scale
(
    /* inputs, not modified */
    int n,
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    /* outputs, not defined on input */
    double Rs [ ],
    /* workspace, not defined on input or output */
    int W [ ]		/* size n */
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
    /* input, and numerical values modified on output */
    klu_numeric *Numeric
) ;

/* -------------------------------------------------------------------------- */
/* KLU Control array */
/* -------------------------------------------------------------------------- */

#define KLU_BTF_CONTROL 20	    /* size of Control array */

/* contents of Control */
#define KLU_BTF_CONTROL_PRL 0		    /* print level */
#define KLU_BTF_CONTROL_PRL_DEFAULT 0

/* used in klu_btf_analyze */
#define KLU_BTF_CONTROL_BTF 1		    /* selecting BTF */
#define KLU_BTF_CONTROL_BTF_DEFAULT 1	    /* (default is to use BTF) */

#define KLU_BTF_CONTROL_AMD_DENSE 2	    /* AMD dense control parameter */
#define KLU_BTF_CONTROL_AMD_DENSE_DEFAULT AMD_DEFAULT_DENSE

#define KLU_BTF_CONTROL_ORDERING 3	    /* AMD: 0, COLAMD: 1. default: AMD*/
#define KLU_BTF_CONTROL_USE_AMD 0
#define KLU_BTF_CONTROL_USE_COLAMD 1
#define KLU_BTF_CONTROL_ORDERING_DEFAULT KLU_BTF_CONTROL_USE_AMD

/* used in klu_btf_factor */
#define KLU_BTF_CONTROL_TOL 4		    /* default pivot tolerance */
#define KLU_BTF_CONTROL_TOL_DEFAULT 0.001

#define KLU_BTF_CONTROL_GROWTH 5	    /* realloc growth size */
#define KLU_BTF_CONTROL_GROWTH_DEFAULT 1.5  /* default realloc growth size */

#define KLU_BTF_CONTROL_INITMEM_AMD 6	    /* initial memory size (w/ AMD) */
#define KLU_BTF_CONTROL_INITMEM_AMD_DEFAULT 1.2	/* 1.2 times nnz(L) + n */

#define KLU_BTF_CONTROL_INITMEM_COLAMD 7    /* initial memory size (w/ COLAMD)*/
#define KLU_BTF_CONTROL_INITMEM_COLAMD_DEFAULT 10	/* 10 times nnz(A) + n*/

/* -------------------------------------------------------------------------- */
/* klu_btf_defaults: sets default control parameters */
/* -------------------------------------------------------------------------- */

void klu_btf_defaults
(
    double Control [KLU_BTF_CONTROL]
) ;

/* -------------------------------------------------------------------------- */
/* KLU Info array */
/* -------------------------------------------------------------------------- */

#define KLU_BTF_INFO 90		    /* size of Info array */

/* returned by all klu_btf_* routines: */
#define KLU_BTF_INFO_STATUS 0	    /* KLU_OK, KLU_OUT_OF_MEMORY, ... */

/* determined in klu_btf_analyze: */
#define KLU_BTF_INFO_N 1	    /* n, dimension of input matrix */
#define KLU_BTF_INFO_NZ 2	    /* # entries in input matrix */
#define KLU_BTF_INFO_NBLOCKS 3	    /* # of blocks in BTF form */
#define KLU_BTF_INFO_MAXBLOCK 4	    /* dimension of largest block */
#define KLU_BTF_INFO_MAXNZ 6	    /* max nz in any block of A */
#define KLU_BTF_INFO_NZOFF 7	    /* nz in off-diagonal blocks of A */
#define KLU_BTF_INFO_SYMMETRY 8	    /* symmetry of largest block */
#define KLU_BTF_INFO_ATIME 9	    /* analyze time */
#define KLU_BTF_INFO_EST_LNZ 10	    /* nz in L, estimated (incl. diagonal) */
#define KLU_BTF_INFO_EST_UNZ 11	    /* nz in U, estimated (incl. diagonal) */
#define KLU_BTF_INFO_EST_FLOPS 12   /* est. factorization flop count */

/* 10..30 unused */

/* determined in klu_btf_factor: */
#define KLU_BTF_INFO_LNZ 30	    /* nz in L, actual (incl. diagonal) */
#define KLU_BTF_INFO_UNZ 31	    /* nz in U, actual (incl. diagonal) */
#define KLU_BTF_INFO_FLOPS 32	    /* actual factorization flop count */
#define KLU_BTF_INFO_UMIN 33	    /* min abs diagonal entry in U */
#define KLU_BTF_INFO_UMAX 34	    /* max abs diagonal entry in U */
#define KLU_BTF_INFO_REALLOC 35	    /* # of reallocations of L and/or U */
#define KLU_BTF_INFO_NOFFDIAG 36    /* number of off-diagonal pivots */
#define KLU_BTF_INFO_FTIME 37	    /* factorize time */

#define KLU_BTF_INFO_F2TIME 60	    /* refactorize time */

#define KLU_BTF_INFO_STIME 80	    /* solve time */

/* 37..89 unused */

#endif
