/* ========================================================================== */
/* === klu_btf include file ================================================= */
/* ========================================================================== */

/* Include file for user programs that call klu_btf_* routines */

#include "klu.h"
#include "amd.h"

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
	maxblock ;	/* size of largest block */

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
    int Ai [ ]		/* size nz, row indices */
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
    double tol,
    klu_symbolic *Symbolic
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
    double Rs [ ]
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
