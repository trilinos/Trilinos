/* ========================================================================== */
/* === klu include file ===================================================== */
/* ========================================================================== */

/* For inclusion in user routines that call klu. */

#ifndef _KLU_H
#define _KLU_H

int klu_factor	/* returns 0 if OK, negative if error */
(
    /* inputs, not modified */
    int n,	    /* A is n-by-n. n must be > 0. */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    int Q [ ],	    /* size n, optional column permutation */
    double Control [ ],	    /* Control parameters (optional) */

    /* outputs, not defined on input */
    int Lp [ ],	    /* Column pointers for L, of size n+1 */
    int **p_Li,	    /* row indices for L */
    double **p_Lx,  /* values of L */
    int Up [ ],	    /* Column pointers for U, of size n+1 */
    int **p_Ui,	    /* row indices for U */
    double **p_Ux,  /* values of U */
    int P [ ],	    /* row permutation, size n */

    /* scalar outputs */
    int *p_noffdiag,	    /* # of off-diagonal pivots chosen */
    double *p_umin,
    double *p_umax,

    /* workspace, undefined on input */
    double *X,	    /* size n double's.  Zero on output */
    int *Work,	    /* size 5n int's */

    /* ---- the following are only used in the BTF case --- */

    /* inputs, not modified on output */
    int k1,	    /* the block of A is from k1 to k2-1 */
    int PSinv [ ],  /* inverse of P from symbolic factorization */
    double Rs [ ],  /* scale factors for A */

    /* inputs, modified on output */
    int Offp [ ],   /* off-diagonal matrix (modified by this routine) */
    int Offi [ ],
    double Offx [ ]
) ;


void klu_free
(
    int **p_Li,
    double **p_Lx,
    int **p_Ui,
    double **p_Ux
) ;


void klu_lsolve
(
    /* inputs, not modified: */
    int n,
    int Lp [ ],
    int Li [ ],
    double Lx [ ],
    int ldim,
    int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    double X [ ]
) ;


void klu_usolve
(
    /* inputs, not modified: */
    int n,
    int Up [ ],
    int Ui [ ],
    double Ux [ ],
    int ldim,
    int nrhs,
    /* right-hand-side on input, solution to Ux=b on output */
    double X [ ]
) ;


void klu_permute
(
    /* inputs, not modified: */
    int n,
    int P [ ],
    int ldim,
    int nrhs,
    double B [ ],
    /* output */
    double X [ ]
) ;


void klu_defaults
(
    double Control [ ]
) ;


/* Control parameters */
#define KLU_TOL 0	    /* partial pivoting tolerance */
#define KLU_LSIZE 1	    /* initial size of L */
#define KLU_USIZE 2	    /* initial size of U */
#define KLU_GROWTH 3	    /* memory growth factor */
#define KLU_CONTROL 10	    /* size of Control array */

/* return values of klu */
#define KLU_OK 0
#define KLU_SINGULAR (-1)
#define KLU_OUT_OF_MEMORY (-2)
#define KLU_INVALID (-3)

#endif
