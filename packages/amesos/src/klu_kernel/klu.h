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

    /* outputs, allocated on output (or NULL if error occurs) */
    int **p_Lp,	    /* Column pointers for L, of size n+1 */
    int **p_Li,	    /* row indices for L */
    double **p_Lx,  /* values of L */
    int **p_Up,	    /* Column pointers for U, of size n+1 */
    int **p_Ui,	    /* row indices for U */
    double **p_Ux,  /* values of U */
    int **p_P,	    /* row permutation */

    /* scalar outputs */
    int *p_noffdiag,	    /* # of off-diagonal pivots chosen */
    double *p_umin,
    double *p_umax,

    /* workspace, ignored if NULL */
    double *X,	    /* size n double's, if present.  Zero on output */
    int *Work	    /* size 5n int's, if present */
) ;


void klu_free
(
    int **p_Lp,
    int **p_Li,
    double **p_Lx,
    int **p_Up,
    int **p_Ui,
    double **p_Ux,
    int **p_P
) ;


void klu_lsolve
(
    /* inputs, not modified: */
    int n,
    int Lp [ ],
    int Li [ ],
    double Lx [ ],
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
    /* right-hand-side on input, solution to Ux=b on output */
    double X [ ]
) ;


void klu_permute
(
    /* inputs, not modified: */
    int n,
    int P [ ],
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
