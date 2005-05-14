/* ========================================================================== */
/* === colamd/symamd prototypes and definitions ============================= */
/* ========================================================================== */

/*
 * CCOLAMD version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/*
    You must include this file (colamd.h) in any routine that uses colamd,
    symamd, or the related macros and definitions.

    Authors:

	The authors of COLAMD are Stefan I. Larimore and Timothy A.
	Davis (davis at cise.ufl.edu), University of Florida.  The algorithm was
	developed in collaboration with John Gilbert, Xerox PARC, and Esmond
	Ng, Oak Ridge National Laboratory.  CCOLAMD is an extension of COLAMD,
	by Timothy A. Davis and Sivasankaran Rajamanickam

    Acknowledgements:

	This work was supported by the National Science Foundation, under
	grants DMS-9504974 and DMS-9803599.
TODO: add Sandia support here.

    Notice:

	Copyright (c) 1998-2005 by the University of Florida.
	All Rights Reserved.

*/

#ifndef CCOLAMD_H
#define CCOLAMD_H

#include <stdlib.h>

/* ========================================================================== */
/* === Knob and statistics definitions ====================================== */
/* ========================================================================== */

/* size of the knobs [ ] array.  Only knobs [0..1] are currently used. */
#define CCOLAMD_KNOBS 20

/* number of output statistics.  Only stats [0..6] are currently used. */
#define CCOLAMD_STATS 20

/* knobs [0] and stats [0]: dense row knob and output statistic. */
#define CCOLAMD_DENSE_ROW 0

/* knobs [1] and stats [1]: dense column knob and output statistic. */
#define CCOLAMD_DENSE_COL 1

/* ------------------------------- */
/* added for ccolamd from UMFPACK: */

/* knobs [2]: aggressive absorption option */
#define CCOLAMD_AGGRESSIVE 2

/* knobs [3]: LU or cholesky factorization option */
#define CCOLAMD_FACT_TYPE 3

/* ------------------------------- */

/* stats [2]: memory defragmentation count output statistic */
#define CCOLAMD_DEFRAG_COUNT 2

/* stats [3]: colamd status:  zero OK, > 0 warning or notice, < 0 error */
#define CCOLAMD_STATUS 3

/* stats [4..6]: error info, or info on jumbled columns */ 
#define CCOLAMD_INFO1 4
#define CCOLAMD_INFO2 5
#define CCOLAMD_INFO3 6

/* ------------------------------- */
/* added for ccolamd from UMFPACK: */
/* stats [7]: number of originally empty rows */
#define CCOLAMD_EMPTY_ROW 7
/* stats [8]: number of originally empty cols */
#define CCOLAMD_EMPTY_COL 8
/* stats [9]: number of rows with entries only in dense cols */
#define CCOLAMD_NEWLY_EMPTY_ROW 9
/* stats [10]: number of cols with entries only in dense rows */
#define CCOLAMD_NEWLY_EMPTY_COL 10

/* ------------------------------- */

/* error codes returned in stats [3]: */
#define CCOLAMD_OK				(0)
#define CCOLAMD_OK_BUT_JUMBLED			(1)
#define CCOLAMD_ERROR_A_not_present		(-1)
#define CCOLAMD_ERROR_p_not_present		(-2)
#define CCOLAMD_ERROR_nrow_negative		(-3)
#define CCOLAMD_ERROR_ncol_negative		(-4)
#define CCOLAMD_ERROR_nnz_negative		(-5)
#define CCOLAMD_ERROR_p0_nonzero		(-6)
#define CCOLAMD_ERROR_A_too_small		(-7)
#define CCOLAMD_ERROR_col_length_negative	(-8)
#define CCOLAMD_ERROR_row_index_out_of_bounds	(-9)
#define CCOLAMD_ERROR_out_of_memory		(-10)
#define CCOLAMD_ERROR_invalid_cset		(-11)
#define CCOLAMD_ERROR_internal_error		(-999)

/* ========================================================================== */
/* === Prototypes of user-callable routines ================================= */
/* ========================================================================== */

int ccolamd_recommended		/* returns recommended value of Alen, */
				/* or (-1) if input arguments are erroneous */
(
    int nnz,			/* nonzeros in A */
    int n_row,			/* number of rows in A */
    int n_col			/* number of columns in A */
) ;

long ccolamd_l_recommended		/* returns recommended value of Alen, */
				/* or (-1) if input arguments are erroneous */
(
    long nnz,			/* nonzeros in A */
    long n_row,			/* number of rows in A */
    long n_col			/* number of columns in A */
) ;

void ccolamd_set_defaults	/* sets default parameters */
(				/* knobs argument is modified on output */
    double knobs [CCOLAMD_KNOBS]	/* parameter settings for colamd */
) ;

void ccolamd_l_set_defaults	/* sets default parameters */
(				/* knobs argument is modified on output */
    double knobs [CCOLAMD_KNOBS]	/* parameter settings for colamd */
) ;

int umf_ccolamd
(				/* A and p arguments are modified on output */
    int n_row,			/* number of rows in A */
    int n_col,			/* number of columns in A */
    int Alen,			/* size of the array A */
    int A [],			/* row indices of A, of size Alen */
    int p [],			/* column pointers of A, of size n_col+1 */
    double knobs [CCOLAMD_KNOBS],/* parameter settings for colamd */
    int stats [CCOLAMD_STATS]	/* colamd output statistics and error codes */

    /* ------------------ */
    /* added for UMFPACK: each Front_ array is of size n_col+1 */
    , int Front_npivcol [ ]     /* # pivot cols in each front */
    , int Front_nrows [ ]       /* # of rows in each front (incl. pivot rows) */
    , int Front_ncols [ ]       /* # of cols in each front (incl. pivot cols) */
    , int Front_parent [ ]      /* parent of each front */
    , int Front_cols [ ]        /* link list of pivot columns for each front */
    , int *p_nfr                /* total number of frontal matrices */
    , int InFront [ ]           /* InFront [row] = f if the original row was */
    /* ------------------ */

    /* ------------------- */
    /* added for ccolamd : */
    , int cset []		/* Constraint set of A */
    /* ------------------- */
    , int max_col_degree
) ;

long umf_ccolamd_l
(				/* A and p arguments are modified on output */
    long n_row,			/* number of rows in A */
    long n_col,			/* number of columns in A */
    long Alen,			/* size of the array A */
    long A [],			/* row indices of A, of size Alen */
    long p [],			/* column pointers of A, of size n_col+1 */
    double knobs [CCOLAMD_KNOBS],/* parameter settings for colamd */
    long stats [CCOLAMD_STATS]	/* colamd output statistics and error codes */

    /* ------------------ */
    /* added for UMFPACK: each Front_ array is of size n_col+1 */
    , long Front_npivcol [ ]     /* # pivot cols in each front */
    , long Front_nrows [ ]       /* # of rows in each front (incl. pivot rows) */
    , long Front_ncols [ ]       /* # of cols in each front (incl. pivot cols) */
    , long Front_parent [ ]      /* parent of each front */
    , long Front_cols [ ]        /* link list of pivot columns for each front */
    , long *p_nfr                /* total number of frontal matrices */
    , long InFront [ ]           /* InFront [row] = f if the original row was */
    /* ------------------ */

    /* ------------------- */
    /* added for ccolamd : */
    , long cset []		/* Constraint set of A */
    /* ------------------- */
    , long max_col_degree
) ;

int ccolamd		/* returns (1) if successful, (0) otherwise*/
(				/* A and p arguments are modified on output */
    int n_row,			/* number of rows in A */
    int n_col,			/* number of columns in A */
    int Alen,			/* size of the array A */
    int A [],			/* row indices of A, of size Alen */
    int p [],			/* column pointers of A, of size n_col+1 */
    double knobs [CCOLAMD_KNOBS],/* parameter settings for colamd */
    int stats [CCOLAMD_STATS]	/* colamd output statistics and error codes */
    /* ------------------- */
    /* added for ccolamd : */
    , int cset []		/* Constraint set of A */
    /* ------------------- */
) ;

long ccolamd_l		/* returns (1) if successful, (0) otherwise*/
(				/* A and p arguments are modified on output */
    long n_row,			/* number of rows in A */
    long n_col,			/* number of columns in A */
    long Alen,			/* size of the array A */
    long A [],			/* row indices of A, of size Alen */
    long p [],			/* column pointers of A, of size n_col+1 */
    double knobs [CCOLAMD_KNOBS],/* parameter settings for colamd */
    long stats [CCOLAMD_STATS]	/* colamd output statistics and error codes */
    /* ------------------- */
    /* added for ccolamd : */
    , long cset []		/* Constraint set of A */
    /* ------------------- */
) ;

int csymamd			/* return (1) if OK, (0) otherwise */
(
    int n,				/* number of rows and columns of A */
    int A [],				/* row indices of A */
    int p [],				/* column pointers of A */
    int perm [],			/* output permutation, size n_col+1 */
    double knobs [CCOLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    int stats [CCOLAMD_STATS],		/* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
    /* ------------------- */
    /* added for ccolamd : */
    , int cset []			/* Constraint set of A */
    , int symmetry			/* 0: both, >0: upper, <0: lower */
    /* ------------------- */
) ;

long csymamd_l			/* return (1) if OK, (0) otherwise */
(
    long n,				/* number of rows and columns of A */
    long A [],				/* row indices of A */
    long p [],				/* column pointers of A */
    long perm [],			/* output permutation, size n_col+1 */
    double knobs [CCOLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    long stats [CCOLAMD_STATS],		/* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
    /* ------------------- */
    /* added for ccolamd : */
    , long cset []			/* Constraint set of A */
    , long symmetry			/* 0: both, >0: upper, <0: lower */
    /* ------------------- */
) ;

void ccolamd_report
(
    int stats [CCOLAMD_STATS]
) ;

void ccolamd_l_report
(
    long stats [CCOLAMD_STATS]
) ;

void csymamd_report
(
    int stats [CCOLAMD_STATS]
) ;

void csymamd_l_report
(
    long stats [CCOLAMD_STATS]
) ;

void ccolamd_apply_order
(
    int Front [ ],	    /* of size nn on input, size nfr on output */
    const int Order [ ],    /* Order [i] = k, i in the range 0..nn-1,
			     * and k in the range 0..nfr-1, means that node
			     * i is the kth node in the postordered tree. */
    int Temp [ ],	    /* workspace of size nfr */
    int nn,		    /* nodes are numbered in the range 0..nn-1 */
    int nfr		    /* the number of nodes actually in use */
) ;

void ccolamd_l_apply_order
(
    long Front [ ],	    /* of size nn on input, size nfr on output */
    const long Order [ ],    /* Order [i] = k, i in the range 0..nn-1,
			     * and k in the range 0..nfr-1, means that node
			     * i is the kth node in the postordered tree. */
    long Temp [ ],	    /* workspace of size nfr */
    long nn,		    /* nodes are numbered in the range 0..nn-1 */
    long nfr		    /* the number of nodes actually in use */
) ;


void ccolamd_fsize
(
    int nn,
    int MaxFsize [ ],
    int Fnrows [ ],
    int Fncols [ ],
    int Parent [ ],
    int Npiv [ ]
) ;

void ccolamd_l_fsize
(
    long nn,
    long MaxFsize [ ],
    long Fnrows [ ],
    long Fncols [ ],
    long Parent [ ],
    long Npiv [ ]
) ;

void ccolamd_postorder
(
    int nn,
    int Parent [ ],
    int Npiv [ ],
    int Fsize [ ],
    int Order [ ],
    int Child [ ],
    int Sibling [ ],
    int Stack [ ],
    int Front_cols [],
    int in_cset []
) ;

void ccolamd_l_postorder
(
    long nn,
    long Parent [ ],
    long Npiv [ ],
    long Fsize [ ],
    long Order [ ],
    long Child [ ],
    long Sibling [ ],
    long Stack [ ],
    long Front_cols [],
    long in_cset []
) ;

int ccolamd_post_tree
(
    int root,
    int k,
    int Child [ ],
    const int Sibling [ ],
    int Order [ ],
    int Stack [ ]
) ;

long ccolamd_l_post_tree
(
    long root,
    long k,
    long Child [ ],
    const long Sibling [ ],
    long Order [ ],
    long Stack [ ]
) ;


#endif /* CCOLAMD_H */
