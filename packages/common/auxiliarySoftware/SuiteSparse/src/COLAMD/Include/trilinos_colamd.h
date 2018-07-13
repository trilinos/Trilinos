/* ========================================================================== */
/* === trilinos_colamd/symamd prototypes and definitions ============================= */
/* ========================================================================== */

/* TRILINOS_COLAMD / SYMAMD include file

    You must include this file (trilinos_colamd.h) in any routine that uses trilinos_colamd,
    symamd, or the related macros and definitions.

    Authors:

	The authors of the code itself are Stefan I. Larimore and Timothy A.
	Davis (davis at cise.ufl.edu), University of Florida.  The algorithm was
	developed in collaboration with John Gilbert, Xerox PARC, and Esmond
	Ng, Oak Ridge National Laboratory.

    Acknowledgements:

	This work was supported by the National Science Foundation, under
	grants DMS-9504974 and DMS-9803599.

    Notice:

	Copyright (c) 1998-2007, Timothy A. Davis, All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

	Permission is hereby granted to use, copy, modify, and/or distribute
	this program, provided that the Copyright, this License, and the
	Availability of the original version is retained on all copies and made
	accessible to the end-user of any code or package that includes TRILINOS_COLAMD
	or any modified version of TRILINOS_COLAMD. 

    Availability:

	The trilinos_colamd/symamd library is available at

	    http://www.cise.ufl.edu/research/sparse/trilinos_colamd/

	This is the http://www.cise.ufl.edu/research/sparse/trilinos_colamd/colamd.h
	file.  It is required by the trilinos_colamd.c, colamdmex.c, and symamdmex.c
	files, and by any C code that calls the routines whose prototypes are
	listed below, or that uses the trilinos_colamd/symamd definitions listed below.

*/

#ifndef TRILINOS_COLAMD_H
#define TRILINOS_COLAMD_H

/* make it easy for C++ programs to include TRILINOS_COLAMD */
#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include <stdlib.h>

/* ========================================================================== */
/* === TRILINOS_COLAMD version ======================================================= */
/* ========================================================================== */

/* TRILINOS_COLAMD Version 2.4 and later will include the following definitions.
 * As an example, to test if the version you are using is 2.4 or later:
 *
 * #ifdef TRILINOS_COLAMD_VERSION
 *	if (TRILINOS_COLAMD_VERSION >= TRILINOS_COLAMD_VERSION_CODE (2,4)) ...
 * #endif
 *
 * This also works during compile-time:
 *
 *  #if defined(TRILINOS_COLAMD_VERSION) && (COLAMD_VERSION >= TRILINOS_COLAMD_VERSION_CODE (2,4))
 *    printf ("This is version 2.4 or later\n") ;
 *  #else
 *    printf ("This is an early version\n") ;
 *  #endif
 *
 * Versions 2.3 and earlier of TRILINOS_COLAMD do not include a #define'd version number.
 */

#define TRILINOS_COLAMD_DATE "May 31, 2007"
#define TRILINOS_COLAMD_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define TRILINOS_COLAMD_MAIN_VERSION 2
#define TRILINOS_COLAMD_SUB_VERSION 7
#define TRILINOS_COLAMD_SUBSUB_VERSION 0
#define TRILINOS_COLAMD_VERSION \
	TRILINOS_COLAMD_VERSION_CODE(TRILINOS_COLAMD_MAIN_VERSION,TRILINOS_COLAMD_SUB_VERSION)

/* ========================================================================== */
/* === Knob and statistics definitions ====================================== */
/* ========================================================================== */

/* size of the knobs [ ] array.  Only knobs [0..1] are currently used. */
#define TRILINOS_COLAMD_KNOBS 20

/* number of output statistics.  Only stats [0..6] are currently used. */
#define TRILINOS_COLAMD_STATS 20

/* knobs [0] and stats [0]: dense row knob and output statistic. */
#define TRILINOS_COLAMD_DENSE_ROW 0

/* knobs [1] and stats [1]: dense column knob and output statistic. */
#define TRILINOS_COLAMD_DENSE_COL 1

/* knobs [2]: aggressive absorption */
#define TRILINOS_COLAMD_AGGRESSIVE 2

/* stats [2]: memory defragmentation count output statistic */
#define TRILINOS_COLAMD_DEFRAG_COUNT 2

/* stats [3]: trilinos_colamd status:  zero OK, > 0 warning or notice, < 0 error */
#define TRILINOS_COLAMD_STATUS 3

/* stats [4..6]: error info, or info on jumbled columns */ 
#define TRILINOS_COLAMD_INFO1 4
#define TRILINOS_COLAMD_INFO2 5
#define TRILINOS_COLAMD_INFO3 6

/* error codes returned in stats [3]: */
#define TRILINOS_COLAMD_OK				(0)
#define TRILINOS_COLAMD_OK_BUT_JUMBLED			(1)
#define TRILINOS_COLAMD_ERROR_A_not_present		(-1)
#define TRILINOS_COLAMD_ERROR_p_not_present		(-2)
#define TRILINOS_COLAMD_ERROR_nrow_negative		(-3)
#define TRILINOS_COLAMD_ERROR_ncol_negative		(-4)
#define TRILINOS_COLAMD_ERROR_nnz_negative		(-5)
#define TRILINOS_COLAMD_ERROR_p0_nonzero			(-6)
#define TRILINOS_COLAMD_ERROR_A_too_small		(-7)
#define TRILINOS_COLAMD_ERROR_col_length_negative	(-8)
#define TRILINOS_COLAMD_ERROR_row_index_out_of_bounds	(-9)
#define TRILINOS_COLAMD_ERROR_out_of_memory		(-10)
#define TRILINOS_COLAMD_ERROR_internal_error		(-999)


/* ========================================================================== */
/* === Prototypes of user-callable routines ================================= */
/* ========================================================================== */

/* define UF_long */
#include "trilinos_UFconfig.h"

size_t trilinos_colamd_recommended	/* returns recommended value of Alen, */
				/* or 0 if input arguments are erroneous */
(
    int nnz,			/* nonzeros in A */
    int n_row,			/* number of rows in A */
    int n_col			/* number of columns in A */
) ;

size_t trilinos_colamd_l_recommended	/* returns recommended value of Alen, */
				/* or 0 if input arguments are erroneous */
(
    UF_long nnz,		/* nonzeros in A */
    UF_long n_row,		/* number of rows in A */
    UF_long n_col		/* number of columns in A */
) ;

void trilinos_colamd_set_defaults	/* sets default parameters */
(				/* knobs argument is modified on output */
    double knobs [TRILINOS_COLAMD_KNOBS]	/* parameter settings for trilinos_colamd */
) ;

void trilinos_colamd_l_set_defaults	/* sets default parameters */
(				/* knobs argument is modified on output */
    double knobs [TRILINOS_COLAMD_KNOBS]	/* parameter settings for trilinos_colamd */
) ;

int trilinos_colamd		/* returns (1) if successful, (0) otherwise*/
(				/* A and p arguments are modified on output */
    int n_row,			/* number of rows in A */
    int n_col,			/* number of columns in A */
    int Alen,			/* size of the array A */
    int A [],			/* row indices of A, of size Alen */
    int p [],			/* column pointers of A, of size n_col+1 */
    double knobs [TRILINOS_COLAMD_KNOBS],/* parameter settings for trilinos_colamd */
    int stats [TRILINOS_COLAMD_STATS]	/* trilinos_colamd output statistics and error codes */
) ;

UF_long trilinos_colamd_l		/* returns (1) if successful, (0) otherwise*/
(				/* A and p arguments are modified on output */
    UF_long n_row,		/* number of rows in A */
    UF_long n_col,		/* number of columns in A */
    UF_long Alen,		/* size of the array A */
    UF_long A [],		/* row indices of A, of size Alen */
    UF_long p [],		/* column pointers of A, of size n_col+1 */
    double knobs [TRILINOS_COLAMD_KNOBS],/* parameter settings for trilinos_colamd */
    UF_long stats [TRILINOS_COLAMD_STATS]/* trilinos_colamd output statistics and error codes */
) ;

int trilinos_symamd		        /* return (1) if OK, (0) otherwise */
(
    int n,				/* number of rows and columns of A */
    int A [],				/* row indices of A */
    int p [],				/* column pointers of A */
    int perm [],			/* output permutation, size n_col+1 */
    double knobs [TRILINOS_COLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    int stats [TRILINOS_COLAMD_STATS],		/* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
) ;

UF_long trilinos_symamd_l			/* return (1) if OK, (0) otherwise */
(
    UF_long n,				/* number of rows and columns of A */
    UF_long A [],			/* row indices of A */
    UF_long p [],			/* column pointers of A */
    UF_long perm [],			/* output permutation, size n_col+1 */
    double knobs [TRILINOS_COLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    UF_long stats [TRILINOS_COLAMD_STATS],	/* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
) ;

void trilinos_colamd_report
(
    int stats [TRILINOS_COLAMD_STATS]
) ;

void trilinos_colamd_l_report
(
    UF_long stats [TRILINOS_COLAMD_STATS]
) ;

void trilinos_symamd_report
(
    int stats [TRILINOS_COLAMD_STATS]
) ;

void trilinos_symamd_l_report
(
    UF_long stats [TRILINOS_COLAMD_STATS]
) ;

#ifndef EXTERN
#define EXTERN extern
#endif

EXTERN int (*trilinos_colamd_printf) (const char *, ...) ;

#ifdef __cplusplus
}
#endif

#endif /* TRILINOS_COLAMD_H */
