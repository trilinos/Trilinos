/* ====================== ccolamd internal definitions ====================== */


/*
 * CCOLAMD version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* --------------------------------------- */
/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include <stdlib.h>

/* ----------------- */
/* Added for ccolamd */
#include <math.h>
#include <limits.h>
/* ----------------- */

/* -------------------------------- */
/* Moved from ccolamd.c for ccolamd */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#else
#include <stdio.h>
#include <assert.h>
#endif /* MATLAB_MEX_FILE */
/* -------------------------------- */

/* ------------------------------- */
/* added for ccolamd integer type: int or long  */

#ifdef DLONG

#define Int long
#define ID "%ld"
#define Int_MAX LONG_MAX
#define Int_MIN LONG_MIN

#define CCOLAMD_f_recommended ccolamd_l_recommended
#define CCOLAMD_set_defaults ccolamd_l_set_defaults
#define CCOLAMD_umf_ccolamd_MAIN umf_ccolamd_l
#define CCOLAMD_MAIN ccolamd_l
#define CCOLAMD_apply_order ccolamd_l_apply_order
#define CCOLAMD_postorder ccolamd_l_postorder
#define CCOLAMD_fsize ccolamd_l_fsize
#define CSYMAMD_MAIN csymamd_l
#define CCOLAMD_report ccolamd_l_report
#define CSYMAMD_report csymamd_l_report

#else /* DLONG */

#define Int int
#define ID "%d"
#define Int_MAX INT_MAX
#define Int_MIN INT_MIN

#define CCOLAMD_f_recommended ccolamd_recommended
#define CCOLAMD_set_defaults ccolamd_set_defaults
#define CCOLAMD_umf_ccolamd_MAIN umf_ccolamd
#define CCOLAMD_MAIN ccolamd
#define CCOLAMD_apply_order ccolamd_apply_order
#define CCOLAMD_postorder ccolamd_postorder
#define CCOLAMD_fsize ccolamd_fsize
#define CSYMAMD_MAIN csymamd
#define CCOLAMD_report ccolamd_report
#define CSYMAMD_report csymamd_report

#endif /* DLONG */

/* ------------------------------- */



/* ========================================================================== */
/* === Row and Column structures ============================================ */
/* ========================================================================== */

/* User code that makes use of the colamd/symamd routines need not directly */
/* reference these structures.  They are used only for the COLAMD_RECOMMENDED */
/* macro. */

typedef struct Colamd_Col_struct
{
    Int start ;		/* index for A of first row in this column, or DEAD */
			/* if column is dead */
    Int length ;	/* number of rows in this column */
    union
    {
	Int thickness ;	/* number of original columns represented by this */
			/* col, if the column is alive */
	Int parent ;	/* parent in parent tree super-column structure, if */
			/* the column is dead */
    } shared1 ;
    union
    {
	Int score ;	
	Int order ;
    } shared2 ; 
    union
    {
	Int headhash ;	/* head of a hash bucket, if col is at the head of */
			/* a degree list */
	Int hash ;	/* hash value, if col is not in a degree list */
	Int prev ;	/* previous column in degree list, if col is in a */
			/* degree list (but not at the head of a degree list) */
    } shared3 ;
    union
    {
	Int degree_next ;	/* next column, if col is in a degree list */
	Int hash_next ;		/* next column, if col is in a hash list */
    } shared4 ;

    /* -------------------------- */
    /* added for ccolamd UMFPACK: */
    Int nextcol ;       /* next column in this supercolumn */
    Int lastcol ;       /* last column in this supercolumn */
    /* -------------------------- */


} Colamd_Col ;

	

typedef struct Colamd_Row_struct
{
    Int start ;		/* index for A of first col in this row */
    Int length ;	/* number of principal columns in this row */
    union
    {
	Int degree ;	/* number of principal & non-principal columns in row */
	Int p ;		/* used as a row pointer in init_rows_cols () */
    } shared1 ;
    union
    {
	Int mark ;	/* for computing set differences and marking dead rows*/
	Int first_column ;/* first column in row (used in garbage collection) */
    } shared2 ;

    /* -------------------------- */
    /* added for ccolamd UMFPACK: */
    Int thickness ;     /* number of original rows represented by this row */
                        /* that are not yet pivotal */
    Int front ;         /* -1 if an original row */
    			/* k if this row represents the kth frontal matrix */
                        /* where k goes from 0 to at most n_col-1 */
    /* -------------------------- */

} Colamd_Row ;

/* ========================================================================== */
/* === Colamd recommended memory size ======================================= */
/* ========================================================================== */

/*
    The recommended length Alen of the array A passed to colamd is given by
    the COLAMD_RECOMMENDED (nnz, n_row, n_col) macro.  It returns -1 if any
    argument is negative.  2*nnz space is required for the row and column
    indices of the matrix. COLAMD_C (n_col) + COLAMD_R (n_row) space is
    required for the Col and Row arrays, respectively, which are internal to
    colamd.  An additional n_col space is the minimal amount of "elbow room",
    and nnz/5 more space is recommended for run time efficiency.

    This macro is not needed when using symamd.

    Explicit typecast to int added Sept. 23, 2002, COLAMD version 2.2, to avoid
    gcc -pedantic warning messages.
*/

#define CCOLAMD_C(n_col) ((Int) (((n_col) + 1) * sizeof (Colamd_Col) / sizeof (Int)))
#define CCOLAMD_R(n_row) ((Int) (((n_row) + 1) * sizeof (Colamd_Row) / sizeof (Int)))

/* -------------------- */
/* modified for ccolamd  */
#define CCOLAMD_RECOMMENDED(nnz, n_row, n_col)                                 \
(                                                                             \
((nnz) < 0 || (n_row) < 0 || (n_col) < 0)                                     \
?                                                                             \
    (-1)                                                                      \
:                                                                             \
    MAX (2 * (nnz), 4 * (n_col)) + \
    CCOLAMD_C (n_col) + CCOLAMD_R (n_row) + (n_col) + ((nnz) / 5) \
    + ((3 * n_col) + 1)  + 5 * (n_col + 1) + (n_row)  \
)
/* -------------------- */

/* ========================================================================== */
/* === [ Definitions added/moved for ccolamd  =============================== */
/* ========================================================================== */

/* Ensure that debugging is turned off: */
#ifndef NDEBUG
#define NDEBUG
#endif

/* turn on debugging by uncommenting the following line
#undef NDEBUG
*/

#define EMPTY (-1)
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* Routines are either PUBLIC (user-callable) or PRIVATE (not user-callable) */
#define GLOBAL 
#define PUBLIC
#define PRIVATE static 

/* ------------------------------- */
/* added for ccolamd from UMFPACK: */

#define UMFPACK_DENSE_DEGREE_THRESHOLD(alpha,n) \
    ((Int) MAX (16.0, (alpha) * 16.0 * sqrt ((double) (n))))

#define IN_CSET(c) ((in_cset == (int *) NULL) ? (0) : (in_cset [c]))

/* ------------------------------- */

/* -------------------------------------------------------------------------- */
/* [ Numerical relop macros for correctly handling the NaN case */
/* Added for ccolamd from umf_version.h */
/* -------------------------------------------------------------------------- */

/*
SCALAR_IS_NAN(x):
    True if x is NaN.  False otherwise.  The commonly-existing isnan(x)
    function could be used, but it's not in Kernighan & Ritchie 2nd edition
    (ANSI C).  It may appear in <math.h>, but I'm not certain about
    portability.  The expression x != x is true if and only if x is NaN,
    according to the IEEE 754 floating-point standard.

SCALAR_IS_ZERO(x):
    True if x is zero.  False if x is nonzero, NaN, or +/- Inf.
    This is (x == 0) if the compiler is IEEE 754 compliant.

SCALAR_IS_NONZERO(x):
    True if x is nonzero, NaN, or +/- Inf.  False if x zero.
    This is (x != 0) if the compiler is IEEE 754 compliant.

SCALAR_IS_LTZERO(x):
    True if x is < zero or -Inf.  False if x is >= 0, NaN, or +Inf.
    This is (x < 0) if the compiler is IEEE 754 compliant.
*/

#if defined (MATHWORKS)

/* The MathWorks has their own macros in util.h that handle NaN's properly. */
#define SCALAR_IS_NAN(x)        (utIsNaN (x))
#define SCALAR_IS_ZERO(x)       (utEQZero (x))
#define SCALAR_IS_NONZERO(x)    (utNEZero (x))
#define SCALAR_IS_LTZERO(x)     (utLTZero (x))

#elif defined (UMF_WINDOWS)

/* Yes, this is exceedingly ugly.  Blame Microsoft, which hopelessly */
/* violates the IEEE 754 floating-point standard in a bizarre way. */
/* If you're using an IEEE 754-compliant compiler, then x != x is true */
/* iff x is NaN.  For Microsoft, (x < x) is true iff x is NaN. */
/* So either way, this macro safely detects a NaN. */
#define SCALAR_IS_NAN(x)        (((x) != (x)) || (((x) < (x))))
#define SCALAR_IS_ZERO(x)       (((x) == 0.) && !SCALAR_IS_NAN(x))
#define SCALAR_IS_NONZERO(x)    (((x) != 0.) || SCALAR_IS_NAN(x))
#define SCALAR_IS_LTZERO(x)     (((x) < 0.) && !SCALAR_IS_NAN(x))

#else

/* These all work properly, according to the IEEE 754 standard ... except on */
/* a PC with windows.  Works fine in Linux on the same PC... */
#define SCALAR_IS_NAN(x)        ((x) != (x))
#define SCALAR_IS_ZERO(x)       ((x) == 0.)
#define SCALAR_IS_NONZERO(x)    ((x) != 0.)
#define SCALAR_IS_LTZERO(x)     ((x) < 0.)

#endif

/* scalar absolute value macro. If x is NaN, the result is NaN: */
#define SCALAR_ABS(x) ((SCALAR_IS_LTZERO (x)) ? -(x) : (x))

/* true if an integer (stored in double x) would overflow (or if x is NaN) */
#define INT_OVERFLOW(x) ((!((x) * (1.0+1e-8) <= (double) Int_MAX)) \
                        || SCALAR_IS_NAN (x))

/* print a scalar (avoid printing "-0" for negative zero).  */
#define PRINT_SCALAR(a) \
{ \
    if (SCALAR_IS_NONZERO (a)) \
    { \
        PRINTF ((" (%g)", (a))) ; \
    } \
    else \
    { \
        PRINTF ((" (0)")) ; \
    } \
}
/* -------------------------------------------------------------------------] */

/* ========================================================================== */
/* === Colamd reporting mechanism =========================================== */
/* ========================================================================== */
    
#ifdef MATLAB_MEX_FILE
    
/* use mexPrintf in a MATLAB mexFunction, for debugging and statistics output */
#define PRINTF mexPrintf
    
/* In MATLAB, matrices are 1-based to the user, but 0-based internally */
#define INDEX(i) ((i)+1)
    
#else
    
/* Use printf in standard C environment, for debugging and statistics output. */
/* Output is generated only if debugging is enabled at compile time, or if */
/* the caller explicitly calls colamd_report or symamd_report. */
#define PRINTF printf
    
/* In C, matrices are 0-based and indices are reported as such in *_report */
#define INDEX(i) (i)

#endif /* MATLAB_MEX_FILE */


/* ========================================================================== */
/* === [ Debugging prototypes and definitions =============================== */
/* ========================================================================== */

#ifndef NDEBUG
    
/* === colamd_debug is retained in ccolamd.c ================================ */
#define DEBUG0(params) { (void) PRINTF params ; }
#define DEBUG1(params) { if (colamd_debug >= 1) (void) PRINTF params ; }
#define DEBUG2(params) { if (colamd_debug >= 2) (void) PRINTF params ; }
#define DEBUG3(params) { if (colamd_debug >= 3) (void) PRINTF params ; }
#define DEBUG4(params) { if (colamd_debug >= 4) (void) PRINTF params ; }
    
#ifdef MATLAB_MEX_FILE
#define ASSERT(expression) (mxAssert ((expression), ""))
#else
#define ASSERT(expression) (assert (expression))
#endif /* MATLAB_MEX_FILE */

PRIVATE void colamd_get_debug   /* gets the debug print level from getenv */
(   
    char *method
) ; 
    
PRIVATE void debug_deg_lists
( 
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int head [],
    Int min_score,
    Int should,
    Int max_deg
) ;

PRIVATE void debug_mark
(
    Int n_row,
    Colamd_Row Row [],
    Int tag_mark,
    Int max_mark
) ;

PRIVATE void debug_matrix
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A []
) ;

PRIVATE void debug_structures
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int in_cset [],
    Int cset_start []
) ;
/* Merged from UMFPACK */
PRIVATE void dump_super
(
    Int super_c,
    Colamd_Col Col [],
    Int n_col
) ;
/* Merged from UMFPACK */

#else /* NDEBUG */

/* === No debugging ========================================================= */

#define DEBUG0(params) ;
#define DEBUG1(params) ;
#define DEBUG2(params) ;
#define DEBUG3(params) ;
#define DEBUG4(params) ;

#define ASSERT(expression)

#endif /* NDEBUG */
