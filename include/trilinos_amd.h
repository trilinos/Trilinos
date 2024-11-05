/* ========================================================================= */
/* === AMD:  approximate minimum degree ordering =========================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD Version 2.2, Copyright (c) 2007 by Timothy A. Davis,                  */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/amd                          */
/* ------------------------------------------------------------------------- */

/* AMD finds a symmetric ordering P of a matrix A so that the Cholesky
 * factorization of P*A*P' has fewer nonzeros and takes less work than the
 * Cholesky factorization of A.  If A is not symmetric, then it performs its
 * ordering on the matrix A+A'.  Two sets of user-callable routines are
 * provided, one for int integers and the other for UF_long integers.
 *
 * The method is based on the approximate minimum degree algorithm, discussed
 * in Amestoy, Davis, and Duff, "An approximate degree ordering algorithm",
 * SIAM Journal of Matrix Analysis and Applications, vol. 17, no. 4, pp.
 * 886-905, 1996.  This package can perform both the AMD ordering (with
 * aggressive absorption), and the AMDBAR ordering (without aggressive
 * absorption) discussed in the above paper.  This package differs from the
 * Fortran codes discussed in the paper:
 *
 *	(1) it can ignore "dense" rows and columns, leading to faster run times
 *	(2) it computes the ordering of A+A' if A is not symmetric
 *	(3) it is followed by a depth-first post-ordering of the assembly tree
 *	    (or supernodal elimination tree)
 *
 * For historical reasons, the Fortran versions, amd.f and amdbar.f, have
 * been left (nearly) unchanged.  They compute the identical ordering as
 * described in the above paper.
 */

#ifndef TRILINOS_AMD_H
#define TRILINOS_AMD_H

/* make it easy for C++ programs to include AMD */
#ifdef __cplusplus
extern "C" {
#endif

/* get the definition of size_t: */
#include <stddef.h>

/* define UF_long */
#include "trilinos_UFconfig.h"

int trilinos_amd_order		    /* returns TRILINOS_AMD_OK, TRILINOS_AMD_OK_BUT_JUMBLED,
			     * TRILINOS_AMD_INVALID, or TRILINOS_AMD_OUT_OF_MEMORY */
(
    int n,		    /* A is n-by-n.  n must be >= 0. */
    const int Ap [ ],	    /* column pointers for A, of size n+1 */
    const int Ai [ ],	    /* row indices of A, of size nz = Ap [n] */
    int P [ ],		    /* output permutation, of size n */
    double Control [ ],	    /* input Control settings, of size TRILINOS_AMD_CONTROL */
    double Info [ ]	    /* output Info statistics, of size TRILINOS_AMD_INFO */
) ;

UF_long trilinos_amd_l_order	    /* see above for description of arguments */
(
    UF_long n,
    const UF_long Ap [ ],
    const UF_long Ai [ ],
    UF_long P [ ],
    double Control [ ],
    double Info [ ]
) ;

/* Input arguments (not modified):
 *
 *	n: the matrix A is n-by-n.
 *	Ap: an int/UF_long array of size n+1, containing column pointers of A.
 *	Ai: an int/UF_long array of size nz, containing the row indices of A,
 *	    where nz = Ap [n].
 *	Control:  a double array of size TRILINOS_AMD_CONTROL, containing control
 *	    parameters.  Defaults are used if Control is NULL.
 *
 * Output arguments (not defined on input):
 *
 *	P: an int/UF_long array of size n, containing the output permutation. If
 *	    row i is the kth pivot row, then P [k] = i.  In MATLAB notation,
 *	    the reordered matrix is A (P,P).
 *	Info: a double array of size TRILINOS_AMD_INFO, containing statistical
 *	    information.  Ignored if Info is NULL.
 *
 * On input, the matrix A is stored in column-oriented form.  The row indices
 * of nonzero entries in column j are stored in Ai [Ap [j] ... Ap [j+1]-1].
 *
 * If the row indices appear in ascending order in each column, and there
 * are no duplicate entries, then trilinos_amd_order is slightly more efficient in
 * terms of time and memory usage.  If this condition does not hold, a copy
 * of the matrix is created (where these conditions do hold), and the copy is
 * ordered.  This feature is new to v2.0 (v1.2 and earlier required this
 * condition to hold for the input matrix).
 * 
 * Row indices must be in the range 0 to
 * n-1.  Ap [0] must be zero, and thus nz = Ap [n] is the number of nonzeros
 * in A.  The array Ap is of size n+1, and the array Ai is of size nz = Ap [n].
 * The matrix does not need to be symmetric, and the diagonal does not need to
 * be present (if diagonal entries are present, they are ignored except for
 * the output statistic Info [TRILINOS_AMD_NZDIAG]).  The arrays Ai and Ap are not
 * modified.  This form of the Ap and Ai arrays to represent the nonzero
 * pattern of the matrix A is the same as that used internally by MATLAB.
 * If you wish to use a more flexible input structure, please see the
 * umfpack_*_triplet_to_col routines in the UMFPACK package, at
 * http://www.cise.ufl.edu/research/sparse/umfpack.
 *
 * Restrictions:  n >= 0.  Ap [0] = 0.  Ap [j] <= Ap [j+1] for all j in the
 *	range 0 to n-1.  nz = Ap [n] >= 0.  Ai [0..nz-1] must be in the range 0
 *	to n-1.  Finally, Ai, Ap, and P must not be NULL.  If any of these
 *	restrictions are not met, AMD returns TRILINOS_AMD_INVALID.
 *
 * AMD returns:
 *
 *	TRILINOS_AMD_OK if the matrix is valid and sufficient memory can be allocated to
 *	    perform the ordering.
 *
 *	TRILINOS_AMD_OUT_OF_MEMORY if not enough memory can be allocated.
 *
 *	TRILINOS_AMD_INVALID if the input arguments n, Ap, Ai are invalid, or if P is
 *	    NULL.
 *
 *	TRILINOS_AMD_OK_BUT_JUMBLED if the matrix had unsorted columns, and/or duplicate
 *	    entries, but was otherwise valid.
 *
 * The AMD routine first forms the pattern of the matrix A+A', and then
 * computes a fill-reducing ordering, P.  If P [k] = i, then row/column i of
 * the original is the kth pivotal row.  In MATLAB notation, the permuted
 * matrix is A (P,P), except that 0-based indexing is used instead of the
 * 1-based indexing in MATLAB.
 *
 * The Control array is used to set various parameters for AMD.  If a NULL
 * pointer is passed, default values are used.  The Control array is not
 * modified.
 *
 *	Control [TRILINOS_AMD_DENSE]:  controls the threshold for "dense" rows/columns.
 *	    A dense row/column in A+A' can cause AMD to spend a lot of time in
 *	    ordering the matrix.  If Control [TRILINOS_AMD_DENSE] >= 0, rows/columns
 *	    with more than Control [TRILINOS_AMD_DENSE] * sqrt (n) entries are ignored
 *	    during the ordering, and placed last in the output order.  The
 *	    default value of Control [TRILINOS_AMD_DENSE] is 10.  If negative, no
 *	    rows/columns are treated as "dense".  Rows/columns with 16 or
 *	    fewer off-diagonal entries are never considered "dense".
 *
 *	Control [TRILINOS_AMD_AGGRESSIVE]: controls whether or not to use aggressive
 *	    absorption, in which a prior element is absorbed into the current
 *	    element if is a subset of the current element, even if it is not
 *	    adjacent to the current pivot element (refer to Amestoy, Davis,
 *	    & Duff, 1996, for more details).  The default value is nonzero,
 *	    which means to perform aggressive absorption.  This nearly always
 *	    leads to a better ordering (because the approximate degrees are
 *	    more accurate) and a lower execution time.  There are cases where
 *	    it can lead to a slightly worse ordering, however.  To turn it off,
 *	    set Control [TRILINOS_AMD_AGGRESSIVE] to 0.
 *
 *	Control [2..4] are not used in the current version, but may be used in
 *	    future versions.
 *
 * The Info array provides statistics about the ordering on output.  If it is
 * not present, the statistics are not returned.  This is not an error
 * condition.
 * 
 *	Info [TRILINOS_AMD_STATUS]:  the return value of AMD, either TRILINOS_AMD_OK,
 *	    TRILINOS_AMD_OK_BUT_JUMBLED, TRILINOS_AMD_OUT_OF_MEMORY, or TRILINOS_AMD_INVALID.
 *
 *	Info [TRILINOS_AMD_N]: n, the size of the input matrix
 *
 *	Info [TRILINOS_AMD_NZ]: the number of nonzeros in A, nz = Ap [n]
 *
 *	Info [TRILINOS_AMD_SYMMETRY]:  the symmetry of the matrix A.  It is the number
 *	    of "matched" off-diagonal entries divided by the total number of
 *	    off-diagonal entries.  An entry A(i,j) is matched if A(j,i) is also
 *	    an entry, for any pair (i,j) for which i != j.  In MATLAB notation,
 *		S = spones (A) ;
 *		B = tril (S, -1) + triu (S, 1) ;
 *		symmetry = nnz (B & B') / nnz (B) ;
 *
 *	Info [TRILINOS_AMD_NZDIAG]: the number of entries on the diagonal of A.
 *
 *	Info [TRILINOS_AMD_NZ_A_PLUS_AT]:  the number of nonzeros in A+A', excluding the
 *	    diagonal.  If A is perfectly symmetric (Info [TRILINOS_AMD_SYMMETRY] = 1)
 *	    with a fully nonzero diagonal, then Info [TRILINOS_AMD_NZ_A_PLUS_AT] = nz-n
 *	    (the smallest possible value).  If A is perfectly unsymmetric
 *	    (Info [TRILINOS_AMD_SYMMETRY] = 0, for an upper triangular matrix, for
 *	    example) with no diagonal, then Info [TRILINOS_AMD_NZ_A_PLUS_AT] = 2*nz
 *	    (the largest possible value).
 *
 *	Info [TRILINOS_AMD_NDENSE]: the number of "dense" rows/columns of A+A' that were
 *	    removed from A prior to ordering.  These are placed last in the
 *	    output order P.
 *
 *	Info [TRILINOS_AMD_MEMORY]: the amount of memory used by AMD, in bytes.  In the
 *	    current version, this is 1.2 * Info  [TRILINOS_AMD_NZ_A_PLUS_AT] + 9*n
 *	    times the size of an integer.  This is at most 2.4nz + 9n.  This
 *	    excludes the size of the input arguments Ai, Ap, and P, which have
 *	    a total size of nz + 2*n + 1 integers.
 *
 *	Info [TRILINOS_AMD_NCMPA]: the number of garbage collections performed.
 *
 *	Info [TRILINOS_AMD_LNZ]: the number of nonzeros in L (excluding the diagonal).
 *	    This is a slight upper bound because mass elimination is combined
 *	    with the approximate degree update.  It is a rough upper bound if
 *	    there are many "dense" rows/columns.  The rest of the statistics,
 *	    below, are also slight or rough upper bounds, for the same reasons.
 *	    The post-ordering of the assembly tree might also not exactly
 *	    correspond to a true elimination tree postordering.
 *
 *	Info [TRILINOS_AMD_NDIV]: the number of divide operations for a subsequent LDL'
 *	    or LU factorization of the permuted matrix A (P,P).
 *
 *	Info [TRILINOS_AMD_NMULTSUBS_LDL]:  the number of multiply-subtract pairs for a
 *	    subsequent LDL' factorization of A (P,P).
 *
 *	Info [TRILINOS_AMD_NMULTSUBS_LU]:  the number of multiply-subtract pairs for a
 *	    subsequent LU factorization of A (P,P), assuming that no numerical
 *	    pivoting is required.
 *
 *	Info [TRILINOS_AMD_DMAX]:  the maximum number of nonzeros in any column of L,
 *	    including the diagonal.
 *
 *	Info [14..19] are not used in the current version, but may be used in
 *	    future versions.
 */    

/* ------------------------------------------------------------------------- */
/* direct interface to AMD */
/* ------------------------------------------------------------------------- */

/* trilinos_amd_2 is the primary AMD ordering routine.  It is not meant to be
 * user-callable because of its restrictive inputs and because it destroys
 * the user's input matrix.  It does not check its inputs for errors, either.
 * However, if you can work with these restrictions it can be faster than
 * trilinos_amd_order and use less memory (assuming that you can create your own copy
 * of the matrix for AMD to destroy).  Refer to AMD/Source/trilinos_amd_2.c for a
 * description of each parameter. */

void trilinos_amd_2
(
    int n,
    int Pe [ ],
    int Iw [ ],
    int Len [ ],
    int iwlen,
    int pfree,
    int Nv [ ],
    int Next [ ], 
    int Last [ ],
    int Head [ ],
    int Elen [ ],
    int Degree [ ],
    int W [ ],
    double Control [ ],
    double Info [ ]
) ;

void trilinos_amd_l2
(
    UF_long n,
    UF_long Pe [ ],
    UF_long Iw [ ],
    UF_long Len [ ],
    UF_long iwlen,
    UF_long pfree,
    UF_long Nv [ ],
    UF_long Next [ ], 
    UF_long Last [ ],
    UF_long Head [ ],
    UF_long Elen [ ],
    UF_long Degree [ ],
    UF_long W [ ],
    double Control [ ],
    double Info [ ]
) ;

/* ------------------------------------------------------------------------- */
/* trilinos_amd_valid */
/* ------------------------------------------------------------------------- */

/* Returns TRILINOS_AMD_OK or TRILINOS_AMD_OK_BUT_JUMBLED if the matrix is valid as input to
 * amd_order; the latter is returned if the matrix has unsorted and/or
 * duplicate row indices in one or more columns.  Returns TRILINOS_AMD_INVALID if the
 * matrix cannot be passed to trilinos_amd_order.  For trilinos_amd_order, the matrix must also
 * be square.  The first two arguments are the number of rows and the number
 * of columns of the matrix.  For its use in AMD, these must both equal n.
 *
 * NOTE: this routine returned TRUE/FALSE in v1.2 and earlier.
 */

int trilinos_amd_valid
(
    int n_row,		    /* # of rows */
    int n_col,		    /* # of columns */
    const int Ap [ ],	    /* column pointers, of size n_col+1 */
    const int Ai [ ]	    /* row indices, of size Ap [n_col] */
) ;

UF_long trilinos_amd_l_valid
(
    UF_long n_row,
    UF_long n_col,
    const UF_long Ap [ ],
    const UF_long Ai [ ]
) ;

/* ------------------------------------------------------------------------- */
/* AMD memory manager and printf routines */
/* ------------------------------------------------------------------------- */

/* The user can redefine these to change the malloc, free, and printf routines
 * that AMD uses. */

#ifndef EXTERN
#define EXTERN extern
#endif

EXTERN void *(*trilinos_amd_malloc) (size_t) ;		    /* pointer to malloc */
EXTERN void (*trilinos_amd_free) (void *) ;		    /* pointer to free */
EXTERN void *(*trilinos_amd_realloc) (void *, size_t) ;	    /* pointer to realloc */
EXTERN void *(*trilinos_amd_calloc) (size_t, size_t) ;	    /* pointer to calloc */
EXTERN int (*trilinos_amd_printf) (const char *, ...) ;	    /* pointer to printf */

/* ------------------------------------------------------------------------- */
/* AMD Control and Info arrays */
/* ------------------------------------------------------------------------- */

/* trilinos_amd_defaults:  sets the default control settings */
void trilinos_amd_defaults   (double Control [ ]) ;
void trilinos_amd_l_defaults (double Control [ ]) ;

/* trilinos_amd_control: prints the control settings */
void trilinos_amd_control    (double Control [ ]) ;
void trilinos_amd_l_control  (double Control [ ]) ;

/* trilinos_amd_info: prints the statistics */
void trilinos_amd_info       (double Info [ ]) ;
void trilinos_amd_l_info     (double Info [ ]) ;

#define TRILINOS_AMD_CONTROL 5	    /* size of Control array */
#define TRILINOS_AMD_INFO 20	    /* size of Info array */

/* contents of Control */
#define TRILINOS_AMD_DENSE 0	    /* "dense" if degree > Control [0] * sqrt (n) */
#define TRILINOS_AMD_AGGRESSIVE 1    /* do aggressive absorption if Control [1] != 0 */

/* default Control settings */
#define TRILINOS_AMD_DEFAULT_DENSE 10.0	    /* default "dense" degree 10*sqrt(n) */
#define TRILINOS_AMD_DEFAULT_AGGRESSIVE 1    /* do aggressive absorption by default */

/* contents of Info */
#define TRILINOS_AMD_STATUS 0	    /* return value of amd_order and amd_l_order */
#define TRILINOS_AMD_N 1		    /* A is n-by-n */
#define TRILINOS_AMD_NZ 2	    /* number of nonzeros in A */ 
#define TRILINOS_AMD_SYMMETRY 3	    /* symmetry of pattern (1 is sym., 0 is unsym.) */
#define TRILINOS_AMD_NZDIAG 4	    /* # of entries on diagonal */
#define TRILINOS_AMD_NZ_A_PLUS_AT 5  /* nz in A+A' */
#define TRILINOS_AMD_NDENSE 6	    /* number of "dense" rows/columns in A */
#define TRILINOS_AMD_MEMORY 7	    /* amount of memory used by AMD */
#define TRILINOS_AMD_NCMPA 8	    /* number of garbage collections in AMD */
#define TRILINOS_AMD_LNZ 9	    /* approx. nz in L, excluding the diagonal */
#define TRILINOS_AMD_NDIV 10	    /* number of fl. point divides for LU and LDL' */
#define TRILINOS_AMD_NMULTSUBS_LDL 11 /* number of fl. point (*,-) pairs for LDL' */
#define TRILINOS_AMD_NMULTSUBS_LU 12  /* number of fl. point (*,-) pairs for LU */
#define TRILINOS_AMD_DMAX 13	     /* max nz. in any column of L, incl. diagonal */

/* ------------------------------------------------------------------------- */
/* return values of AMD */
/* ------------------------------------------------------------------------- */

#define TRILINOS_AMD_OK 0		/* success */
#define TRILINOS_AMD_OUT_OF_MEMORY -1	/* malloc failed, or problem too large */
#define TRILINOS_AMD_INVALID -2		/* input arguments are not valid */
#define TRILINOS_AMD_OK_BUT_JUMBLED 1	/* input matrix is OK for amd_order, but
    * columns were not sorted, and/or duplicate entries were present.  AMD had
    * to do extra work before ordering the matrix.  This is a warning, not an
    * error.  */

/* ========================================================================== */
/* === AMD version ========================================================== */
/* ========================================================================== */

/* AMD Version 1.2 and later include the following definitions.
 * As an example, to test if the version you are using is 1.2 or later:
 *
 * #ifdef TRILINOS_AMD_VERSION
 *	if (TRILINOS_AMD_VERSION >= TRILINOS_AMD_VERSION_CODE (1,2)) ...
 * #endif
 *
 * This also works during compile-time:
 *
 *	#if defined(TRILINOS_AMD_VERSION) && (AMD_VERSION >= TRILINOS_AMD_VERSION_CODE (1,2))
 *	    printf ("This is version 1.2 or later\n") ;
 *	#else
 *	    printf ("This is an early version\n") ;
 *	#endif
 *
 * Versions 1.1 and earlier of AMD do not include a #define'd version number.
 */

#define TRILINOS_AMD_DATE "May 31, 2007"
#define TRILINOS_AMD_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define TRILINOS_AMD_MAIN_VERSION 2
#define TRILINOS_AMD_SUB_VERSION 2
#define TRILINOS_AMD_SUBSUB_VERSION 0
#define TRILINOS_AMD_VERSION TRILINOS_AMD_VERSION_CODE(TRILINOS_AMD_MAIN_VERSION,TRILINOS_AMD_SUB_VERSION)

#ifdef __cplusplus
}
#endif

#endif
