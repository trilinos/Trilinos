/* ========================================================================== */
/* === colamd/symamd - a sparse matrix column ordering algorithm ============ */
/* ========================================================================== */

/*
 * CCOLAMD version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/*
    colamd:  an approximate minimum degree column ordering algorithm,
    	for LU factorization of symmetric or unsymmetric matrices,
	QR factorization, least squares, interior point methods for
	linear programming problems, and other related problems.

    symamd:  an approximate minimum degree ordering algorithm for Cholesky
    	factorization of symmetric matrices.

    Purpose:

	Colamd computes a permutation Q such that the Cholesky factorization of
	(AQ)'(AQ) has less fill-in and requires fewer floating point operations
	than A'A.  This also provides a good ordering for sparse partial
	pivoting methods, P(AQ) = LU, where Q is computed prior to numerical
	factorization, and P is computed during numerical factorization via
	conventional partial pivoting with row interchanges.  Colamd is the
	column ordering method used in SuperLU, part of the ScaLAPACK library.
	It is also available as built-in function in MATLAB Version 6,
	available from MathWorks, Inc. (http://www.mathworks.com).  This
	routine can be used in place of colmmd in MATLAB.

    	Symamd computes a permutation P of a symmetric matrix A such that the
	Cholesky factorization of PAP' has less fill-in and requires fewer
	floating point operations than A.  Symamd constructs a matrix M such
	that M'M has the same nonzero pattern of A, and then orders the columns
	of M using colmmd.  The column ordering of M is then returned as the
	row and column ordering P of A.

    Authors:

	The authors of the code itself are Stefan I. Larimore and Timothy A.
	Davis (davis@cise.ufl.edu), University of Florida.  The algorithm was
	developed in collaboration with John Gilbert, Xerox PARC, and Esmond
	Ng, Oak Ridge National Laboratory.

    Date:

	September 8, 2003.  Version 2.3.

    Acknowledgements:

	This work was supported by the National Science Foundation, under
	grants DMS-9504974 and DMS-9803599.

    Copyright and License:

	Copyright (c) 1998-2003 by the University of Florida.
	All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

	Permission is hereby granted to use, copy, modify, and/or distribute
	this program, provided that the Copyright, this License, and the
	Availability of the original version is retained on all copies and made
	accessible to the end-user of any code or package that includes COLAMD
	or any modified version of COLAMD.

    Availability:

	The colamd/symamd library is available at

	    http://www.cise.ufl.edu/research/sparse/colamd/

	This is the http://www.cise.ufl.edu/research/sparse/colamd/colamd.c
	file.  It requires the colamd.h file.  It is required by the colamdmex.c
	and symamdmex.c files, for the MATLAB interface to colamd and symamd.

    See the ChangeLog file for changes since Version 1.0.

*/

/* ========================================================================== */
/* === Description of user-callable routines ================================ */
/* ========================================================================== */

/*
    ----------------------------------------------------------------------------
    colamd_recommended:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    int colamd_recommended (int nnz, int n_row, int n_col) ;

	    or as a C macro

	    #include "colamd.h"
	    Alen = COLAMD_RECOMMENDED (int nnz, int n_row, int n_col) ;

	Purpose:

	    Returns recommended value of Alen for use by colamd.  Returns -1
	    if any input argument is negative.  The use of this routine
	    or macro is optional.  Note that the macro uses its arguments
	    more than once, so be careful for side effects, if you pass
	    expressions as arguments to COLAMD_RECOMMENDED.  Not needed for
	    symamd, which dynamically allocates its own memory.

	Arguments (all input arguments):

	    int nnz ;		Number of nonzeros in the matrix A.  This must
				be the same value as p [n_col] in the call to
				colamd - otherwise you will get a wrong value
				of the recommended memory to use.

	    int n_row ;		Number of rows in the matrix A.

	    int n_col ;		Number of columns in the matrix A.

    ----------------------------------------------------------------------------
    colamd_set_defaults:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    colamd_set_defaults (double knobs [COLAMD_KNOBS]) ;

	Purpose:

	    Sets the default parameters.  The use of this routine is optional.

	Arguments:

	    double knobs [COLAMD_KNOBS] ;	Output only.

		Colamd: rows with more than (knobs [COLAMD_DENSE_ROW] * n_col)
		entries are removed prior to ordering.  Columns with more than
		(knobs [COLAMD_DENSE_COL] * n_row) entries are removed prior to
		ordering, and placed last in the output column ordering.

		Symamd: uses only knobs [COLAMD_DENSE_ROW], which is knobs [0].
		Rows and columns with more than (knobs [COLAMD_DENSE_ROW] * n)
		entries are removed prior to ordering, and placed last in the
		output ordering.

		COLAMD_DENSE_ROW and COLAMD_DENSE_COL are defined as 0 and 1,
		respectively, in colamd.h.  Default values of these two knobs
		are both 0.5.  Currently, only knobs [0] and knobs [1] are
		used, but future versions may use more knobs.  If so, they will
		be properly set to their defaults by the future version of
		colamd_set_defaults, so that the code that calls colamd will
		not need to change, assuming that you either use
		colamd_set_defaults, or pass a (double *) NULL pointer as the
		knobs array to colamd or symamd.

    ----------------------------------------------------------------------------
    colamd:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    int colamd (int n_row, int n_col, int Alen, int *A, int *p,
	    	double knobs [COLAMD_KNOBS], int stats [CCOLAMD_STATS]) ;

	Purpose:

	    Computes a column ordering (Q) of A such that P(AQ)=LU or
	    (AQ)'AQ=LL' have less fill-in and require fewer floating point
	    operations than factorizing the unpermuted matrix A or A'A,
	    respectively.

	Returns:

	    TRUE (1) if successful, FALSE (0) otherwise.

	Arguments:

	    int n_row ;		Input argument.

		Number of rows in the matrix A.
		Restriction:  n_row >= 0.
		Colamd returns FALSE if n_row is negative.

	    int n_col ;		Input argument.

		Number of columns in the matrix A.
		Restriction:  n_col >= 0.
		Colamd returns FALSE if n_col is negative.

	    int Alen ;		Input argument.

		Restriction (see note):
		Alen >= 2*nnz + 6*(n_col+1) + 4*(n_row+1) + n_col
		Colamd returns FALSE if these conditions are not met.

		Note:  this restriction makes an modest assumption regarding
		the size of the two typedef's structures in colamd.h.
		We do, however, guarantee that

			Alen >= colamd_recommended (nnz, n_row, n_col)

		or equivalently as a C preprocessor macro:

			Alen >= COLAMD_RECOMMENDED (nnz, n_row, n_col)

		will be sufficient.

	    int A [Alen] ;	Input argument, undefined on output.

		A is an integer array of size Alen.  Alen must be at least as
		large as the bare minimum value given above, but this is very
		low, and can result in excessive run time.  For best
		performance, we recommend that Alen be greater than or equal to
		colamd_recommended (nnz, n_row, n_col), which adds
		nnz/5 to the bare minimum value given above.

		On input, the row indices of the entries in column c of the
		matrix are held in A [(p [c]) ... (p [c+1]-1)].  The row indices
		in a given column c need not be in ascending order, and
		duplicate row indices may be be present.  However, colamd will
		work a little faster if both of these conditions are met
		(Colamd puts the matrix into this format, if it finds that the
		the conditions are not met).

		The matrix is 0-based.  That is, rows are in the range 0 to
		n_row-1, and columns are in the range 0 to n_col-1.  Colamd
		returns FALSE if any row index is out of range.

		The contents of A are modified during ordering, and are
		undefined on output.

	    int p [n_col+1] ;	Both input and output argument.

		p is an integer array of size n_col+1.  On input, it holds the
		"pointers" for the column form of the matrix A.  Column c of
		the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
		entry, p [0], must be zero, and p [c] <= p [c+1] must hold
		for all c in the range 0 to n_col-1.  The value p [n_col] is
		thus the total number of entries in the pattern of the matrix A.
		Colamd returns FALSE if these conditions are not met.

		On output, if colamd returns TRUE, the array p holds the column
		permutation (Q, for P(AQ)=LU or (AQ)'(AQ)=LL'), where p [0] is
		the first column index in the new ordering, and p [n_col-1] is
		the last.  That is, p [k] = j means that column j of A is the
		kth pivot column, in AQ, where k is in the range 0 to n_col-1
		(p [0] = j means that column j of A is the first column in AQ).

		If colamd returns FALSE, then no permutation is returned, and
		p is undefined on output.

	    double knobs [COLAMD_KNOBS] ;	Input argument.

		See colamd_set_defaults for a description.

	    int stats [CCOLAMD_STATS] ;		Output argument.

		Statistics on the ordering, and error status.
		See colamd.h for related definitions.
		Colamd returns FALSE if stats is not present.

		stats [0]:  number of dense or empty rows ignored.

		stats [1]:  number of dense or empty columns ignored (and
				ordered last in the output permutation p)
				Note that a row can become "empty" if it
				contains only "dense" and/or "empty" columns,
				and similarly a column can become "empty" if it
				only contains "dense" and/or "empty" rows.

		stats [2]:  number of garbage collections performed.
				This can be excessively high if Alen is close
				to the minimum required value.

		stats [3]:  status code.  < 0 is an error code.
			    > 1 is a warning or notice.

			0	OK.  Each column of the input matrix contained
				row indices in increasing order, with no
				duplicates.

			1	OK, but columns of input matrix were jumbled
				(unsorted columns or duplicate entries).  Colamd
				had to do some extra work to sort the matrix
				first and remove duplicate entries, but it
				still was able to return a valid permutation
				(return value of colamd was TRUE).

					stats [4]: highest numbered column that
						is unsorted or has duplicate
						entries.
					stats [5]: last seen duplicate or
						unsorted row index.
					stats [6]: number of duplicate or
						unsorted row indices.

			-1	A is a null pointer

			-2	p is a null pointer

			-3 	n_row is negative

					stats [4]: n_row

			-4	n_col is negative

					stats [4]: n_col

			-5	number of nonzeros in matrix is negative

					stats [4]: number of nonzeros, p [n_col]

			-6	p [0] is nonzero

					stats [4]: p [0]

			-7	A is too small

					stats [4]: required size
					stats [5]: actual size (Alen)

			-8	a column has a negative number of entries

					stats [4]: column with < 0 entries
					stats [5]: number of entries in col

			-9	a row index is out of bounds

					stats [4]: column with bad row index
					stats [5]: bad row index
					stats [6]: n_row, # of rows of matrx

			-10	(unused; see symamd.c)

			-999	(unused; see symamd.c)

		Future versions may return more statistics in the stats array.

	Example:

	    See http://www.cise.ufl.edu/research/sparse/colamd/example.c
	    for a complete example.

	    To order the columns of a 5-by-4 matrix with 11 nonzero entries in
	    the following nonzero pattern

	    	x 0 x 0
		x 0 x x
		0 x x 0
		0 0 x x
		x x 0 0

	    with default knobs and no output statistics, do the following:

		#include "colamd.h"
		#define ALEN COLAMD_RECOMMENDED (11, 5, 4)
		int A [ALEN] = {1, 2, 5, 3, 5, 1, 2, 3, 4, 2, 4} ;
		int p [ ] = {0, 3, 5, 9, 11} ;
		int stats [CCOLAMD_STATS] ;
		colamd (5, 4, ALEN, A, p, (double *) NULL, stats) ;

	    The permutation is returned in the array p, and A is destroyed.

    ----------------------------------------------------------------------------
    symamd:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    int symamd (int n, int *A, int *p, int *perm,
	    	double knobs [COLAMD_KNOBS], int stats [CCOLAMD_STATS],
		void (*allocate) (size_t, size_t), void (*release) (void *)) ;

	Purpose:

    	    The symamd routine computes an ordering P of a symmetric sparse
	    matrix A such that the Cholesky factorization PAP' = LL' remains
	    sparse.  It is based on a column ordering of a matrix M constructed
	    so that the nonzero pattern of M'M is the same as A.  The matrix A
	    is assumed to be symmetric; only the strictly lower triangular part
	    is accessed.  You must pass your selected memory allocator (usually
	    calloc/free or mxCalloc/mxFree) to symamd, for it to allocate
	    memory for the temporary matrix M.

	    TODO: uses lower, upper, or both.  Reword

	Returns:

	    TRUE (1) if successful, FALSE (0) otherwise.

	Arguments:

	    int n ;		Input argument.

	    	Number of rows and columns in the symmetrix matrix A.
		Restriction:  n >= 0.
		Symamd returns FALSE if n is negative.

	    int A [nnz] ;	Input argument.

	    	A is an integer array of size nnz, where nnz = p [n].

		The row indices of the entries in column c of the matrix are
		held in A [(p [c]) ... (p [c+1]-1)].  The row indices in a
		given column c need not be in ascending order, and duplicate
		row indices may be present.  However, symamd will run faster
		if the columns are in sorted order with no duplicate entries.

		The matrix is 0-based.  That is, rows are in the range 0 to
		n-1, and columns are in the range 0 to n-1.  Symamd
		returns FALSE if any row index is out of range.

		The contents of A are not modified.

	    int p [n+1] ;   	Input argument.

		p is an integer array of size n+1.  On input, it holds the
		"pointers" for the column form of the matrix A.  Column c of
		the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
		entry, p [0], must be zero, and p [c] <= p [c+1] must hold
		for all c in the range 0 to n-1.  The value p [n] is
		thus the total number of entries in the pattern of the matrix A.
		Symamd returns FALSE if these conditions are not met.

		The contents of p are not modified.

	    int perm [n+1] ;   	Output argument.

		On output, if symamd returns TRUE, the array perm holds the
		permutation P, where perm [0] is the first index in the new
		ordering, and perm [n-1] is the last.  That is, perm [k] = j
		means that row and column j of A is the kth column in PAP',
		where k is in the range 0 to n-1 (perm [0] = j means
		that row and column j of A are the first row and column in
		PAP').  The array is used as a workspace during the ordering,
		which is why it must be of length n+1, not just n.

	    double knobs [COLAMD_KNOBS] ;	Input argument.

		See colamd_set_defaults for a description.

	    int stats [CCOLAMD_STATS] ;		Output argument.

		Statistics on the ordering, and error status.
		See colamd.h for related definitions.
		Symamd returns FALSE if stats is not present.

		stats [0]:  number of dense or empty row and columns ignored
				(and ordered last in the output permutation
				perm).  Note that a row/column can become
				"empty" if it contains only "dense" and/or
				"empty" columns/rows.

		stats [1]:  (same as stats [0])

		stats [2]:  number of garbage collections performed.

		stats [3]:  status code.  < 0 is an error code.
			    > 1 is a warning or notice.

			0	OK.  Each column of the input matrix contained
				row indices in increasing order, with no
				duplicates.

			1	OK, but columns of input matrix were jumbled
				(unsorted columns or duplicate entries).  Symamd
				had to do some extra work to sort the matrix
				first and remove duplicate entries, but it
				still was able to return a valid permutation
				(return value of symamd was TRUE).

					stats [4]: highest numbered column that
						is unsorted or has duplicate
						entries.
					stats [5]: last seen duplicate or
						unsorted row index.
					stats [6]: number of duplicate or
						unsorted row indices.

			-1	A is a null pointer

			-2	p is a null pointer

			-3	(unused, see colamd.c)

			-4 	n is negative

					stats [4]: n

			-5	number of nonzeros in matrix is negative

					stats [4]: # of nonzeros (p [n]).

			-6	p [0] is nonzero

					stats [4]: p [0]

			-7	(unused)

			-8	a column has a negative number of entries

					stats [4]: column with < 0 entries
					stats [5]: number of entries in col

			-9	a row index is out of bounds

					stats [4]: column with bad row index
					stats [5]: bad row index
					stats [6]: n_row, # of rows of matrx

			-10	out of memory (unable to allocate temporary
				workspace for M or count arrays using the
				"allocate" routine passed into symamd).

		Future versions may return more statistics in the stats array.

	    void * (*allocate) (size_t, size_t)

	    	A pointer to a function providing memory allocation.  The
		allocated memory must be returned initialized to zero.  For a
		C application, this argument should normally be a pointer to
		calloc.  For a MATLAB mexFunction, the routine mxCalloc is
		passed instead.

	    void (*release) (size_t, size_t)

	    	A pointer to a function that frees memory allocated by the
		memory allocation routine above.  For a C application, this
		argument should normally be a pointer to free.  For a MATLAB
		mexFunction, the routine mxFree is passed instead.


    ----------------------------------------------------------------------------
    colamd_report:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    colamd_report (int stats [CCOLAMD_STATS]) ;

	Purpose:

	    Prints the error status and statistics recorded in the stats
	    array on the standard error output (for a standard C routine)
	    or on the MATLAB output (for a mexFunction).

	Arguments:

	    int stats [CCOLAMD_STATS] ;	Input only.  Statistics from colamd.


    ----------------------------------------------------------------------------
    symamd_report:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    symamd_report (int stats [CCOLAMD_STATS]) ;

	Purpose:

	    Prints the error status and statistics recorded in the stats
	    array on the standard error output (for a standard C routine)
	    or on the MATLAB output (for a mexFunction).

	Arguments:

	    int stats [CCOLAMD_STATS] ;	Input only.  Statistics from symamd.


*/

/* ========================================================================== */
/* === Scaffolding code definitions  ======================================== */
/* ========================================================================== */

/* ===== DEBUG definition moved ccolamd.h for ccolamd ======================= */

/*
   Our "scaffolding code" philosophy:  In our opinion, well-written library
   code should keep its "debugging" code, and just normally have it turned off
   by the compiler so as not to interfere with performance.  This serves
   several purposes:

   (1) assertions act as comments to the reader, telling you what the code
	expects at that point.  All assertions will always be true (unless
	there really is a bug, of course).

   (2) leaving in the scaffolding code assists anyone who would like to modify
	the code, or understand the algorithm (by reading the debugging output,
	one can get a glimpse into what the code is doing).

   (3) (gasp!) for actually finding bugs.  This code has been heavily tested
	and "should" be fully functional and bug-free ... but you never know...

    To enable debugging, comment out the "#define NDEBUG" above.  For a MATLAB
    mexFunction, you will also need to modify mexopts.sh to remove the -DNDEBUG
    definition.  The code will become outrageously slow when debugging is
    enabled.  To control the level of debugging output, set an environment
    variable D to 0 (little), 1 (some), 2, 3, or 4 (lots).  When debugging,
    you should see the following message on the standard output:

    	colamd: debug version, D = 1 (THIS WILL BE SLOW!)

    or a similar message for symamd.  If you don't, then debugging has not
    been enabled.

*/

/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include "ccolamd.h"
#include "ccolamd_internal.h"
#include <limits.h>

/* ========================================================================== */
/* === Definitions ========================================================== */
/* ========================================================================== */

/* === PUBLIC and PRIVATE moved to ccolamd.h for ccolamd ==================== */

/* === MAX & MIN moved to ccolamd.h for ccolamd ============================= */

#define ONES_COMPLEMENT(r) (-(r)-1)

/* -------------------------------------------------------------------------- */
/* Change for version 2.1:  define TRUE and FALSE only if not yet defined */
/* -------------------------------------------------------------------------- */

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

/* -------------------------------------------------------------------------- */

/* === EMPTY moved to ccolamd.h for ccolamd ================================= */

/* Row and column status */
#define ALIVE	(0)
#define DEAD	(-1)

/* Column status */
#define DEAD_PRINCIPAL		(-1)
#define DEAD_NON_PRINCIPAL	(-2)

/* Macros for row and column status update and checking. */
#define ROW_IS_DEAD(r)			ROW_IS_MARKED_DEAD (Row[r].shared2.mark)
#define ROW_IS_MARKED_DEAD(row_mark)	(row_mark < ALIVE)
#define ROW_IS_ALIVE(r)			(Row [r].shared2.mark >= ALIVE)
#define COL_IS_DEAD(c)			(Col [c].start < ALIVE)
#define COL_IS_ALIVE(c)			(Col [c].start >= ALIVE)
#define COL_IS_DEAD_PRINCIPAL(c)	(Col [c].start == DEAD_PRINCIPAL)
#define KILL_ROW(r)			{ Row [r].shared2.mark = DEAD ; }
#define KILL_PRINCIPAL_COL(c)		{ Col [c].start = DEAD_PRINCIPAL ; }
#define KILL_NON_PRINCIPAL_COL(c)	{ Col [c].start = DEAD_NON_PRINCIPAL ; }

/* ========================================================================== */
/* === Colamd reporting mechanism moved to ccolamd.h for ccolamd ============ */
/* ========================================================================== */

/* ========================================================================== */
/* === Prototypes of PRIVATE routines ======================================= */
/* ========================================================================== */

PRIVATE Int init_rows_cols
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int p [],
    Int stats [CCOLAMD_STATS]
) ;

PRIVATE void init_scoring
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    double knobs [CCOLAMD_KNOBS],
    Int *p_n_row2,
    Int *p_n_col2,
    Int *p_max_deg
    /* ------------------ */
    /* Added for ccolamd */
    , Int in_cset []
    , Int n_cset
    , Int cset_start []
    , Int dead_cols []
    /* added for ccolamd from UMFPACK */
    , Int *p_ndense_row         /* number of dense rows */
    , Int *p_nempty_row         /* number of original empty rows */
    , Int *p_nnewlyempty_row    /* number of newly empty rows */
    , Int *p_ndense_col         /* number of dense cols (excl "empty" cols) */
    , Int *p_nempty_col         /* number of original empty cols */
    , Int *p_nnewlyempty_col    /* number of newly empty cols */
    /* ------------------ */
    , Int max_col_degree
) ;

PRIVATE Int find_ordering
(
    Int n_row,
    Int n_col,
    Int Alen,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
#ifndef NDEBUG
    Int n_col2,
#endif
    Int max_deg,
    Int pfree
    /* ------------------ */
    /*[ Added for ccolamd */
    , Int cset []
    , Int cset_start []
#ifndef NDEBUG
    , Int n_cset
#endif
    , Int in_cset []

    /* [ Merged from UMFPACK */
    , Int Front_npivcol []
    , Int Front_nrows []
    , Int Front_ncols []
    , Int Front_parent []
    , Int Front_cols []
    , Int *p_nfr
    , Int aggressive
    , Int InFront []
    , Int fact_type
    /* Merged from UMFPACK ] */
    /* ------------------  ] */
) ;


/* ============== order_children deleted for ccolamd =========================*/

PRIVATE void detect_super_cols
(

#ifndef NDEBUG
    Int n_col,
    Colamd_Row Row [],
#endif /* NDEBUG */

    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int row_start,
    Int row_length

    /* ------------------ */
    /* Added for ccolamd */
    , Int in_set []
    /* ------------------ */
) ;

PRIVATE Int garbage_collection
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int *pfree
) ;

PRIVATE Int clear_mark
(
    Int tag_mark,
    Int max_mark,
    Int n_row,
    Colamd_Row Row []
) ;

PRIVATE void print_report
(
    char *method,
    Int stats [CCOLAMD_STATS]
) ;

/* ========================================================================== */
/* === Debugging prototypes and definitions moved to ccolamd.h for ccolamd == */
/* ========================================================================== */

#ifndef NDEBUG

/* colamd_debug is the *ONLY* global variable, and is only */
/* present when debugging */

/* ------------------- */
/* changed for ccolamd */
Int colamd_debug ;      /* debug print level */
/* ------------------- */

#endif

/* ========================================================================== */
/* === USER-CALLABLE ROUTINES: ============================================== */
/* ========================================================================== */


/* ========================================================================== */
/* === colamd_recommended =================================================== */
/* ========================================================================== */

/*
    The colamd_recommended routine returns the suggested size for Alen.  This
    value has been determined to provide good balance between the number of
    garbage collections and the memory requirements for colamd.  If any
    argument is negative, a -1 is returned as an error condition.  This
    function is also available as a macro defined in colamd.h, so that you
    can use it for a statically-allocated array size.
*/

PUBLIC Int CCOLAMD_f_recommended	/* returns recommended value of Alen. */
(
    /* === Parameters ====================================================== */

    Int nnz,			/* number of nonzeros in A */
    Int n_row,			/* number of rows in A */
    Int n_col			/* number of columns in A */
)
{
    return (CCOLAMD_RECOMMENDED (nnz, n_row, n_col)) ;
}


/* ========================================================================== */
/* === colamd_set_defaults ================================================== */
/* ========================================================================== */

/*
    The colamd_set_defaults routine sets the default values of the user-
    controllable parameters for colamd:

	knobs [0]	rows with knobs[0]*n_col entries or more are removed
			prior to ordering in colamd.  Rows and columns with
			knobs[0]*n_col entries or more are removed prior to
			ordering in symamd and placed last in the output
			ordering.

	knobs [1]	columns with knobs[1]*n_row entries or more are removed
			prior to ordering in colamd, and placed last in the
			column permutation.  Symamd ignores this knob.

	knobs [2..19]	unused, but future versions might use this
*/

PUBLIC void CCOLAMD_set_defaults
(
    /* === Parameters ======================================================= */

    double knobs [CCOLAMD_KNOBS]		/* knob array */
)
{
    /* === Local variables ================================================== */

    Int i ;

    if (!knobs)
    {
	return ;			/* no knobs to initialize */
    }
    for (i = 0 ; i < CCOLAMD_KNOBS ; i++)
    {
	knobs [i] = 0 ;
    }
    /* TODO: determine proper defaults */

#ifdef COLAMD_DENSE_KNOBS
    /* COLAMD */
    knobs [CCOLAMD_DENSE_ROW] = 0.5 ;	/* colamd default */
    knobs [CCOLAMD_DENSE_COL] = 0.5 ;	/* colamd default */
#else
    /* UMFPACK is 0.2 */
    knobs [CCOLAMD_DENSE_ROW] = 0.625 ;	/* to match AMD */
    knobs [CCOLAMD_DENSE_COL] = 0.625 ;	/* to match AMD */
#endif
    /* ----------------- */
    /* Added for ccolamd */
    knobs [CCOLAMD_AGGRESSIVE] = TRUE ;	/* default is to do aggressive */
    					/* absorption*/
    knobs [CCOLAMD_FACT_TYPE] = 1 ;	/* default is to do LU */
    /* ----------------- */
}


/* ========================================================================== */
/* === symamd =============================================================== */
/* ========================================================================== */

PUBLIC Int CSYMAMD_MAIN		/* return TRUE if OK, FALSE otherwise */
(
    /* === Parameters ======================================================= */

    Int n,				/* number of rows and columns of A */
    Int A [],				/* row indices of A */
    Int p [],				/* column pointers of A */
    Int perm [],			/* output permutation, size n+1 */
    double knobs [CCOLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    Int stats [CCOLAMD_STATS],		/* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
    /* ------------------ */
    /* Added for ccolamd */
    , Int cset [], Int symmetry
    /* ------------------ */
)
{
    /* === Local variables ================================================== */

    Int *count ;		/* length of each column of M, and col pointer*/
    Int *mark ;			/* mark array for finding duplicate entries */
    Int *M ;			/* row indices of matrix M */
    Int Mlen ;			/* length of M */
    Int n_row ;			/* number of rows in M */
    Int nnz ;			/* number of entries in A */
    Int i ;			/* row index of A */
    Int j ;			/* column index of A */
    Int k ;			/* row index of M */
    Int mnz ;			/* number of nonzeros in M */
    Int pp ;			/* index into a column of A */
    Int last_row ;		/* last row seen in the current column */
    Int length ;		/* number of nonzeros in a column */

    Int both ;			/* TRUE if ordering A+A' */
    Int upper ;			/* TRUE if ordering triu(A)+triu(A)' */
    Int lower ;			/* TRUE if ordering tril(A)+tril(A)' */

    double cknobs [CCOLAMD_KNOBS] ;		/* knobs for colamd */
    double default_knobs [CCOLAMD_KNOBS] ;	/* default knobs for colamd */
    Int cstats [CCOLAMD_STATS] ;			/* colamd stats */

#ifndef NDEBUG
    colamd_get_debug ("csymamd") ;
#endif /* NDEBUG */

    both = (symmetry == 0) ;
    upper = (symmetry > 0) ;
    lower = (symmetry < 0) ;

    /* === Check the input arguments ======================================== */

    if (!stats)
    {
	DEBUG0 (("csymamd: stats not present\n")) ;
	return (FALSE) ;
    }
    for (i = 0 ; i < CCOLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [CCOLAMD_STATUS] = CCOLAMD_OK ;
    stats [CCOLAMD_INFO1] = -1 ;
    stats [CCOLAMD_INFO2] = -1 ;

    if (!A)
    {
    	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_A_not_present ;
	DEBUG0 (("csymamd: A not present\n")) ;
	return (FALSE) ;
    }

    if (!p)		/* p is not present */
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_p_not_present ;
	DEBUG0 (("csymamd: p not present\n")) ;
    	return (FALSE) ;
    }

    if (n < 0)		/* n must be >= 0 */
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_ncol_negative ;
	stats [CCOLAMD_INFO1] = n ;
	DEBUG0 (("csymamd: n negative "ID" \n", n)) ;
    	return (FALSE) ;
    }

    nnz = p [n] ;
    if (nnz < 0)	/* nnz must be >= 0 */
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_nnz_negative ;
	stats [CCOLAMD_INFO1] = nnz ;
	DEBUG0 (("csymamd: number of entries negative "ID" \n", nnz)) ;
	return (FALSE) ;
    }

    if (p [0] != 0)
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_p0_nonzero ;
	stats [CCOLAMD_INFO1] = p [0] ;
	DEBUG0 (("csymamd: p[0] not zero "ID"\n", p [0])) ;
	return (FALSE) ;
    }

    /* === If no knobs, set default knobs =================================== */

    if (!knobs)
    {
	ccolamd_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }

    /* === Allocate count and mark ========================================== */

    count = (Int *) ((*allocate) (n+1, sizeof (Int))) ;

    if (!count)
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_out_of_memory ;
	DEBUG0 (("csymamd: allocate count (size "ID") failed\n", n+1)) ;
	return (FALSE) ;
    }

    mark = (Int *) ((*allocate) (n+1, sizeof (Int))) ;

    if (!mark)
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_out_of_memory ;
	(*release) ((void *) count) ;
	DEBUG0 (("csymamd: allocate mark (size "ID") failed\n", n+1)) ;
	return (FALSE) ;
    }

    /* === Compute column counts of M, check if A is valid ================== */

    stats [CCOLAMD_INFO3] = 0 ; /* number of duplicate or unsorted row indices*/

    for (i = 0 ; i < n ; i++)
    {
    	mark [i] = -1 ;
    }

    for (j = 0 ; j < n ; j++)
    {
	last_row = -1 ;

	length = p [j+1] - p [j] ;
	if (length < 0)
	{
	    /* column pointers must be non-decreasing */
	    stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_col_length_negative ;
	    stats [CCOLAMD_INFO1] = j ;
	    stats [CCOLAMD_INFO2] = length ;
	    (*release) ((void *) count) ;
	    (*release) ((void *) mark) ;
	    DEBUG0 (("csymamd: col "ID" negative length "ID"\n", j, length)) ;
	    return (FALSE) ;
	}

	for (pp = p [j] ; pp < p [j+1] ; pp++)
	{
	    i = A [pp] ;
	    if (i < 0 || i >= n)
	    {
		/* row index i, in column j, is out of bounds */
		stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_row_index_out_of_bounds ;
		stats [CCOLAMD_INFO1] = j ;
		stats [CCOLAMD_INFO2] = i ;
		stats [CCOLAMD_INFO3] = n ;
		(*release) ((void *) count) ;
		(*release) ((void *) mark) ;
		DEBUG0 (("csymamd: row "ID" col "ID" out of bounds\n", i, j)) ;
		return (FALSE) ;
	    }

	    if (i <= last_row || mark [i] == j)
	    {
		/* row index is unsorted or repeated (or both), thus col */
		/* is jumbled.  This is a notice, not an error condition. */
		stats [CCOLAMD_STATUS] = CCOLAMD_OK_BUT_JUMBLED ;
		stats [CCOLAMD_INFO1] = j ;
		stats [CCOLAMD_INFO2] = i ;
		(stats [CCOLAMD_INFO3]) ++ ;
		DEBUG1 (("csymamd: row "ID" col "ID" unsorted/dupl.\n", i, j)) ;
	    }

	    if (mark [i] != j)
	    {
		if ((both && i != j) || (lower && i > j) || (upper && i < j))
		{
		    /* row k of M will contain column indices i and j */
		    count [i]++ ;
		    count [j]++ ;
		}
	    }

	    /* mark the row as having been seen in this column */
	    mark [i] = j ;

	    last_row = i ;
	}
    }

    /* === Compute column pointers of M ===================================== */

    /* use output permutation, perm, for column pointers of M */
    perm [0] = 0 ;
    for (j = 1 ; j <= n ; j++)
    {
	perm [j] = perm [j-1] + count [j-1] ;
    }
    for (j = 0 ; j < n ; j++)
    {
	count [j] = perm [j] ;
    }

    /* === Construct M ====================================================== */

    mnz = perm [n] ;
    n_row = mnz / 2 ;
    Mlen = ccolamd_recommended (mnz, n_row, n) ;
    M = (Int *) ((*allocate) (Mlen, sizeof (Int))) ;
    DEBUG0 (("csymamd: M is "ID"-by-"ID" with "ID" entries, Mlen = "ID"\n",
    	n_row, n, mnz, Mlen)) ;

    if (!M)
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_out_of_memory ;
	(*release) ((void *) count) ;
	(*release) ((void *) mark) ;
	DEBUG0 (("csymamd: allocate M (size "ID") failed\n", Mlen)) ;
	return (FALSE) ;
    }

    k = 0 ;

    if (stats [CCOLAMD_STATUS] == CCOLAMD_OK)
    {
	/* Matrix is OK */
	for (j = 0 ; j < n ; j++)
	{
	    ASSERT (p [j+1] - p [j] >= 0) ;
	    for (pp = p [j] ; pp < p [j+1] ; pp++)
	    {
		i = A [pp] ;
		ASSERT (i >= 0 && i < n) ;
		if ((both && i != j) || (lower && i > j) || (upper && i < j))
		{
		    /* row k of M contains column indices i and j */
		    M [count [i]++] = k ;
		    M [count [j]++] = k ;
		    k++ ;
		}
	    }
	}
    }
    else
    {
	/* Matrix is jumbled.  Do not add duplicates to M.  Unsorted cols OK. */
	DEBUG0 (("csymamd: Duplicates in A.\n")) ;
	for (i = 0 ; i < n ; i++)
	{
	    mark [i] = -1 ;
	}
	for (j = 0 ; j < n ; j++)
	{
	    ASSERT (p [j+1] - p [j] >= 0) ;
	    for (pp = p [j] ; pp < p [j+1] ; pp++)
	    {
		i = A [pp] ;
		ASSERT (i >= 0 && i < n) ;
		if (mark [i] != j)
		{
		    if ((both && i != j) || (lower && i > j) || (upper && i < j))
		    {
			/* row k of M contains column indices i and j */
			M [count [i]++] = k ;
			M [count [j]++] = k ;
			k++ ;
			mark [i] = j ;
		    }
		}
	    }
	}
    }

    /* count and mark no longer needed */
    (*release) ((void *) mark) ;
    (*release) ((void *) count) ;
    ASSERT (k == n_row) ;

    /* === Adjust the knobs for M =========================================== */

    for (i = 0 ; i < CCOLAMD_KNOBS ; i++)
    {
	cknobs [i] = knobs [i] ;
    }

    /* there are no dense rows in M */
    cknobs [CCOLAMD_DENSE_ROW] = 1.0 ;

    if (n_row != 0 && n < n_row)
    {
#ifdef COLAMD_DENSE_KNOBS
	/* On input, the knob is a fraction of 1..n, the number of rows of A. */
	/* Convert it to a fraction of 1..n_row, of the number of rows of M. */
    	cknobs [CCOLAMD_DENSE_COL] = (knobs [CCOLAMD_DENSE_ROW] * n) / n_row ;
#else
	/* no change to cknobs */
#endif
    }
    else
    {
	/* no dense columns in M */
    	cknobs [CCOLAMD_DENSE_COL] = 1.0 ;
    }

    DEBUG0 (("csymamd: dense col knob for M: %g\n", cknobs [CCOLAMD_DENSE_COL]));

    /* === Order the columns of M =========================================== */

    /* ------------------ */
    /* Added cset as parameter for ccolamd */

    (void) CCOLAMD_umf_ccolamd_MAIN (n_row, n, Mlen, M, perm, cknobs, cstats,
             (Int *) NULL, (Int *) NULL, (Int *) NULL, (Int *) NULL,
             (Int *) NULL, (Int *) NULL, (Int *) NULL
             , cset, n) ;

    /* Note that the output permutation is now in perm */

    /* === get the statistics for symamd from colamd ======================== */

    /* note that a dense column in colamd means a dense row and col in symamd */
    stats [CCOLAMD_DENSE_ROW]    = cstats [CCOLAMD_DENSE_COL] ;
    stats [CCOLAMD_DENSE_COL]    = cstats [CCOLAMD_DENSE_COL] ;
    stats [CCOLAMD_DEFRAG_COUNT] = cstats [CCOLAMD_DEFRAG_COUNT] ;

    /* === Free M =========================================================== */

    (*release) ((void *) M) ;
    DEBUG0 (("csymamd: done.\n")) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === ccolamd ============================================================== */
/* ========================================================================== */

PUBLIC Int CCOLAMD_MAIN
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows in A */
    Int n_col,			/* number of columns in A */
    Int Alen,			/* length of A */
    Int A [],			/* row indices of A */
    Int p [],			/* pointers to columns in A */
    double knobs [CCOLAMD_KNOBS],/* parameters (uses defaults if NULL) */
    Int stats [CCOLAMD_STATS]	/* output statistics and error codes */
    /* ------------------ */
    /* Added for ccolamd */
    , Int in_cset []		/* constraint set of A */
    /* ------------------ */
)
{
    /* The Front arrays will be constructed in ccolamd within A whereas in
     * UMFPACK they will be passed to the ccolamd.
     */

     return (CCOLAMD_umf_ccolamd_MAIN (n_row, n_col, Alen, A, p, knobs, stats,
             (Int *) NULL, (Int *) NULL, (Int *) NULL, (Int *) NULL,
             (Int *) NULL, (Int *) NULL, (Int *) NULL
             , in_cset, n_row)) ;
}

/* ========================================================================== */
/* === umf_ccolamd ========================================================== */
/* ========================================================================== */

/*
    The colamd routine computes a column ordering Q of a sparse matrix
    A such that the LU factorization P(AQ) = LU remains sparse, where P is
    selected via partial pivoting.   The routine can also be viewed as
    providing a permutation Q such that the Cholesky factorization
    (AQ)'(AQ) = LL' remains sparse.
*/

PUBLIC Int CCOLAMD_umf_ccolamd_MAIN	/* returns TRUE if successful, FALSE otherwise*/
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows in A */
    Int n_col,			/* number of columns in A */
    Int Alen,			/* length of A */
    Int A [],			/* row indices of A */
    Int p [],			/* pointers to columns in A */
    double knobs [CCOLAMD_KNOBS],/* parameters (uses defaults if NULL) */
    Int stats [CCOLAMD_STATS]	/* output statistics and error codes */

    /* ------------------------------ */
    /* Added from UMFPACK for ccolamd */
    /* each Front array is of size n_col+1. */
    , Int Front_npivcol [ ]     /* # pivot cols in each front */
    , Int Front_nrows [ ]       /* # of rows in each front (incl. pivot rows) */
    , Int Front_ncols [ ]       /* # of cols in each front (incl. pivot cols) */
    , Int Front_parent [ ]      /* parent of each front */
    , Int Front_cols [ ]        /* link list of pivot columns for each front */
    , Int *p_nfr                /* total number of frontal matrices */
    , Int InFront [ ]           /* InFront [row] = f if the original row was
                                 * absorbed into front f.  EMPTY if the row was
				 * empty, dense, or not absorbed.  This array
				 * has size n_row+1 */
    /* ------------------------------ */

    /* ------------------ */
    /* Added for ccolamd */
    , Int in_cset []		/* constraint set of A */
    , Int max_col_degree	/* n_row for ccolamd, n for csymamd */
    /* ------------------ */
)
{
    /* === Local variables ================================================== */

    Int i ;			/* loop index */
    Int nnz ;			/* nonzeros in A */
    Int Row_size ;		/* size of Row [], in integers */
    Int Col_size ;		/* size of Col [], in integers */
    Int need ;			/* minimum required length of A */
    Colamd_Row *Row ;		/* pointer into A of Row [0..n_row] array */
    Colamd_Col *Col ;		/* pointer into A of Col [0..n_col] array */
    Int n_col2 ;		/* number of non-dense, non-empty columns */
    Int n_row2 ;		/* number of non-dense, non-empty rows */
    Int ngarbage ;		/* number of garbage collections performed */
    Int max_deg ;		/* maximum row degree */
    double default_knobs [CCOLAMD_KNOBS] ;	/* default knobs array */

    /* ------------------ */
    /* [ Added for ccolamd */

    Int n_cset ;		/* number of constraint sets */
    Int *cset ;			/* cset of A */
    Int *cset_start ;		/* pointer into cset */
    Int *temp_cstart ;		/* temp pointer to start of cset */
    Int *csize ;		/* temp pointer to cset size */
    Int ap ;			/* column index */
    Int fact_type ;		/* lu or cholesky ? */

    Int ndense_row, nempty_row, parent, ndense_col,
    	nempty_col, k, col, nfr, *Front_child, *Front_sibling, *Front_stack,
    	*Front_order, *Front_size ;
    Int nnewlyempty_col, nnewlyempty_row ;
    Int aggressive ;
    Int row ;
    Int *dead_cols ;
    Int set1 ;
    Int set2 ;
    Int cs ;

    /* ------------------ ] */


#ifndef NDEBUG
    colamd_get_debug ("ccolamd") ;
#endif /* NDEBUG */

    /* === Check the input arguments ======================================== */

    if (!stats)
    {
	DEBUG0 (("ccolamd: stats not present\n")) ;
	return (FALSE) ;
    }
    for (i = 0 ; i < CCOLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [CCOLAMD_STATUS] = CCOLAMD_OK ;
    stats [CCOLAMD_INFO1] = -1 ;
    stats [CCOLAMD_INFO2] = -1 ;

    if (!A)		/* A is not present */
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_A_not_present ;
	DEBUG0 (("ccolamd: A not present\n")) ;
	return (FALSE) ;
    }

    if (!p)		/* p is not present */
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_p_not_present ;
	DEBUG0 (("ccolamd: p not present\n")) ;
    	return (FALSE) ;
    }

    if (n_row < 0)	/* n_row must be >= 0 */
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_nrow_negative ;
	stats [CCOLAMD_INFO1] = n_row ;
	DEBUG0 (("ccolamd: nrow negative "ID"\n", n_row)) ;
    	return (FALSE) ;
    }

    if (n_col < 0)	/* n_col must be >= 0 */
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_ncol_negative ;
	stats [CCOLAMD_INFO1] = n_col ;
	DEBUG0 (("ccolamd: ncol negative "ID"\n", n_col)) ;
    	return (FALSE) ;
    }

    nnz = p [n_col] ;
    if (nnz < 0)	/* nnz must be >= 0 */
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_nnz_negative ;
	stats [CCOLAMD_INFO1] = nnz ;
	DEBUG0 (("ccolamd: number of entries negative "ID"\n", nnz)) ;
	return (FALSE) ;
    }

    if (p [0] != 0)
    {
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_p0_nonzero	;
	stats [CCOLAMD_INFO1] = p [0] ;
	DEBUG0 (("ccolamd: p[0] not zero "ID"\n", p [0])) ;
	return (FALSE) ;
    }

    /* === If no knobs, set default knobs =================================== */

    if (!knobs)
    {
	ccolamd_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }

    /* ------------------------------- */
    /* added for ccolamd from UMFPACK: */
    aggressive = (knobs [CCOLAMD_AGGRESSIVE] != 0) ;
    fact_type = knobs [CCOLAMD_FACT_TYPE] ;
    /* ------------------------------- */

    /* === Allocate the Row and Col arrays from array A ===================== */

    Col_size = CCOLAMD_C (n_col) ;
    Row_size = CCOLAMD_R (n_row) ;
    /* ---------------------- */
    /* [ Modified for ccolamd */
    /* min size of A is 2nnz+ncol.  cset and cset_start are of size 2ncol+1 */
    /* Each of the 5 fronts is of size n_col + 1. InFront is of size nrow.  */
    need = MAX(2*nnz, 4*n_col) + n_col +
    		Col_size + Row_size +
		(3*n_col+1) + (5*(n_col+1)) + n_row ;
    /*  Modified for ccolamd ] */
    /* ----------------------- */

    if (need > Alen)
    {
	/* not enough space in array A to perform the ordering */
	stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_A_too_small ;
	stats [CCOLAMD_INFO1] = need ;
	stats [CCOLAMD_INFO2] = Alen ;
	DEBUG0 (("ccolamd: Need Alen >= "ID", given only Alen = "ID"\n", need,Alen));
	return (FALSE) ;
    }

    /* =================[ Added for ccolamd ================================= */

    Alen -= Col_size + Row_size + 3*n_col + 1 + 5*(n_col+1) + n_row ;

    /* Size of A is now Alen >= MAX (2*nnz, 4*n_col) + n_col.  The ordering
     * requires Alen >= 2*nnz + n_col, and the postorder requires
     * Alen >= 5*n_col. */

    ap = Alen ;

    if ( !Front_npivcol || !Front_nrows || !Front_ncols || !Front_parent ||
         !Front_cols || !Front_cols || !InFront )
    {
	Front_npivcol = &A [ap] ;
	ap += (n_col + 1) ;
	Front_nrows = &A [ap] ;
	ap += (n_col + 1) ;
	Front_ncols = &A [ap] ;
	ap += (n_col + 1) ;
	Front_parent = &A [ap] ;
	ap += (n_col + 1) ;
	Front_cols = &A [ap] ;
	ap += (n_col + 1) ;
	InFront = &A [ap] ;
	ap += (n_row) ;
    }
    else
    {
	/* Fronts are present. Leave the additional space as elbow room. */
    	ap += 5*(n_col+1) + n_row ;
	ap = Alen ;
    }

    /* cset_start is of size n_col + 1 */
    cset_start = &A [ap] ;
    ap += n_col+ 1 ;

    /* dead_col is of size n_col */
    dead_cols = &A [ap] ;
    ap += n_col ;

    /* cset is of size n_col */
    cset = &A [ap] ;
    ap += n_col ;

    /* Col is of size Col_size.  The space is shared by temp_cstart and csize */
    Col = (Colamd_Col *) &A [ap] ;
    temp_cstart = (Int *) Col ;		/* [ temp_cstart is of size n_col+1 */
    csize = temp_cstart + (n_col+1) ;	/* csize is of size n_col */
    ap += Col_size ;
    ASSERT (Col_size >= 2*n_col+1) ;

    /* Row is of size Row_size */
    Row = (Colamd_Row *) &A [ap] ;
    ap += Row_size ;

    /* Initialize csize & dead_cols to zero */
    for ( i = 0 ; i < n_col ; i++ )
    {
    	csize [i] = 0 ;
	dead_cols [i] = 0 ;
    }

    /* === Construct the constraint set ===================================== */

    if (in_cset == (int *) NULL)
    {
	/* no constraint set; all columns belong to set zero */
	n_cset = 1 ;
	csize [0] = n_col ;
	DEBUG1 (("no in_cset present\n")) ;
    }
    else
    {
	for ( i = 0, n_cset = 0 ; i < n_col ; i++ )
	{
	    if (in_cset[i] < 0 || in_cset[i] >= n_col )
	    {
		stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_invalid_cset ;
		DEBUG0 (("ccolamd: malformed in_cset \n")) ;
		return (FALSE) ;
	    }
	    n_cset = MAX( n_cset, in_cset[i] ) ;
	    csize[in_cset[i]]++;
	}
	/* cset is zero based */
	n_cset++;
    }

    ASSERT ((n_cset > 0) && (n_cset<=n_col)) ;

    cset_start[0] = temp_cstart[0] = 0 ;
    for ( i = 1 ; i <= n_cset ; i++ )
    {
	cset_start[i] = cset_start[i-1] + csize[i-1] ;
	DEBUG4 ((" cset_start["ID"] = "ID" \n", i , cset_start[i])) ;
	temp_cstart[i] = cset_start[i];
    }

    /* do in reverse order to encourage natural tie-breaking */
    if (in_cset == (int *) NULL)
    {
	for (i = n_col-1 ; i >= 0 ; i-- )
	{
	    cset[ temp_cstart[0]++ ] = i ;
	}
    }
    else
    {
	for (i = n_col-1 ; i >= 0 ; i-- )
	{
	    cset[ temp_cstart[in_cset[i]]++ ] = i ;
	}
    }

    /* ] temp_cstart and csize are no longer used */

    /* ==================================================================== ] */

    /* === Construct the row and column data structures ===================== */

    if (!init_rows_cols (n_row, n_col, Row, Col, A, p, stats))
    {
	/* input matrix is invalid */
	DEBUG0 (("ccolamd: Matrix invalid\n")) ;
	return (FALSE) ;
    }

    /* [---------------- */
    /* Added for ccolamd */
    /* === Initialize front info ============================== */
    for (col = 0 ; col < n_col ; col++)
    {
    	Front_npivcol [col] = 0 ;
    	Front_nrows [col] = 0 ;
    	Front_ncols [col] = 0 ;
    	Front_parent [col] = EMPTY ;
    	Front_cols [col] = EMPTY ;
    }
    /* ------------------ ] */

    /* === Initialize scores, kill dense rows/columns ======================= */

    init_scoring (n_row, n_col, Row, Col, A, p, knobs,
	&n_row2, &n_col2, &max_deg
	/* [------------------ */
	/* Added for ccolamd */
	, in_cset, n_cset, cset_start, dead_cols
	/* Added from UMFPACK for ccolamd */
	, &ndense_row, &nempty_row, &nnewlyempty_row
	, &ndense_col, &nempty_col, &nnewlyempty_col
	/* added for symamd/colamd dense column detection */
	, max_col_degree
	) ;
	/*  ------------------]  */

    ASSERT (n_row2 == n_row - nempty_row - nnewlyempty_row - ndense_row) ;
    ASSERT (n_col2 == n_col - nempty_col - nnewlyempty_col - ndense_col) ;
    DEBUG1 (("# dense rows %d cols %d\n", ndense_row, ndense_col)) ;

    /* === Order the supercolumns =========================================== */

    ngarbage = find_ordering (n_row, n_col, Alen, Row, Col, A, p,
#ifndef NDEBUG
	n_col2,
#endif
	max_deg, 2*nnz
	/* [ ------------------ */
	/* Added for ccolamd */
	, cset, cset_start,
#ifndef NDEBUG
	n_cset,
#endif
	in_cset
	/* Added from UMFPACK for ccolamd */
	, Front_npivcol, Front_nrows, Front_ncols, Front_parent, Front_cols
	, &nfr, aggressive, InFront, fact_type
	) ;
	/*  ------------------]  */

	ASSERT( Alen >= 5*n_col ) ;
    /* === order_children deleted for ccolamd =============================== */

    /* [------------------------------ */
    /* added for ccolamd from UMFPACK: */

    /* A is no longer needed, so use A [0..5*nfr-1] as workspace [ [ */
    /* This step requires Alen >= 5*n_col */
    Front_child   = A ;
    Front_sibling = Front_child + nfr ;
    Front_stack   = Front_sibling + nfr ;
    Front_order   = Front_stack + nfr ;
    Front_size    = Front_order + nfr ;

    CCOLAMD_fsize (nfr, Front_size, Front_nrows, Front_ncols,
            Front_parent, Front_npivcol) ;

    /* --------------- */
    /* in_cset and Front_cols added for ccolamd */
    CCOLAMD_postorder (nfr, Front_parent, Front_npivcol, Front_size,
        Front_order, Front_child, Front_sibling, Front_stack, Front_cols,
	in_cset) ;
    /* --------------- */

    /* Front_size, Front_stack, Front_child, Front_sibling no longer needed ] */

    /* use A [0..nfr-1] as workspace */
    CCOLAMD_apply_order (Front_npivcol, Front_order, A, nfr, nfr) ;
    CCOLAMD_apply_order (Front_nrows,   Front_order, A, nfr, nfr) ;
    CCOLAMD_apply_order (Front_ncols,   Front_order, A, nfr, nfr) ;
    CCOLAMD_apply_order (Front_parent,  Front_order, A, nfr, nfr) ;
    CCOLAMD_apply_order (Front_cols,    Front_order, A, nfr, nfr) ;

    /* fix the parent to refer to the new numbering */
    for (i = 0 ; i < nfr ; i++)
    {
        parent = Front_parent [i] ;
        if (parent != EMPTY)
        {
            Front_parent [i] = Front_order [parent] ;
        }
    }

    /* fix InFront to refer to the new numbering */
    for (row = 0 ; row < n_row ; row++)
    {
        i = InFront [row] ;
        ASSERT (i >= EMPTY && i < nfr) ;
        if (i != EMPTY)
        {
            InFront [row] = Front_order [i] ;
        }
    }

    /* Front_order longer needed ] */

    /* === Order the columns in the fronts ================================== */

    /* use A [0..n_col-1] as inverse permutation */
    for (i = 0 ; i < n_col ; i++)
    {
        A [i] = EMPTY ;
    }

    k = 0 ;
    set1 = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
        ASSERT (Front_npivcol [i] > 0) ;

        /* ----------------- */
        /* Added for ccolamd */

	set2 = IN_CSET (Front_cols[i]) ;
        while ( set1 < set2 )
        {
            k += dead_cols[set1] ;
            DEBUG3 (("Skip null/dense columns of set "ID"\n",set1)) ;
            set1++ ;
        }
        set1 = set2 ;

        /* ----------------- */

        for (col = Front_cols [i] ; col != EMPTY ; col = Col [col].nextcol)
        {
            ASSERT (col >= 0 && col < n_col) ;
            DEBUG1 (("Colamd output ordering: k "ID" col "ID"\n", k, col)) ;
            p [k] = col ;
            ASSERT (A [col] == EMPTY) ;

            /* ----------------- */
            /* Added for ccolamd */
	    cs = IN_CSET(col) ;
            ASSERT (k >= cset_start[cs] && k < cset_start[cs+1]) ;
            /* ----------------- */

            A [col] = k ;
            k++ ;
        }
    }


    /* === Order the "dense" and null columns =============================== */

    /* ----------------- */
    /* ccolamd : no longer valid */
    /* ASSERT (k == n_col2) ; */
    /* ----------------- */

    if (n_col2 < n_col)
    {
        for (col = 0 ; col < n_col ; col++)
        {
            if (A [col] == EMPTY)
            {
                k = Col [col].shared2.order ;

                /* ----------------- */
                /* ccolamd : no longer valid */
                /* ASSERT (k >= n_col2 && k < n_col) ; */
                /* ----------------- */

                /* ----------------- */
                /* Added for ccolamd */
		cs = IN_CSET(col) ;
#ifndef NDEBUG
                dead_cols[cs]-- ;
#endif
                ASSERT (k >= cset_start[cs] && k < cset_start[cs+1]) ;
                /* ----------------- */

                DEBUG1 (("Colamd output ordering: k "ID" col "ID
                    " (dense or null col)\n", k, col)) ;
                p [k] = col ;
                A [col] = k ;
            }
        }
    }

    /* ----------------- */
    /* Added for ccolamd */
#ifndef NDEBUG
    for (i = 0 ; i < n_cset ; i++)
    {
    	ASSERT (dead_cols[i] == 0) ;
    }
#endif
    /* ----------------- */

    /* -----------------] */
    /* === Return statistics in stats ======================================= */

    /* ------------------ */
    /* modified for UMFPACK */
    stats [CCOLAMD_DENSE_ROW] = ndense_row ;
    stats [CCOLAMD_DENSE_COL] = nempty_row ;
    stats [CCOLAMD_NEWLY_EMPTY_ROW] = nnewlyempty_row ;
    stats [CCOLAMD_DENSE_COL] = ndense_col ;
    stats [CCOLAMD_EMPTY_COL] = nempty_col ;
    stats [CCOLAMD_NEWLY_EMPTY_COL] = nnewlyempty_col ;
    ASSERT (ndense_col + nempty_col + nnewlyempty_col == n_col - n_col2) ;
    if ( p_nfr )
    {
    	*p_nfr = nfr ;
    }
    /* ------------------ */
    stats [CCOLAMD_DEFRAG_COUNT] = ngarbage ;
    DEBUG0 (("ccolamd: done.\n")) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === colamd_report ======================================================== */
/* ========================================================================== */

PUBLIC void CCOLAMD_report
(
    Int stats [CCOLAMD_STATS]
)
{
    print_report ("ccolamd", stats) ;
}


/* ========================================================================== */
/* === symamd_report ======================================================== */
/* ========================================================================== */

PUBLIC void CSYMAMD_report
(
    Int stats [CCOLAMD_STATS]
)
{
    print_report ("csymamd", stats) ;
}



/* ========================================================================== */
/* === NON-USER-CALLABLE ROUTINES: ========================================== */
/* ========================================================================== */

/* There are no user-callable routines beyond this point in the file */


/* ========================================================================== */
/* === init_rows_cols ======================================================= */
/* ========================================================================== */

/*
    Takes the column form of the matrix in A and creates the row form of the
    matrix.  Also, row and column attributes are stored in the Col and Row
    structs.  If the columns are un-sorted or contain duplicate row indices,
    this routine will also sort and remove duplicate row indices from the
    column form of the matrix.  Returns FALSE if the matrix is invalid,
    TRUE otherwise.  Not user-callable.
*/

PRIVATE Int init_rows_cols	/* returns TRUE if OK, or FALSE otherwise */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Colamd_Row Row [],		/* of size n_row+1 */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* row indices of A, of size Alen */
    Int p [],			/* pointers to columns in A, of size n_col+1 */
    Int stats [CCOLAMD_STATS]	/* colamd statistics */
)
{
    /* === Local variables ================================================== */

    Int col ;			/* a column index */
    Int row ;			/* a row index */
    Int *cp ;			/* a column pointer */
    Int *cp_end ;		/* a pointer to the end of a column */
    Int *rp ;			/* a row pointer */
    Int *rp_end ;		/* a pointer to the end of a row */
    Int last_row ;		/* previous row */

    /* === Initialize columns, and check column pointers ==================== */

    for (col = 0 ; col < n_col ; col++)
    {
	Col [col].start = p [col] ;
	Col [col].length = p [col+1] - p [col] ;

	if (Col [col].length < 0)
	{
	    /* column pointers must be non-decreasing */
	    stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_col_length_negative ;
	    stats [CCOLAMD_INFO1] = col ;
	    stats [CCOLAMD_INFO2] = Col [col].length ;
	    DEBUG0 (("ccolamd: col "ID" length "ID" < 0\n", col, Col [col].length)) ;
	    return (FALSE) ;
	}

	Col [col].shared1.thickness = 1 ;
	Col [col].shared2.score = 0 ;
	Col [col].shared3.prev = EMPTY ;
	Col [col].shared4.degree_next = EMPTY ;

        /* ------------------ */
        /* added for UMFPACK: */
        Col [col].nextcol = EMPTY ;
        Col [col].lastcol = col ;
        /* ------------------ */

    }

    /* p [0..n_col] no longer needed, used as "head" in subsequent routines */

    /* === Scan columns, compute row degrees, and check row indices ========= */

    stats [CCOLAMD_INFO3] = 0 ;	/* number of duplicate or unsorted row indices*/

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].length = 0 ;
	Row [row].shared2.mark = -1 ;
        /* ------------------ */
        /* added for ccolamd & UMFPACK: */
        Row [row].thickness = 1 ;
        Row [row].front = EMPTY ;
        /* ------------------ */
    }

    for (col = 0 ; col < n_col ; col++)
    {
	DEBUG1 (("\nCcolamd input column "ID":\n", col)) ;
	last_row = -1 ;

	cp = &A [p [col]] ;
	cp_end = &A [p [col+1]] ;

	while (cp < cp_end)
	{
	    row = *cp++ ;
	    DEBUG1 (("row: "ID"\n", row)) ;

	    /* make sure row indices within range */
	    if (row < 0 || row >= n_row)
	    {
		stats [CCOLAMD_STATUS] = CCOLAMD_ERROR_row_index_out_of_bounds ;
		stats [CCOLAMD_INFO1] = col ;
		stats [CCOLAMD_INFO2] = row ;
		stats [CCOLAMD_INFO3] = n_row ;
		DEBUG0 (("ccolamd: row "ID" col "ID" out of bounds\n", row, col)) ;
		return (FALSE) ;
	    }

	    if (row <= last_row || Row [row].shared2.mark == col)
	    {
		/* row index are unsorted or repeated (or both), thus col */
		/* is jumbled.  This is a notice, not an error condition. */
		stats [CCOLAMD_STATUS] = CCOLAMD_OK_BUT_JUMBLED ;
		stats [CCOLAMD_INFO1] = col ;
		stats [CCOLAMD_INFO2] = row ;
		(stats [CCOLAMD_INFO3]) ++ ;
		DEBUG1 (("ccolamd: row "ID" col "ID" unsorted/duplicate\n",row,col));
	    }

	    if (Row [row].shared2.mark != col)
	    {
		Row [row].length++ ;
	    }
	    else
	    {
		/* this is a repeated entry in the column, */
		/* it will be removed */
		Col [col].length-- ;
	    }

	    /* mark the row as having been seen in this column */
	    Row [row].shared2.mark = col ;

	    last_row = row ;
	}
    }

    /* === Compute row pointers ============================================= */

    /* row form of the matrix starts directly after the column */
    /* form of matrix in A */
    Row [0].start = p [n_col] ;
    Row [0].shared1.p = Row [0].start ;
    Row [0].shared2.mark = -1 ;
    for (row = 1 ; row < n_row ; row++)
    {
	Row [row].start = Row [row-1].start + Row [row-1].length ;
	Row [row].shared1.p = Row [row].start ;
	Row [row].shared2.mark = -1 ;
    }

    /* === Create row form ================================================== */

    if (stats [CCOLAMD_STATUS] == CCOLAMD_OK_BUT_JUMBLED)
    {
	/* if cols jumbled, watch for repeated row indices */
 	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		row = *cp++ ;
		if (Row [row].shared2.mark != col)
		{
		    A [(Row [row].shared1.p)++] = col ;
		    Row [row].shared2.mark = col ;
		}
	    }
	}
    }
    else
    {
	/* if cols not jumbled, we don't need the mark (this is faster) */
	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		A [(Row [*cp++].shared1.p)++] = col ;
	    }
	}
    }

    /* === Clear the row marks and set row degrees ========================== */

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].shared2.mark = 0 ;
	Row [row].shared1.degree = Row [row].length ;
    }

    /* === See if we need to re-create columns ============================== */

    if (stats [CCOLAMD_STATUS] == CCOLAMD_OK_BUT_JUMBLED)
    {
    	DEBUG0 (("ccolamd: reconstructing column form, matrix jumbled\n")) ;

#ifndef NDEBUG
	/* make sure column lengths are correct */
 	for (col = 0 ; col < n_col ; col++)
	{
	    p [col] = Col [col].length ;
	}
	for (row = 0 ; row < n_row ; row++)
	{
	    rp = &A [Row [row].start] ;
	    rp_end = rp + Row [row].length ;
	    while (rp < rp_end)
	    {
		p [*rp++]-- ;
	    }
	}
	for (col = 0 ; col < n_col ; col++)
	{
	    ASSERT (p [col] == 0) ;
	}
	/* now p is all zero (different than when debugging is turned off) */
#endif /*  NDEBUG */

	/* === Compute col pointers ========================================= */

	/* col form of the matrix starts at A [0]. */
	/* Note, we may have a gap between the col form and the row */
	/* form if there were duplicate entries, if so, it will be */
	/* removed upon the first garbage collection */
	Col [0].start = 0 ;
	p [0] = Col [0].start ;
	for (col = 1 ; col < n_col ; col++)
	{
	    /* note that the lengths here are for pruned columns, i.e. */
	    /* no duplicate row indices will exist for these columns */
	    Col [col].start = Col [col-1].start + Col [col-1].length ;
	    p [col] = Col [col].start ;
	}

	/* === Re-create col form =========================================== */

	for (row = 0 ; row < n_row ; row++)
	{
	    rp = &A [Row [row].start] ;
	    rp_end = rp + Row [row].length ;
	    while (rp < rp_end)
	    {
		A [(p [*rp++])++] = row ;
	    }
	}
    }

    /* === Done.  Matrix is not (or no longer) jumbled ====================== */


    return (TRUE) ;
}


/* ========================================================================== */
/* === init_scoring ========================================================= */
/* ========================================================================== */

/*
    Kills dense or empty columns and rows, calculates an initial score for
    each column, and places all columns in the degree lists.  Not user-callable.
*/

PRIVATE void init_scoring
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Colamd_Row Row [],		/* of size n_row+1 */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* column form and row form of A */
    Int head [],		/* of size n_col+1 */
    double knobs [CCOLAMD_KNOBS],/* parameters */
    Int *p_n_row2,		/* number of non-dense, non-empty rows */
    Int *p_n_col2,		/* number of non-dense, non-empty columns */
    Int *p_max_deg		/* maximum row degree */
    /* ------------------ */
    /* added for ccolamd */
    , Int in_cset []
    , Int n_cset
    , Int cset_start []
    , Int dead_cols []
    /* added from UMFPACK for ccolamd */
    , Int *p_ndense_row         /* number of dense rows */
    , Int *p_nempty_row         /* number of original empty rows */
    , Int *p_nnewlyempty_row    /* number of newly empty rows */
    , Int *p_ndense_col         /* number of dense cols (excl "empty" cols) */
    , Int *p_nempty_col         /* number of original empty cols */
    , Int *p_nnewlyempty_col    /* number of newly empty cols */
    /* ------------------ */
    , Int max_col_degree
)
{
/* === Local variables ================================================== */

    Int c ;			/* a column index */
    Int r, row ;		/* a row index */
    Int *cp ;			/* a column pointer */
    Int deg ;			/* degree of a row or column */
    Int *cp_end ;		/* a pointer to the end of a column */
    Int *new_cp ;		/* new column pointer */
    Int col_length ;		/* length of pruned column */
    Int score ;			/* current column score */
    Int n_col2 ;		/* number of non-dense, non-empty columns */
    Int n_row2 ;		/* number of non-dense, non-empty rows */
    Int dense_row_count ;	/* remove rows with more entries than this */
    Int dense_col_count ;	/* remove cols with more entries than this */
    Int max_deg ;		/* maximum row degree */
    /* ------------------ */
    /* added for ccolamd*/
    Int s ;			/* a cset index */
    /* added from UMFPACK for ccolamd */
    Int ndense_row ;            /* number of dense rows */
    Int nempty_row ;            /* number of empty rows */
    Int nnewlyempty_row ;       /* number of newly empty rows */
    Int ndense_col ;            /* number of dense cols (excl "empty" cols) */
    Int nempty_col ;            /* number of original empty cols */
    Int nnewlyempty_col ;       /* number of newly empty cols */
    Int ne ;
    /* ------------------ */


#ifndef NDEBUG
    Int debug_count ;		/* debug only. */
#endif /* NDEBUG */

    /* === Extract knobs ==================================================== */

    /* --------------------- */
    /* modified for ccolamd*/
    /* old dense row/column knobs:
     */

    /* TODO determine knob function for ccolamd */

    /* COLAMD: */
#ifdef COLAMD_DENSE_KNOBS
    dense_row_count = MAX (0, MIN (knobs [CCOLAMD_DENSE_ROW] * n_col, n_col)) ;
    dense_col_count = MAX (0, MIN (knobs [CCOLAMD_DENSE_COL] * n_row, n_row)) ;
#else
    /* new, for ccolamd from UMFPACK: */
    /* Note: if knobs contains a NaN, this is undefined: */
    dense_row_count =
        UMFPACK_DENSE_DEGREE_THRESHOLD (knobs [CCOLAMD_DENSE_ROW], n_col) ;
    dense_col_count =
        UMFPACK_DENSE_DEGREE_THRESHOLD (knobs [CCOLAMD_DENSE_COL],
	    max_col_degree) ;
    /* Make sure dense_*_count is between 0 and n: */
    dense_row_count = MAX (0, MIN (dense_row_count, n_col)) ;
    dense_col_count = MAX (0, MIN (dense_col_count, max_col_degree)) ;
    /* --------------------- */
#endif

    DEBUG1 (("ccolamd: densecount: "ID" "ID"\n", dense_row_count, dense_col_count)) ;
    max_deg = 0 ;

    n_col2 = n_col ;
    n_row2 = n_row ;

    /* [ --------------- */
    /* added for ccolamd */
    /* Set the head array for bookkeeping of dense and empty columns. */
    /* This will be used as hash buckets later. */
    for (s = 0 ; s < n_cset ; s++)
    {
   	 head [s] = cset_start [s+1] ;
    }

    /* added for ccolamd from UMFPACK */
    ndense_col = 0 ;
    nempty_col = 0 ;
    nnewlyempty_col = 0 ;
    ndense_row = 0 ;
    nempty_row = 0 ;
    nnewlyempty_row = 0 ;
    /* -------------- ] */



    /* === Kill empty columns =============================================== */

    /* Put the empty columns at the end in their natural order, so that LU */
    /* factorization can proceed as far as possible. */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	deg = Col [c].length ;
	if (deg == 0)
	{
	    /* this is a empty column, kill and order it last of its cset */
	    /* --------------------- */
	    /* modified for ccolamd*/
	    Col [c].shared2.order = --head [IN_CSET(c)] ;
	    --n_col2 ;
	    dead_cols[IN_CSET(c)] ++ ;
	    /* added for ccolamd from UMFPACK */
	    nempty_col++ ;
    	    /* --------------------- */
	    KILL_PRINCIPAL_COL (c) ;
	}
    }
    DEBUG1 (("ccolamd: null columns killed: "ID"\n", n_col - n_col2)) ;

    /* === Kill dense columns =============================================== */

    /* Put the dense columns at the end, in their natural order */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* skip any dead columns */
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	deg = Col [c].length ;
	if (deg > dense_col_count)
	{
	    /* this is a dense column, kill and order it last of its cset */
   	    /* --------------------- */
   	    /* modified for ccolamd*/
	    Col [c].shared2.order = --head [IN_CSET(c)] ;
	    --n_col2 ;
	    dead_cols[IN_CSET(c)] ++ ;
            /* added for ccolamd from UMFPACK */
            ndense_col++ ;
	    /* --------------------- */
	    /* decrement the row degrees */
	    cp = &A [Col [c].start] ;
	    cp_end = cp + Col [c].length ;
	    while (cp < cp_end)
	    {
		Row [*cp++].shared1.degree-- ;
	    }
	    KILL_PRINCIPAL_COL (c) ;
	}
    }
    DEBUG1 (("ccolamd: Dense and null columns killed: "ID"\n", n_col - n_col2)) ;

    /* === Kill dense and empty rows ======================================== */

    /* Note that there can now be empty rows, since dense columns have
     * been deleted.  These are "newly" empty rows. */

    ne = 0 ;
    for (r = 0 ; r < n_row ; r++)
    {
	deg = Row [r].shared1.degree ;
	ASSERT (deg >= 0 && deg <= n_col) ;
        /* --------------------- */
        /* added for ccolamd UMFPACK */
        if (deg > dense_row_count)
        {
            /* There is at least one dense row.  Continue ordering, but */
            /* symbolic factorization will be redone after ccolamd is done.*/
            ndense_row++ ;
        }
        if (deg == 0)
        {
            /* this is a newly empty row, or original empty row */
            ne++ ;
        }
        /* --------------------- */

	if (deg > dense_row_count || deg == 0)
	{
	    /* kill a dense or empty row */
	    KILL_ROW (r) ;
   	    /* --------------------- */
   	    /* added for ccolamd from UMFPACK */
	    Row [r].thickness = 0 ;
            /* --------------------- */
	    --n_row2 ;
	}
	else
	{
	    /* keep track of max degree of remaining rows */
	    max_deg = MAX (max_deg, deg) ;
	}
    }
    nnewlyempty_row = ne - nempty_row ;
    DEBUG1 (("ccolamd: Dense and null rows killed: "ID"\n", n_row - n_row2)) ;

    /* === Compute initial column scores ==================================== */

    /* At this point the row degrees are accurate.  They reflect the number */
    /* of "live" (non-dense) columns in each row.  No empty rows exist. */
    /* Some "live" columns may contain only dead rows, however.  These are */
    /* pruned in the code below. */

    /* now find the initial COLMMD score for each column */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* skip dead column */
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	score = 0 ;
	cp = &A [Col [c].start] ;
	new_cp = cp ;
	cp_end = cp + Col [c].length ;
	while (cp < cp_end)
	{
	    /* get a row */
	    row = *cp++ ;
	    /* skip if dead */
	    if (ROW_IS_DEAD (row))
	    {
		continue ;
	    }
	    /* compact the column */
	    *new_cp++ = row ;
	    /* add row's external degree */
	    score += Row [row].shared1.degree - 1 ;
	    /* guard against integer overflow */
	    score = MIN (score, n_col) ;
	}
	/* determine pruned column length */
	col_length = (Int) (new_cp - &A [Col [c].start]) ;
	if (col_length == 0)
	{
	    /* a newly-made null column (all rows in this col are "dense" */
	    /* and have already been killed) */
	    DEBUG1 (("Newly null killed: "ID"\n", c)) ;
	    /* --------------------- */
	    /* modified for ccolamd*/
	    Col [c].shared2.order = -- head [IN_CSET(c)] ;
	    --n_col2 ;
	    dead_cols[IN_CSET(c)] ++ ;
            /* added for ccolamd from UMFPACK */
            nnewlyempty_col++ ;
    	    /* --------------------- */
	    KILL_PRINCIPAL_COL (c) ;
	}
	else
	{
	    /* set column length and set score */
	    ASSERT (score >= 0) ;
	    ASSERT (score <= n_col) ;
	    Col [c].length = col_length ;
	    Col [c].shared2.score = score ;

	}
    }
    DEBUG1 (("ccolamd: Dense, null, and newly-null columns killed: "ID"\n",
    	n_col-n_col2)) ;

    /* At this point, all empty rows and columns are dead.  All live columns */
    /* are "clean" (containing no dead rows) and simplicial (no supercolumns */
    /* yet).  Rows may contain dead columns, but all live rows contain at */
    /* least one live column. */

#ifndef NDEBUG
    debug_count = 0 ;
#endif /* NDEBUG */

    /* clear the hash buckets */
    for (c = 0 ; c <= n_col ; c++)
    {
	head [c] = EMPTY ;
    }

/* === Initialize degree lists removed for ccolamd ========================== */
#ifndef NDEBUG

    /* --------------------- */
    /* modified for ccolamd*/
    debug_structures (n_row, n_col, Row, Col, A, in_cset, cset_start) ;
    /* ASSERT (debug_count == n_col2) ; TBD : Add similar in find_ordering*/
    /* --------------------- */

#endif /* NDEBUG */

    /* === Return number of remaining columns, and max row degree =========== */

    *p_n_col2 = n_col2 ;
    *p_n_row2 = n_row2 ;
    *p_max_deg = max_deg ;

    /* --------------------- */
    /* added for ccolamd from UMFPACK */
    *p_ndense_row = ndense_row ;
    *p_nempty_row = nempty_row ;        /* original empty rows */
    *p_nnewlyempty_row = nnewlyempty_row ;
    *p_ndense_col = ndense_col ;
    *p_nempty_col = nempty_col ;        /* original empty cols */
    *p_nnewlyempty_col = nnewlyempty_col ;
    /* --------------------- */

}


/* ========================================================================== */
/* === find_ordering ======================================================== */
/* ========================================================================== */

/*
    Order the principal columns of the supercolumn form of the matrix
    (no supercolumns on input).  Uses a minimum approximate column minimum
    degree ordering method.  Not user-callable.
*/

PRIVATE Int find_ordering	/* return the number of garbage collections */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Int Alen,			/* size of A, 2*nnz + n_col or larger */
    Colamd_Row Row [],		/* of size n_row+1 */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* column form and row form of A */
    Int head [],		/* of size n_col+1 */
#ifndef NDEBUG
    Int n_col2,			/* Remaining columns to order */
#endif
    Int max_deg,		/* Maximum row degree */
    Int pfree			/* index of first free slot (2*nnz on entry) */
    /* --------------------- */
    /* added for ccolamd */
    , Int cset []		/* constraint set of A */
    , Int cset_start []		/* pointer to the start of every cset */
#ifndef NDEBUG
    , Int n_cset		/* number of csets */
#endif
    , Int in_cset []		/* col -> cset mapping */
    /* added for ccolamd from UMFPACK */
    , Int Front_npivcol [ ]
    , Int Front_nrows [ ]
    , Int Front_ncols [ ]
    , Int Front_parent [ ]
    , Int Front_cols [ ]
    , Int *p_nfr                /* number of fronts */
    , Int aggressive
    , Int InFront [ ]
    , Int fact_type
    /* ------------------ */
)
{
    /* === Local variables ================================================== */

    Int k ;			/* current pivot ordering step */
    Int pivot_col ;		/* current pivot column */
    Int *cp ;			/* a column pointer */
    Int *rp ;			/* a row pointer */
    Int pivot_row ;		/* current pivot row */
    Int *new_cp ;		/* modified column pointer */
    Int *new_rp ;		/* modified row pointer */
    Int pivot_row_start ;	/* pointer to start of pivot row */
    Int pivot_row_degree ;	/* number of columns in pivot row */
    Int pivot_row_length ;	/* number of supercolumns in pivot row */
    Int pivot_col_score ;	/* score of pivot column */
    Int needed_memory ;		/* free space needed for pivot row */
    Int *cp_end ;		/* pointer to the end of a column */
    Int *rp_end ;		/* pointer to the end of a row */
    Int row ;			/* a row index */
    Int col ;			/* a column index */
    Int max_score ;		/* maximum possible score */
    Int cur_score ;		/* score of current column */
    unsigned Int hash ;		/* hash value for supernode detection */
    Int head_column ;		/* head of hash bucket */
    Int first_col ;		/* first column in hash bucket */
    Int tag_mark ;		/* marker value for mark array */
    Int row_mark ;		/* Row [row].shared2.mark */
    Int set_difference ;	/* set difference size of row with pivot row */
    Int min_score ;		/* smallest column score */
    Int col_thickness ;		/* "thickness" (no. of columns in a supercol) */
    Int max_mark ;		/* maximum value of tag_mark */
    Int pivot_col_thickness ;	/* number of columns represented by pivot col */
    Int prev_col ;		/* Used by Dlist operations. */
    Int next_col ;		/* Used by Dlist operations. */
    Int ngarbage ;		/* number of garbage collections performed */
    /* --------------------- */
    /* added for ccolamd */
    Int current_set ;		/* consraint set that is being ordered */
    Int score ;			/* score of a column */
    Int colstart ;		/* pointer to first column in current cset */
    Int colend ;		/* pointer to last column in current cset */
    Int deadcol ;		/* number of dense & null columns in a cset */
    /* --------------------- */

#ifndef NDEBUG
    Int debug_d ;		/* debug loop counter */
    Int debug_step = 0 ;	/* debug loop counter */
    Int cols_thickness = 0 ;	/* the thickness of the columns in current */
    				/* cset degreelist and in pivot row pattern. */
#endif /* NDEBUG */

    /* ------------------ */
    /* added for ccolamd from UMFPACK : */
    Int pivot_row_thickness ;   /* number of rows represented by pivot row */
    Int nfr = 0 ;               /* number of fronts */
    Int child ;
    /* ------------------ */


    /* === Initialization and clear mark ==================================== */

    max_mark = Int_MAX - n_col ;	/* Int_MAX defined in <limits.h> */
    tag_mark = clear_mark (0, max_mark, n_row, Row) ;
    min_score = 0 ;
    ngarbage = 0 ;
    current_set = -1 ;
    deadcol = 0 ;
    DEBUG1 (("ccolamd: Ordering, n_col2="ID"\n", n_col2)) ;

    /* ------------------ */
    /* added for ccolamd from UMFPACK : */
    for (row = 0 ; row < n_row ; row++)
    {
        InFront [row] = EMPTY ;
    }
    /* ------------------ */

    /* === Order the columns ================================================ */

    /* loop terminating condition modified for ccolamd */
    for (k = 0 ; k < n_col ; /* 'k' is incremented below */)
    {

	/* make sure degree list isn't empty */
	ASSERT (min_score >= 0) ;
	ASSERT (min_score <= n_col) ;
	ASSERT (head [min_score] >= EMPTY) ;

#ifndef NDEBUG
	for (debug_d = 0 ; debug_d < min_score ; debug_d++)
	{
	    ASSERT (head [debug_d] == EMPTY) ;
	}
#endif /* NDEBUG */


	/* ------------------- */
	/* [ added for ccolamd */
	/* Initialize the degree list with columns from next non-empty cset */

	while ( (k+deadcol) == cset_start[current_set+1])
	{
	    current_set++ ;
	    DEBUG1 (("\n\n\n============ CSET: "ID"\n", current_set)) ;
	    k += deadcol ; /* jump to start of next cset */
  	    deadcol = 0 ; /* reset dead column count */

	    ASSERT( (current_set == n_cset) == (k == n_col) ) ;
	    /* return if all columns are ordered. */
	    if ( k == n_col )
	    {
		*p_nfr = nfr ;
	    	return (ngarbage) ;
	    }

#ifndef NDEBUG
	    for (col = 0 ; col <= n_col ; col++)
	    {
	        ASSERT (head [col] == EMPTY) ;
	    }
#endif /* NDEBUG */

	    min_score = n_col ;
	    colstart =  cset_start[current_set] ;
	    colend =  cset_start[current_set+1] ;

	    while ( colstart < colend )
	    {
		col = cset[colstart++] ;

		if (COL_IS_DEAD(col))
		{
		    DEBUG1 (("Column "ID" is dead\n", col)) ;
		    /* count dense and null columns */
		    if ( Col [col].shared2.order != EMPTY )
		    {
			    deadcol++ ;
		    }
		    continue ;
		}

		/* only add principal columns in current set to degree lists */
		ASSERT (IN_CSET(col) == current_set) ;

		score = Col [col].shared2.score ;
		DEBUG1 (("Column "ID" is alive, score %d\n", col, score)) ;

		ASSERT (min_score >= 0) ;
		ASSERT (min_score <= n_col) ;
		ASSERT (score >= 0) ;
		ASSERT (score <= n_col) ;
		ASSERT (head [score] >= EMPTY) ;

		/* now add this column to dList at proper score location */
		next_col = head [score] ;
		Col [col].shared3.prev = EMPTY ;
		Col [col].shared4.degree_next = next_col ;

		/* if there already was a column with the same score, set its */
		/* previous pointer to this new column */
		if (next_col != EMPTY)
		{
			Col [next_col].shared3.prev = col ;
		}
		head [score] = col ;

		/* see if this score is less than current min */
		min_score = MIN (min_score, score) ;

	    }

#ifndef NDEBUG
	    DEBUG1 (("degree lists initialized \n")) ;
	    debug_deg_lists (n_row, n_col, Row, Col, head, min_score,
	    ((cset_start[current_set+1]-cset_start[current_set])-deadcol),
	    max_deg) ;
#endif /* NDEBUG */
	}

#ifndef NDEBUG
	if (debug_step % 100 == 0)
	{
	    DEBUG2 (("\n...       Step k: "ID" out of n_col2: "ID"\n", k, n_col2)) ;
	}
	else
	{
	    DEBUG3 (("\n----------Step k: "ID" out of n_col2: "ID"\n", k, n_col2)) ;
	}
	debug_step++ ;
	DEBUG1 (("start of step k="ID": ", k)) ;
	/* fix should here */
	debug_deg_lists (n_row, n_col, Row, Col, head,
	     min_score, cset_start[current_set+1]-(k+deadcol), max_deg) ;
	debug_matrix (n_row, n_col, Row, Col, A) ;
#endif /* NDEBUG */
	/* ------------------- ] */

	/* === Select pivot column, and order it ============================ */

	while (head [min_score] == EMPTY && min_score < n_col)
	{
	    min_score++ ;
	}

	pivot_col = head [min_score] ;

	ASSERT (pivot_col >= 0 && pivot_col <= n_col) ;
	next_col = Col [pivot_col].shared4.degree_next ;
	head [min_score] = next_col ;
	if (next_col != EMPTY)
	{
	    Col [next_col].shared3.prev = EMPTY ;
	}

	ASSERT (COL_IS_ALIVE (pivot_col)) ;
	DEBUG3 (("Pivot col: "ID"\n", pivot_col)) ;

	/* remember score for defrag check */
	pivot_col_score = Col [pivot_col].shared2.score ;

	/* the pivot column is the kth column in the pivot order */
	Col [pivot_col].shared2.order = k ;

	/* increment order count by column thickness */
	pivot_col_thickness = Col [pivot_col].shared1.thickness ;
	k += pivot_col_thickness ;
	ASSERT (pivot_col_thickness > 0) ;

	/* === Garbage_collection, if necessary ============================= */

	needed_memory = MIN (pivot_col_score, n_col - k) ;
	if (pfree + needed_memory >= Alen)
	{
	    pfree = garbage_collection (n_row, n_col, Row, Col, A, &A [pfree]) ;
	    ngarbage++ ;
	    /* after garbage collection we will have enough */
	    ASSERT (pfree + needed_memory < Alen) ;
	    /* garbage collection has wiped out the Row[].shared2.mark array */
	    tag_mark = clear_mark (0, max_mark, n_row, Row) ;

#ifndef NDEBUG
	    debug_matrix (n_row, n_col, Row, Col, A) ;
#endif /* NDEBUG */
	}

	/* === Compute pivot row pattern ==================================== */

	/* get starting location for this new merged row */
	pivot_row_start = pfree ;

	/* initialize new row counts to zero */
	pivot_row_degree = 0 ;

        /* ------------------ */
        /* added for ccolamd from UMFPACK: */
        pivot_row_thickness = 0 ;
        /* ------------------ */

	/* tag pivot column as having been visited so it isn't included */
	/* in merged pivot row */
	Col [pivot_col].shared1.thickness = -pivot_col_thickness ;

	/* pivot row is the union of all rows in the pivot column pattern */
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    /* get a row */
	    row = *cp++ ;
	    ASSERT (row >= 0 && row < n_row) ;
	    DEBUG4 (("Pivcol pattern "ID" "ID"\n", ROW_IS_ALIVE (row), row)) ;
	    /* skip if row is dead */
	    if (ROW_IS_ALIVE (row))
	    {

		/* ------------------ */
		/* added for ccolamd from UMFPACK: */
		/* sum the thicknesses of all the rows */
		/* ASSERT (Row [row].thickness > 0) ; */
		pivot_row_thickness += Row [row].thickness ;
		/* ------------------ */

		rp = &A [Row [row].start] ;
		rp_end = rp + Row [row].length ;
		while (rp < rp_end)
		{
		    /* get a column */
		    col = *rp++ ;
		    /* add the column, if alive and untagged */
		    col_thickness = Col [col].shared1.thickness ;
		    if (col_thickness > 0 && COL_IS_ALIVE (col))
		    {
			/* tag column in pivot row */
			Col [col].shared1.thickness = -col_thickness ;
			ASSERT (pfree < Alen) ;
			/* place column in pivot row */
			A [pfree++] = col ;
			pivot_row_degree += col_thickness ;
			/* ------------------------------ */
			/* added for ccolamd from UMFPACK: */
			DEBUG4 (("\t\t\tNew live col in pivrow: "ID"\n",col)) ;
			/* ------------------------------ */
		    }
		    /* ------------------------------ */
		    /* added for ccolamd from UMFPACK */
#ifndef NDEBUG
		    if (col_thickness < 0 && COL_IS_ALIVE (col))
		    {
			DEBUG4 (("\t\t\tOld live col in pivrow: "ID"\n",col)) ;
		    }
#endif
		    /* ------------------------------ */
		}
	    }
	}

        /* ------------------------------ */
        /* added for ccolamd from UMFPACK: */
        /* pivot_row_thickness is the number of rows in frontal matrix */
        /* both pivotal rows and nonpivotal rows */
        /* ------------------------------ */


	/* clear tag on pivot column */
	Col [pivot_col].shared1.thickness = pivot_col_thickness ;
	max_deg = MAX (max_deg, pivot_row_degree) ;

#ifndef NDEBUG
	DEBUG3 (("check2\n")) ;
	debug_mark (n_row, Row, tag_mark, max_mark) ;
#endif /* NDEBUG */

	/* === Kill all rows used to construct pivot row ==================== */

	/* also kill pivot row, temporarily */
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    /* may be killing an already dead row */
	    row = *cp++ ;
	    DEBUG3 (("Kill row in pivot col: "ID"\n", row)) ;
	    ASSERT (row >= 0 && row < n_row) ;
	    /* ------------------------------ */
	    /* added for ccolamd from UMFPACK: */
            if (ROW_IS_ALIVE (row))
            {
                if (Row [row].front != EMPTY)
                {
                    /* This row represents a frontal matrix. */
                    /* Row [row].front is a child of current front */
                    child = Row [row].front ;
                    Front_parent [child] = nfr ;
                    DEBUG1 (("Front "ID" => front "ID", normal\n", child, nfr));
                }
                else
                {
                    /* This is an original row.  Keep track of which front
                     * is its parent in the row-merge tree. */
                    InFront [row] = nfr ;
                    DEBUG1 (("Row "ID" => front "ID", normal\n", row, nfr)) ;
                }
            }
	    /* ------------------------------ */

            KILL_ROW (row) ;

           /* ------------------------------ */
           /* added for ccolamd from UMFPACK: */
            Row [row].thickness = 0 ;
	}

	/* === Select a row index to use as the new pivot row =============== */

	pivot_row_length = pfree - pivot_row_start ;
	if (pivot_row_length > 0)
	{
	    /* pick the "pivot" row arbitrarily (first row in col) */
	    pivot_row = A [Col [pivot_col].start] ;
	    DEBUG3 (("Pivotal row is "ID"\n", pivot_row)) ;
	}
	else
	{
	    /* there is no pivot row, since it is of zero length */
	    pivot_row = EMPTY ;
	    ASSERT (pivot_row_length == 0) ;
	}
	ASSERT (Col [pivot_col].length > 0 || pivot_row_length == 0) ;

	/* === Approximate degree computation =============================== */

	/* Here begins the computation of the approximate degree.  The column */
	/* score is the sum of the pivot row "length", plus the size of the */
	/* set differences of each row in the column minus the pattern of the */
	/* pivot row itself.  The column ("thickness") itself is also */
	/* excluded from the column score (we thus use an approximate */
	/* external degree). */

	/* The time taken by the following code (compute set differences, and */
	/* add them up) is proportional to the size of the data structure */
	/* being scanned - that is, the sum of the sizes of each column in */
	/* the pivot row.  Thus, the amortized time to compute a column score */
	/* is proportional to the size of that column (where size, in this */
	/* context, is the column "length", or the number of row indices */
	/* in that column).  The number of row indices in a column is */
	/* monotonically non-decreasing, from the length of the original */
	/* column on input to colamd. */

	/* === Compute set differences ====================================== */

	DEBUG3 (("** Computing set differences phase. **\n")) ;

	/* pivot row is currently dead - it will be revived later. */

	DEBUG3 (("Pivot row: ")) ;
	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    DEBUG3 (("Col: "ID"\n", col)) ;

	    /* clear tags used to construct pivot row pattern */
	    col_thickness = -Col [col].shared1.thickness ;
	    ASSERT (col_thickness > 0) ;
	    Col [col].shared1.thickness = col_thickness ;

	    /* === Remove column from degree list =========================== */

           /* ------------------------ */
           /* added check for ccolamd  */
	   /* only columns in current_set will be in degree list */
	    if ( IN_CSET(col) == current_set )
	    {
#ifndef NDEBUG
		cols_thickness += col_thickness ;
#endif /* NDEBUG */
		cur_score = Col [col].shared2.score ;
		prev_col = Col [col].shared3.prev ;
		next_col = Col [col].shared4.degree_next ;
		DEBUG3 (("            cur_score "ID" prev_col "ID" next_col "ID"\n",
			cur_score, prev_col, next_col)) ;
		ASSERT (cur_score >= 0) ;
		ASSERT (cur_score <= n_col) ;
		ASSERT (cur_score >= EMPTY) ;
		if (prev_col == EMPTY)
		{
		    head [cur_score] = next_col ;
		}
		else
		{
		    Col [prev_col].shared4.degree_next = next_col ;
		}
		if (next_col != EMPTY)
		{
		    Col [next_col].shared3.prev = prev_col ;
		}

	    }
           /* ------------------------- */

	    /* === Scan the column ========================================== */

	    cp = &A [Col [col].start] ;
	    cp_end = cp + Col [col].length ;
	    while (cp < cp_end)
	    {
		/* get a row */
		row = *cp++ ;
		row_mark = Row [row].shared2.mark ;
		/* skip if dead */
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    continue ;
		}
		ASSERT (row != pivot_row) ;
		set_difference = row_mark - tag_mark ;
		/* check if the row has been seen yet */
		if (set_difference < 0)
		{
		    ASSERT (Row [row].shared1.degree <= max_deg) ;
		    set_difference = Row [row].shared1.degree ;
		}
		/* subtract column thickness from this row's set difference */
		set_difference -= col_thickness ;
		ASSERT (set_difference >= 0) ;
		/* [ -------------------------------- */
		/* modified for ccolamd from UMFPACK: */
		/* absorb this row if the set difference becomes zero */
		if (set_difference == 0 && aggressive)
		{
		    DEBUG3 (("aggressive absorption. Row: "ID"\n", row)) ;

                    if (Row [row].front != EMPTY)
                    {
                        /* Row [row].front is a child of current front. */
                        child = Row [row].front ;
                        Front_parent [child] = nfr ;
                        DEBUG1 (("Front "ID" => front "ID", aggressive\n",
                                    child, nfr)) ;
                    }
                    else
                    {
                        /* this is an original row.  Keep track of which front
                         * assembles it, for the row-merge tree */
                        InFront [row] = nfr ;
                        DEBUG1 (("Row "ID" => front "ID", aggressive\n",
                                    row, nfr)) ;
                    }

                    KILL_ROW (row) ;

                    /* sum the thicknesses of all the rows */
                    /* ASSERT (Row [row].thickness > 0) ; */
                    pivot_row_thickness += Row [row].thickness ;
                    Row [row].thickness = 0 ;
		/* -------------------------------- ] */

		}
		else
		{
		    /* save the new mark */
		    Row [row].shared2.mark = set_difference + tag_mark ;
		}
	    }
	}

#ifndef NDEBUG
	debug_deg_lists (n_row, n_col, Row, Col, head, min_score,
	cset_start[current_set+1]-(k+deadcol)-(cols_thickness),
		max_deg) ;
	cols_thickness = 0 ;
#endif /* NDEBUG */

	/* === Add up set differences for each column ======================= */

	DEBUG3 (("** Adding set differences phase. **\n")) ;

	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    /* get a column */
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    hash = 0 ;
	    cur_score = 0 ;
	    cp = &A [Col [col].start] ;
	    /* compact the column */
	    new_cp = cp ;
	    cp_end = cp + Col [col].length ;

	    DEBUG4 (("Adding set diffs for Col: "ID".\n", col)) ;

	    while (cp < cp_end)
	    {
		/* get a row */
		row = *cp++ ;
		ASSERT(row >= 0 && row < n_row) ;
		row_mark = Row [row].shared2.mark ;
		/* skip if dead */
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		   /* ---------------------------------- */
		   /* modified for ccolamd from UMFPACK: */
		    DEBUG4 ((" Row "ID", dead\n", row)) ;
		   /* ---------------------------------- */
		    continue ;
		}
		/* ---------------------------------- */
		/* modified for ccolamd from UMFPACK: */
		/* ASSERT (row_mark > tag_mark) ; */
		DEBUG4 ((" Row "ID", set diff "ID"\n", row, row_mark-tag_mark));
                ASSERT (row_mark >= tag_mark) ;
		/* ---------------------------------- */
		/* compact the column */
		*new_cp++ = row ;
		/* compute hash function */
		hash += row ;
		/* add set difference */
		cur_score += row_mark - tag_mark ;
		/* integer overflow... */
		cur_score = MIN (cur_score, n_col) ;
	    }

	    /* recompute the column's length */
	    Col [col].length = (Int) (new_cp - &A [Col [col].start]) ;

	    /* === Further mass elimination ================================= */

	    /*  modified for ccolamd */
	    if (Col [col].length == 0 && IN_CSET(col) == current_set)
	    {
		DEBUG4 (("further mass elimination. Col: "ID"\n", col)) ;
		/* nothing left but the pivot row in this column */
		KILL_PRINCIPAL_COL (col) ;
		pivot_row_degree -= Col [col].shared1.thickness ;
		ASSERT (pivot_row_degree >= 0) ;
		/* order it */
		Col [col].shared2.order = k ;
		/* increment order count by column thickness */
		k += Col [col].shared1.thickness ;

		/* ---------------------------------- */
		/* added for ccolamd from UMFPACK: */
                pivot_col_thickness += Col [col].shared1.thickness ;

                /* add to column list of front ... */
#ifndef NDEBUG
                DEBUG1 (("Mass")) ;
                dump_super (col, Col, n_col) ;
#endif
                Col [Col [col].lastcol].nextcol = Front_cols [nfr] ;
                Front_cols [nfr] = col ;
		/* ---------------------------------- */
	    }
	    else
	    {
		/* === Prepare for supercolumn detection ==================== */

		DEBUG4 (("Preparing supercol detection for Col: "ID".\n", col)) ;

		/* save score so far */
		Col [col].shared2.score = cur_score ;

		/* add column to hash table, for supercolumn detection */
		hash %= n_col + 1 ;

		DEBUG4 ((" Hash = "ID", n_col = "ID".\n", hash, n_col)) ;
		ASSERT (hash <= n_col) ;

		head_column = head [hash] ;
		if (head_column > EMPTY)
		{
		    /* degree list "hash" is non-empty, use prev (shared3) of */
		    /* first column in degree list as head of hash bucket */
		    first_col = Col [head_column].shared3.headhash ;
		    Col [head_column].shared3.headhash = col ;
		}
		else
		{
		    /* degree list "hash" is empty, use head as hash bucket */
		    first_col = - (head_column + 2) ;
		    head [hash] = - (col + 2) ;
		}
		Col [col].shared4.hash_next = first_col ;

		/* save hash function in Col [col].shared3.hash */
		Col [col].shared3.hash = (Int) hash ;
		ASSERT (COL_IS_ALIVE (col)) ;
	    }
	}

	/* The approximate external column degree is now computed.  */

	/* === Supercolumn detection ======================================== */

	DEBUG3 (("** Supercolumn detection phase. **\n")) ;

	detect_super_cols (

#ifndef NDEBUG
		n_col, Row,
#endif /* NDEBUG */

		Col, A, head, pivot_row_start, pivot_row_length,

		/* ------------------- */
		/* [ added for ccolamd */
		in_cset) ;
		/* ------------------- */

	/* === Kill the pivotal column ====================================== */

	DEBUG1 ((" KILLING column detect supercols "ID" \n", pivot_col ));
	KILL_PRINCIPAL_COL (pivot_col) ;

	/* ---------------------------------- */
	/* added for ccolamd from UMFPACK: */
	/* add columns to column list of front */
#ifndef NDEBUG
	DEBUG1 (("Pivot")) ;
	dump_super (pivot_col, Col, n_col) ;
#endif
	Col [Col [pivot_col].lastcol].nextcol = Front_cols [nfr] ;
	Front_cols [nfr] = pivot_col ;
	/* ---------------------------------- */

	/* === Clear mark =================================================== */

	tag_mark = clear_mark (tag_mark+max_deg+1, max_mark, n_row, Row) ;

#if 0
	tag_mark += (max_deg + 1) ;
	if (tag_mark >= max_mark)
	{
	    DEBUG2 (("clearing tag_mark\n")) ;
	    tag_mark = clear_mark (n_row, Row) ;
	}
#endif

#ifndef NDEBUG
	DEBUG3 (("check3\n")) ;
	debug_mark (n_row, Row, tag_mark, max_mark) ;
#endif /* NDEBUG */

	/* === Finalize the new pivot row, and column scores ================ */

	DEBUG3 (("** Finalize scores phase. **\n")) ;

	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	/* compact the pivot row */
	new_rp = rp ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    /* skip dead columns */
	    if (COL_IS_DEAD (col))
	    {
		continue ;
	    }
	    *new_rp++ = col ;
	    /* add new pivot row to column */
	    A [Col [col].start + (Col [col].length++)] = pivot_row ;

	    /* retrieve score so far and add on pivot row's degree. */
	    /* (we wait until here for this in case the pivot */
	    /* row's degree was reduced due to mass elimination). */
	    cur_score = Col [col].shared2.score + pivot_row_degree ;

	    /* calculate the max possible score as the number of */
	    /* external columns minus the 'k' value minus the */
	    /* columns thickness */
	    max_score = n_col - k - Col [col].shared1.thickness ;

	    /* make the score the external degree of the union-of-rows */
	    cur_score -= Col [col].shared1.thickness ;

	    /* make sure score is less or equal than the max score */
	    cur_score = MIN (cur_score, max_score) ;
	    ASSERT (cur_score >= 0) ;

	    /* store updated score */
	    Col [col].shared2.score = cur_score ;

	    /* === Place column back in degree list ========================= */

	    /* ---------------------- */
	    /* [ modified for ccolamd */
	    if ( IN_CSET(col) == current_set )
	    {
		ASSERT (min_score >= 0) ;
		ASSERT (min_score <= n_col) ;
		ASSERT (cur_score >= 0) ;
		ASSERT (cur_score <= n_col) ;
		ASSERT (head [cur_score] >= EMPTY) ;
		next_col = head [cur_score] ;
		Col [col].shared4.degree_next = next_col ;
		Col [col].shared3.prev = EMPTY ;
		if (next_col != EMPTY)
		{
		    Col [next_col].shared3.prev = col ;
		}
		head [cur_score] = col ;

		/* see if this score is less than current min */
		min_score = MIN (min_score, cur_score) ;
	    }
	    else
	    {
		Col [col].shared4.degree_next = EMPTY ;
		Col [col].shared3.prev = EMPTY ;
	    }
	    /* ---------------------- ] */

	}

#ifndef NDEBUG
	/* ---------------------- */
	/* [ modified for ccolamd */
	debug_deg_lists (n_row, n_col, Row, Col, head,
		min_score, cset_start[current_set+1]-(k+deadcol), max_deg) ;
	/* ---------------------- ] */
#endif /* NDEBUG */

	/* ---------------------------------- */
	/* [ added for ccolamd from UMFPACK: */
	/* frontal matrix can have more pivot cols than pivot rows for */
	/* singular matrices. */

	/* number of candidate pivot columns */
	Front_npivcol [nfr] = pivot_col_thickness ;

	/* all rows (not just size of contrib. block) */
	Front_nrows [nfr] = pivot_row_thickness ;

	/* all cols */
	Front_ncols [nfr] = pivot_col_thickness + pivot_row_degree ;

	Front_parent [nfr] = EMPTY ;

	pivot_row_thickness -= pivot_col_thickness ;
	DEBUG1 (("Front "ID" Pivot_row_thickness after pivot cols elim: "ID"\n",
	     nfr, pivot_row_thickness)) ;
	pivot_row_thickness = MAX (0, pivot_row_thickness) ;
	/* ---------------------------------- ] */

	/* === Resurrect the new pivot row ================================== */

	if ((pivot_row_degree > 0 && pivot_row_thickness > 0 && fact_type == 1)
	   || (pivot_row_degree > 0 && fact_type == 2))
	{
	    /* update pivot row length to reflect any cols that were killed */
	    /* during super-col detection and mass elimination */
	    Row [pivot_row].start  = pivot_row_start ;
	    Row [pivot_row].length = (Int) (new_rp - &A[pivot_row_start]) ;
	    Row [pivot_row].shared1.degree = pivot_row_degree ;
	    Row [pivot_row].shared2.mark = 0 ;
	    /* ---------------------------------- */
	    /*   added for ccolamd from UMFPACK: */
	    Row [pivot_row].thickness = pivot_row_thickness ;
	    Row [pivot_row].front = nfr ;
	    /* ---------------------------------- */
	    /* pivot row is no longer dead */
	    DEBUG1 (("Resurrect Pivot_row %d deg: %d\n",
			pivot_row, pivot_row_degree)) ;
	}

	/* ---------------------------------- */
	/* [ added for ccolamd from UMFPACK: */

#ifndef NDEBUG
	DEBUG1 (("Front "ID" : "ID" "ID" "ID" ", nfr,
		 Front_npivcol [nfr], Front_nrows [nfr], Front_ncols [nfr])) ;
	DEBUG1 ((" cols:[ ")) ;
	debug_d = 0 ;
	for (col = Front_cols [nfr] ; col != EMPTY ; col = Col [col].nextcol)
	{
		DEBUG1 ((" "ID, col)) ;
		ASSERT (col >= 0 && col < n_col) ;
		ASSERT (COL_IS_DEAD (col)) ;
		debug_d++ ;
		ASSERT (debug_d <= pivot_col_thickness) ;
	}
	ASSERT (debug_d == pivot_col_thickness) ;
	DEBUG1 ((" ]\n ")) ;
#endif
	 nfr++ ; /* one more front */
	/* ---------------------------------- ] */

    }


    /* === All principal columns have now been ordered ====================== */

   /* -------------------------------- */
   /*  added for ccolamd from UMFPACK: */
   *p_nfr = nfr ;
   /* -------------------------------- */

    return (ngarbage) ;
}


/* ========================================================================== */
/* === order_children deleted for cccolamd ================================== */
/* ========================================================================== */


/* ========================================================================== */
/* === detect_super_cols ==================================================== */
/* ========================================================================== */

/*
    Detects supercolumns by finding matches between columns in the hash buckets.
    Check amongst columns in the set A [row_start ... row_start + row_length-1].
    The columns under consideration are currently *not* in the degree lists,
    and have already been placed in the hash buckets.

    The hash bucket for columns whose hash function is equal to h is stored
    as follows:

	if head [h] is >= 0, then head [h] contains a degree list, so:

		head [h] is the first column in degree bucket h.
		Col [head [h]].headhash gives the first column in hash bucket h.

	otherwise, the degree list is empty, and:

		-(head [h] + 2) is the first column in hash bucket h.

    For a column c in a hash bucket, Col [c].shared3.prev is NOT a "previous
    column" pointer.  Col [c].shared3.hash is used instead as the hash number
    for that column.  The value of Col [c].shared4.hash_next is the next column
    in the same hash bucket.

    Assuming no, or "few" hash collisions, the time taken by this routine is
    linear in the sum of the sizes (lengths) of each column whose score has
    just been computed in the approximate degree computation.
    Not user-callable.
*/

PRIVATE void detect_super_cols
(
    /* === Parameters ======================================================= */

#ifndef NDEBUG
    /* these two parameters are only needed when debugging is enabled: */
    Int n_col,			/* number of columns of A */
    Colamd_Row Row [],		/* of size n_row+1 */
#endif /* NDEBUG */

    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* row indices of A */
    Int head [],		/* head of degree lists and hash buckets */
    Int row_start,		/* pointer to set of columns to check */
    Int row_length		/* number of columns to check */
    /* ----------------- */
    /* added for ccolamd */
    , Int in_cset[]		/* col -> cset mapping */
    /* ----------------- */
)
{
    /* === Local variables ================================================== */

    Int hash ;			/* hash value for a column */
    Int *rp ;			/* pointer to a row */
    Int c ;			/* a column index */
    Int super_c ;		/* column index of the column to absorb into */
    Int *cp1 ;			/* column pointer for column super_c */
    Int *cp2 ;			/* column pointer for column c */
    Int length ;		/* length of column super_c */
    Int prev_c ;		/* column preceding c in hash bucket */
    Int i ;			/* loop counter */
    Int *rp_end ;		/* pointer to the end of the row */
    Int col ;			/* a column index in the row to check */
    Int head_column ;		/* first column in hash bucket or degree list */
    Int first_col ;		/* first column in hash bucket */

    /* === Consider each column in the row ================================== */

    rp = &A [row_start] ;
    rp_end = rp + row_length ;
    while (rp < rp_end)
    {
	col = *rp++ ;
	if (COL_IS_DEAD (col))
	{
	    continue ;
	}

	/* get hash number for this column */
	hash = Col [col].shared3.hash ;
	ASSERT (hash <= n_col) ;

	/* === Get the first column in this hash bucket ===================== */

	head_column = head [hash] ;
	if (head_column > EMPTY)
	{
	    first_col = Col [head_column].shared3.headhash ;
	}
	else
	{
	    first_col = - (head_column + 2) ;
	}

	/* === Consider each column in the hash bucket ====================== */

	for (super_c = first_col ; super_c != EMPTY ;
	    super_c = Col [super_c].shared4.hash_next)
	{
	    ASSERT (COL_IS_ALIVE (super_c)) ;
	    ASSERT (Col [super_c].shared3.hash == hash) ;
	    length = Col [super_c].length ;

	    /* prev_c is the column preceding column c in the hash bucket */
	    prev_c = super_c ;

	    /* === Compare super_c with all columns after it ================ */

	    for (c = Col [super_c].shared4.hash_next ;
		 c != EMPTY ; c = Col [c].shared4.hash_next)
	    {
		ASSERT (c != super_c) ;
		ASSERT (COL_IS_ALIVE (c)) ;
		ASSERT (Col [c].shared3.hash == hash) ;

		/* not identical if lengths or scores are different, */
		/* or if in different constraint sets */
		if (Col [c].length != length ||
		    Col [c].shared2.score != Col [super_c].shared2.score
		    /* ---------------------- */
		    /* [ modified for ccolamd */
		    || IN_CSET(c)!= IN_CSET(super_c) )
		    /* ---------------------- */
		{
		    prev_c = c ;
		    continue ;
		}

		/* compare the two columns */
		cp1 = &A [Col [super_c].start] ;
		cp2 = &A [Col [c].start] ;

		for (i = 0 ; i < length ; i++)
		{
		    /* the columns are "clean" (no dead rows) */
		    ASSERT (ROW_IS_ALIVE (*cp1))  ;
		    ASSERT (ROW_IS_ALIVE (*cp2))  ;
		    /* row indices will same order for both supercols, */
		    /* no gather scatter nessasary */
		    if (*cp1++ != *cp2++)
		    {
			break ;
		    }
		}

		/* the two columns are different if the for-loop "broke" */
	        /* super columns should belong to the same constraint set */
		if (i != length)
		{
		    prev_c = c ;
		    continue ;
		}

		/* === Got it!  two columns are identical =================== */

		ASSERT (Col [c].shared2.score == Col [super_c].shared2.score) ;

		Col [super_c].shared1.thickness += Col [c].shared1.thickness ;
		Col [c].shared1.parent = super_c ;
		KILL_NON_PRINCIPAL_COL (c) ;
		/* order c later, in order_children() */
		Col [c].shared2.order = EMPTY ;
		/* remove c from hash bucket */
		Col [prev_c].shared4.hash_next = Col [c].shared4.hash_next ;

		/* --------------------------------- */
		/* [ added for ccolamd from UMFPACK: */
		/* add c to end of list of super_c */
		ASSERT (Col [super_c].lastcol >= 0) ;
		ASSERT (Col [super_c].lastcol < n_col) ;
		Col [Col [super_c].lastcol].nextcol = c ;
		Col [super_c].lastcol = Col [c].lastcol ;
#ifndef NDEBUG
		/* dump the supercolumn */
		DEBUG1 (("Super")) ;
		dump_super (super_c, Col, n_col) ;
#endif
		/* -------------------------------- ] */

	    }
	}

	/* === Empty this hash bucket ======================================= */

	if (head_column > EMPTY)
	{
	    /* corresponding degree list "hash" is not empty */
	    Col [head_column].shared3.headhash = EMPTY ;
	}
	else
	{
	    /* corresponding degree list "hash" is empty */
	    head [hash] = EMPTY ;
	}
    }
}


/* ========================================================================== */
/* === garbage_collection =================================================== */
/* ========================================================================== */

/*
    Defragments and compacts columns and rows in the workspace A.  Used when
    all avaliable memory has been used while performing row merging.  Returns
    the index of the first free position in A, after garbage collection.  The
    time taken by this routine is linear is the size of the array A, which is
    itself linear in the number of nonzeros in the input matrix.
    Not user-callable.
*/

PRIVATE Int garbage_collection  /* returns the new value of pfree */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows */
    Int n_col,			/* number of columns */
    Colamd_Row Row [],		/* row info */
    Colamd_Col Col [],		/* column info */
    Int A [],			/* A [0 ... Alen-1] holds the matrix */
    Int *pfree			/* &A [0] ... pfree is in use */
)
{
    /* === Local variables ================================================== */

    Int *psrc ;			/* source pointer */
    Int *pdest ;		/* destination pointer */
    Int j ;			/* counter */
    Int r ;			/* a row index */
    Int c ;			/* a column index */
    Int length ;		/* length of a row or column */

#ifndef NDEBUG
    Int debug_rows ;
    DEBUG2 (("Defrag..\n")) ;
    for (psrc = &A[0] ; psrc < pfree ; psrc++) ASSERT (*psrc >= 0) ;
    debug_rows = 0 ;
#endif /* NDEBUG */

    /* === Defragment the columns =========================================== */

    pdest = &A[0] ;
    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    psrc = &A [Col [c].start] ;

	    /* move and compact the column */
	    ASSERT (pdest <= psrc) ;
	    Col [c].start = (Int) (pdest - &A [0]) ;
	    length = Col [c].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		r = *psrc++ ;
		if (ROW_IS_ALIVE (r))
		{
		    *pdest++ = r ;
		}
	    }
	    Col [c].length = (Int) (pdest - &A [Col [c].start]) ;
	}
    }

    /* === Prepare to defragment the rows =================================== */

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_DEAD (r) || (Row [r].length == 0))
	{
	    /* This row is already dead, or is of zero length.  Cannot compact
	     * a row of zero length, so kill it.  NOTE: in the current version,
	     * there are no zero-length live rows.  Kill the row (for the first
	     * time, or again) just to be safe. */
	    KILL_ROW (r) ;
	}
	else
	{
	    /* save first column index in Row [r].shared2.first_column */
	    psrc = &A [Row [r].start] ;
	    Row [r].shared2.first_column = *psrc ;
	    ASSERT (ROW_IS_ALIVE (r)) ;
	    /* flag the start of the row with the one's complement of row */
	    *psrc = ONES_COMPLEMENT (r) ;
#ifndef NDEBUG
	    debug_rows++ ;
#endif /* NDEBUG */
	}
    }

    /* === Defragment the rows ============================================== */

    psrc = pdest ;
    while (psrc < pfree)
    {
	/* find a negative number ... the start of a row */
	if (*psrc++ < 0)
	{
	    psrc-- ;
	    /* get the row index */
	    r = ONES_COMPLEMENT (*psrc) ;
	    ASSERT (r >= 0 && r < n_row) ;
	    /* restore first column index */
	    *psrc = Row [r].shared2.first_column ;
	    ASSERT (ROW_IS_ALIVE (r)) ;

	    /* move and compact the row */
	    ASSERT (pdest <= psrc) ;
	    Row [r].start = (Int) (pdest - &A [0]) ;
	    length = Row [r].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		c = *psrc++ ;
		if (COL_IS_ALIVE (c))
		{
		    *pdest++ = c ;
		}
	    }
	    Row [r].length = (Int) (pdest - &A [Row [r].start]) ;

#ifndef NDEBUG
	    debug_rows-- ;
#endif /* NDEBUG */

	}
    }
    /* ensure we found all the rows */
    ASSERT (debug_rows == 0) ;

    /* === Return the new value of pfree ==================================== */

    return ((Int) (pdest - &A [0])) ;
}


/* ========================================================================== */
/* === clear_mark =========================================================== */
/* ========================================================================== */

/*
    Clears the Row [].shared2.mark array, and returns the new tag_mark.
    Return value is the new tag_mark.  Not user-callable.
*/

PRIVATE Int clear_mark	/* return the new value for tag_mark */
(
    /* === Parameters ======================================================= */

    Int tag_mark,	/* new value of tag_mark */
    Int max_mark,	/* max allowed value of tag_mark */

    Int n_row,		/* number of rows in A */
    Colamd_Row Row []	/* Row [0 ... n_row-1].shared2.mark is set to zero */
)
{
    /* === Local variables ================================================== */

    Int r ;

    if (tag_mark <= 0 || tag_mark >= max_mark)
    {
	for (r = 0 ; r < n_row ; r++)
	{
	    if (ROW_IS_ALIVE (r))
	    {
		Row [r].shared2.mark = 0 ;
	    }
	}
	tag_mark = 1 ;
    }

    return (tag_mark) ;
}


/* ========================================================================== */
/* === print_report ========================================================= */
/* ========================================================================== */

PRIVATE void print_report
(
    char *method,
    Int stats [CCOLAMD_STATS]
)
{

    Int i1, i2, i3 ;

    if (!stats)
    {
    	PRINTF ("%s: No statistics available.\n", method) ;
	return ;
    }

    i1 = stats [CCOLAMD_INFO1] ;
    i2 = stats [CCOLAMD_INFO2] ;
    i3 = stats [CCOLAMD_INFO3] ;

    if (stats [CCOLAMD_STATUS] >= 0)
    {
    	PRINTF ("%s: OK.  ", method) ;
    }
    else
    {
    	PRINTF ("%s: ERROR.  ", method) ;
    }

    switch (stats [CCOLAMD_STATUS])
    {

	case CCOLAMD_OK_BUT_JUMBLED:

	    PRINTF ("Matrix has unsorted or duplicate row indices.\n") ;

	    PRINTF ("%s: number of duplicate or out-of-order row indices: "ID"\n",
	    method, i3) ;

	    PRINTF ("%s: last seen duplicate or out-of-order row index:   "ID"\n",
	    method, INDEX (i2)) ;

	    PRINTF ("%s: last seen in column:                             "ID"",
	    method, INDEX (i1)) ;

	    /* no break - fall through to next case instead */

	case CCOLAMD_OK:

	    PRINTF ("\n") ;

 	    PRINTF ("%s: number of dense or empty rows ignored:           "ID"\n",
	    method, stats [CCOLAMD_DENSE_ROW]) ;

	    PRINTF ("%s: number of dense or empty columns ignored:        "ID"\n",
	    method, stats [CCOLAMD_DENSE_COL]) ;

	    PRINTF ("%s: number of garbage collections performed:         "ID"\n",
	    method, stats [CCOLAMD_DEFRAG_COUNT]) ;
	    break ;

	case CCOLAMD_ERROR_A_not_present:

	    PRINTF ("Array A (row indices of matrix) not present.\n") ;
	    break ;

	case CCOLAMD_ERROR_p_not_present:

	    PRINTF ("Array p (column pointers for matrix) not present.\n") ;
	    break ;

	case CCOLAMD_ERROR_nrow_negative:

	    PRINTF ("Invalid number of rows ("ID").\n", i1) ;
	    break ;

	case CCOLAMD_ERROR_ncol_negative:

	    PRINTF ("Invalid number of columns ("ID").\n", i1) ;
	    break ;

	case CCOLAMD_ERROR_nnz_negative:

	    PRINTF ("Invalid number of nonzero entries ("ID").\n", i1) ;
	    break ;

	case CCOLAMD_ERROR_p0_nonzero:

	    PRINTF ("Invalid column pointer, p [0] = "ID", must be zero.\n", i1) ;
	    break ;

	case CCOLAMD_ERROR_A_too_small:

	    PRINTF ("Array A too small.\n") ;
	    PRINTF ("        Need Alen >= "ID", but given only Alen = "ID".\n",
	    i1, i2) ;
	    break ;

	case CCOLAMD_ERROR_col_length_negative:

	    PRINTF
	    ("Column "ID" has a negative number of nonzero entries ("ID").\n",
	    INDEX (i1), i2) ;
	    break ;

	case CCOLAMD_ERROR_row_index_out_of_bounds:

	    PRINTF
	    ("Row index (row "ID") out of bounds ("ID" to "ID") in column "ID".\n",
	    INDEX (i2), INDEX (0), INDEX (i3-1), INDEX (i1)) ;
	    break ;

	case CCOLAMD_ERROR_out_of_memory:

	    PRINTF ("Out of memory.\n") ;
	    break ;

	case CCOLAMD_ERROR_invalid_cset:

	    PRINTF ("Cset invalid\n") ;
	    break ;
    }
}




/* ========================================================================== */
/* === colamd debugging routines ============================================ */
/* ========================================================================== */

/* When debugging is disabled, the remainder of this file is ignored. */

#ifndef NDEBUG


/* ========================================================================== */
/* === debug_structures ===================================================== */
/* ========================================================================== */

/*
    At this point, all empty rows and columns are dead.  All live columns
    are "clean" (containing no dead rows) and simplicial (no supercolumns
    yet).  Rows may contain dead columns, but all live rows contain at
    least one live column.
*/

PRIVATE void debug_structures
(
    /* === Parameters ======================================================= */

    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A []
    /* ----------------- */
    /* Added for ccolamd */
    , Int in_cset []
    , Int cset_start []
    /* ----------------- */
)
{
    /* === Local variables ================================================== */

    Int i ;
    Int c ;
    Int *cp ;
    Int *cp_end ;
    Int len ;
    Int score ;
    Int r ;
    Int *rp ;
    Int *rp_end ;
    Int deg ;
    /* ----------------- */
    /* Added for ccolamd */
    Int cs ;
    /* ----------------- */

    /* === Check A, Row, and Col ============================================ */

    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    len = Col [c].length ;
	    score = Col [c].shared2.score ;
	    DEBUG4 (("initial live col %5d %5d %5d\n", c, len, score)) ;
	    ASSERT (len > 0) ;
	    ASSERT (score >= 0) ;
	    ASSERT (Col [c].shared1.thickness == 1) ;
	    cp = &A [Col [c].start] ;
	    cp_end = cp + len ;
	    while (cp < cp_end)
	    {
		r = *cp++ ;
		ASSERT (ROW_IS_ALIVE (r)) ;
	    }
	}
	else
	{
	    /* -------------------- */
	    /* Modified for ccolamd */
	    i = Col [c].shared2.order ;
	    cs = IN_CSET(c) ;
	    ASSERT (i >= cset_start[cs] && i < cset_start[cs+1]) ;
	    /* -------------------- */
	}
    }

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_ALIVE (r))
	{
	    i = 0 ;
	    len = Row [r].length ;
	    deg = Row [r].shared1.degree ;
	    ASSERT (len > 0) ;
	    ASSERT (deg > 0) ;
	    rp = &A [Row [r].start] ;
	    rp_end = rp + len ;
	    while (rp < rp_end)
	    {
		c = *rp++ ;
		if (COL_IS_ALIVE (c))
		{
		    i++ ;
		}
	    }
	    ASSERT (i > 0) ;
	}
    }
}


/* ========================================================================== */
/* === debug_deg_lists ====================================================== */
/* ========================================================================== */

/*
    Prints the contents of the degree lists.  Counts the number of columns
    in the degree list and compares it to the total it should have.  Also
    checks the row degrees.
*/

PRIVATE void debug_deg_lists
(
    /* === Parameters ======================================================= */

    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int head [],
    Int min_score,
    Int should,
    Int max_deg
)
{
    /* === Local variables ================================================== */

    Int deg ;
    Int col ;
    Int have ;
    Int row ;

    /* === Check the degree lists =========================================== */

    if (n_col > 10000 && colamd_debug <= 0)
    {
	return ;
    }
    have = 0 ;
    DEBUG4 (("Degree lists: "ID"\n", min_score)) ;
    for (deg = 0 ; deg <= n_col ; deg++)
    {
	col = head [deg] ;
	if (col == EMPTY)
	{
	    continue ;
	}
	/* DEBUG4 ((""ID": (prev "ID")", deg, Col [col].shared3.prev)) ; */
	ASSERT (Col [col].shared3.prev == EMPTY) ;
	while (col != EMPTY)
	{
	    DEBUG4 ((" "ID"", col)) ;
	    have += Col [col].shared1.thickness ;
	    ASSERT (COL_IS_ALIVE (col)) ;
	    col = Col [col].shared4.degree_next ;
	}
	DEBUG4 (("\n")) ;
    }
    DEBUG4 (("should "ID" have "ID"\n", should, have)) ;
    ASSERT (should == have) ;

    /* === Check the row degrees ============================================ */

    if (n_row > 10000 && colamd_debug <= 0)
    {
	return ;
    }
    for (row = 0 ; row < n_row ; row++)
    {
	if (ROW_IS_ALIVE (row))
	{
	    ASSERT (Row [row].shared1.degree <= max_deg) ;
	}
    }
}


/* ========================================================================== */
/* === debug_mark =========================================================== */
/* ========================================================================== */

/*
    Ensures that the tag_mark is less that the maximum and also ensures that
    each entry in the mark array is less than the tag mark.
*/

PRIVATE void debug_mark
(
    /* === Parameters ======================================================= */

    Int n_row,
    Colamd_Row Row [],
    Int tag_mark,
    Int max_mark
)
{
    /* === Local variables ================================================== */

    Int r ;

    /* === Check the Row marks ============================================== */

    ASSERT (tag_mark > 0 && tag_mark <= max_mark) ;
    if (n_row > 10000 && colamd_debug <= 0)
    {
	return ;
    }
    for (r = 0 ; r < n_row ; r++)
    {
	ASSERT (Row [r].shared2.mark < tag_mark) ;
    }
}


/* ========================================================================== */
/* === debug_matrix ========================================================= */
/* ========================================================================== */

/*
    Prints out the contents of the columns and the rows.
*/

PRIVATE void debug_matrix
(
    /* === Parameters ======================================================= */

    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A []
)
{
    /* === Local variables ================================================== */

    Int r ;
    Int c ;
    Int *rp ;
    Int *rp_end ;
    Int *cp ;
    Int *cp_end ;

    /* === Dump the rows and columns of the matrix ========================== */

    if (colamd_debug < 3)
    {
	return ;
    }
    DEBUG3 (("DUMP MATRIX:\n")) ;
    for (r = 0 ; r < n_row ; r++)
    {
	DEBUG3 (("Row "ID" alive? "ID"\n", r, ROW_IS_ALIVE (r))) ;
	if (ROW_IS_DEAD (r))
	{
	    continue ;
	}

	/* -------------------- */
	/* Modified for ccolamd */
	DEBUG3 (("start "ID" length "ID" degree "ID" thickness "ID"\n",
		Row [r].start, Row [r].length, Row [r].shared1.degree,
		Row [r].thickness
		)) ;
	/* -------------------- */

	rp = &A [Row [r].start] ;
	rp_end = rp + Row [r].length ;
	while (rp < rp_end)
	{
	    c = *rp++ ;
	    DEBUG4 (("	"ID" col "ID"\n", COL_IS_ALIVE (c), c)) ;
	}
    }

    for (c = 0 ; c < n_col ; c++)
    {
	DEBUG3 (("Col "ID" alive? "ID"\n", c, COL_IS_ALIVE (c))) ;
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	DEBUG3 (("start "ID" length "ID" shared1 "ID" shared2 "ID"\n",
		Col [c].start, Col [c].length,
		Col [c].shared1.thickness, Col [c].shared2.score)) ;
	cp = &A [Col [c].start] ;
	cp_end = cp + Col [c].length ;
	while (cp < cp_end)
	{
	    r = *cp++ ;
	    DEBUG4 (("	"ID" row "ID"\n", ROW_IS_ALIVE (r), r)) ;
	}
    }
}

/* ------------------------------------------ */
/* dump_super added for ccolamd from UMFPACK: */
PRIVATE void dump_super
(
    Int super_c,
    Colamd_Col Col [],
    Int n_col
)
{
    Int col, ncols ;

    DEBUG1 ((" =[ ")) ;
    ncols = 0 ;
    for (col = super_c ; col != EMPTY ; col = Col [col].nextcol)
    {
        DEBUG1 ((" "ID, col)) ;
        ASSERT (col >= 0 && col < n_col) ;
        if (col != super_c)
        {
            ASSERT (COL_IS_DEAD (col)) ;
        }
        if (Col [col].nextcol == EMPTY)
        {
            ASSERT (col == Col [super_c].lastcol) ;
        }
        ncols++ ;
        ASSERT (ncols <= Col [super_c].shared1.thickness) ;
    }
    ASSERT (ncols == Col [super_c].shared1.thickness) ;
    DEBUG1 (("]\n")) ;
}
/* ------------------------------------------ */


PRIVATE void colamd_get_debug
(
    char *method
)
{
    FILE *debug_file ;
    colamd_debug = 0 ;		/* no debug printing */

    /* -------------------- */
    /* Modified for ccolamd. */
    /* Read debug info from the debug file. */
    /* get the print level from the debug file. */
    debug_file = fopen("debug", "r") ;
    if (debug_file)
    {
	(void) fscanf (debug_file, ""ID"", &colamd_debug) ;
	(void) fclose (debug_file) ;
    }

    DEBUG0 (("%s: debug version, D = "ID" (THIS WILL BE SLOW!)\n",
    	method, colamd_debug)) ;
    DEBUG0 ((" Debug printing level: "ID"\n", colamd_debug)) ;
    /* -------------------- */
}

#endif /* NDEBUG */

