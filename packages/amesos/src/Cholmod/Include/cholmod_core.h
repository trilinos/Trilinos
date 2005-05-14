/* ========================================================================== */
/* === Include/cholmod_core.h =============================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* CHOLMOD Core module: basic CHOLMOD objects and routines.
 * Required by all CHOLMOD modules.  Requires no other module or package.
 *
 * The CHOLMOD modules are:
 *
 * Core		    basic data structures and definitions
 * Check	    check/print the 5 CHOLMOD objects, & 3 kinds of int vectors
 * Cholesky	    sparse Cholesky factorization
 * Modify	    sparse Cholesky update/downdate/row-add/row-delete
 * MatrixOps	    sparse matrix functions (add, multiply, norm, ...)
 * Supernodal	    supernodal sparse Cholesky factorization
 * Partition	    graph-partitioning based orderings
 *
 * The CHOLMOD objects:
 * --------------------
 *
 * cholmod_common   parameters, statistics, and workspace
 * cholmod_sparse   a sparse matrix in compressed column form
 * cholmod_factor   an LL' or LDL' factorization
 * cholmod_dense    a dense matrix
 * cholmod_triplet  a sparse matrix in "triplet" form
 *
 * The Core module described here defines the CHOLMOD data structures, and
 * basic operations on them.  To create and solve a sparse linear system Ax=b,
 * the user must create A and b, populate them with values, and then pass them
 * to the routines in the CHOLMOD Cholesky module.  There are two primary
 * methods for creating A: (1) allocate space for a column-oriented sparse
 * matrix and fill it with pattern and values, or (2) create a triplet form
 * matrix and convert it to a sparse matrix.  The latter option is simpler.
 *
 * The matrices b and x are typically dense matrices, but can also be sparse.
 * You can allocate and free them as dense matrices with the
 * cholmod_allocate_dense and cholmod_free_dense routines.
 *
 * The cholmod_factor object contains the symbolic and numeric LL' or LDL'
 * factorization of sparse symmetric matrix.  The matrix must be positive
 * definite for an LL' factorization.  It need only be symmetric and have well-
 * conditioned leading submatrices for it to have an LDL' factorization
 * (CHOLMOD does not pivot for numerical stability).  It is typically created
 * with the cholmod_factorize routine in the Cholesky module, but can also
 * be initialized to L=D=I in the Core module and then modified by the Modify
 * module.  It must be freed with cholmod_free_factor, defined below.
 *
 * The Core routines for each object are described below.  Each list is split
 * into two parts: the primary routines and secondary routines.
 *
 * ============================================================================
 * === cholmod_common =========================================================
 * ============================================================================
 *
 * The Common object contains control parameters, statistics, and 
 * You must call cholmod_start before calling any other CHOLMOD routine, and
 * must call cholmod_finish as your last call to CHOLMOD (with two exceptions;
 * you may call cholmod_print_common and cholmod_check_common in the Check
 * module after calling cholmod_finish).
 *
 * cholmod_start		first call to CHOLMOD
 * cholmod_finish		last call to CHOLMOD
 *
 * cholmod_defaults		restore default parameters
 * cholmod_maxrank		column dimension of real workspace in Common
 * cholmod_allocate_work	allocate workspace in Common
 * cholmod_free_work		free workspace in Common
 * cholmod_clear_flag		clear Flag workspace in Common
 * cholmod_dmin			for internal use in CHOLMOD only
 * 
 * ============================================================================
 * === cholmod_sparse =========================================================
 * ============================================================================
 *
 * A sparse matrix is held in compressed column form.  In the basic type
 * ("packed", which corresponds to a MATLAB sparse matrix), an n-by-n matrix
 * with nz entries is held in three arrays: p of size n+1, i of size nz, and x
 * of size nz.  Row indices of column j are held in i [p [j] ... p [j+1]-1] and
 * in the same locations in x.  There may be no duplicate entries in a column.
 * Row indices in each column may be sorted or unsorted (CHOLMOD keeps track).
 *
 * cholmod_allocate_sparse	allocate a sparse matrix
 * cholmod_free_sparse		free a sparse matrix
 *
 * cholmod_reallocate_sparse	change the size (# entries) of sparse matrix
 * cholmod_nnz			number of nonzeros in a sparse matrix
 * cholmod_speye		sparse identity matrix
 * cholmod_zero			sparse zero matrix
 * cholmod_transpose_unsym	transpose unsymmetric sparse matrix
 * cholmod_transpose_sym	transpose symmetric sparse matrix
 * cholmod_transpose		transpose sparse matrix
 * cholmod_sort			sort row indices in each column of sparse matrix
 * cholmod_band			C = tril (triu (A,k1), k2)
 * cholmod_band_inplace		A = tril (triu (A,k1), k2)
 * cholmod_aat			C = A*A'
 * cholmod_copy_sparse		C = A, create an exact copy of a sparse matrix
 * cholmod_copy			C = A, with possible change of symmetry
 * cholmod_copy_sym_to_unsym	C = A, convert from symmetric to unsymmetric
 *
 * cholmod_add			C = alpha*A + beta*B
 *
 * ============================================================================
 * === cholmod_factor =========================================================
 * ============================================================================
 *
 * The data structure for an LL' or LDL' factorization is too complex to
 * describe in one sentence.  This object can hold the symbolic analysis alone,
 * or in combination with a "simplicial" (similar to a sparse matrix) or
 * "supernodal" form of the numerical factorization.  Only the routine to free
 * a factor is primary, since a factor object is created by the factorization
 * routine (cholmod_factorize).  It must be freed with cholmod_free_factor.
 *
 * cholmod_free_factor		free a factor
 *
 * cholmod_allocate_factor	allocate a factor (LL' or LDL')
 * cholmod_reallocate_factor	change the # entries in a factor 
 * cholmod_change_ftype		change the type of factor (e.g., LDL' to LL')
 * cholmod_pack_factor		pack the columns of a factor
 * cholmod_reallocate_column	resize a single column of a factor
 * cholmod_factor_to_sparse	create a sparse matrix copy of a factor
 * cholmod_copy_factor		create a copy of a factor
 *
 * Note that there is no cholmod_sparse_to_factor routine to create a factor
 * as a copy of a sparse matrix.  It could be done, after a fashion, but a
 * lower triangular sparse matrix would not necessarily have a chordal graph,
 * which would break the many CHOLMOD routines that rely on this property.
 *
 * ============================================================================
 * === cholmod_dense ==========================================================
 * ============================================================================
 *
 * The solve routines and some of the MatrixOps and Modify routines use dense
 * matrices as inputs.  These are held in column-major order.  With a leading
 * dimension of d, the entry in row i and column j is held in x [i+j*d].
 *
 * cholmod_allocate_dense	allocate a dense matrix
 * cholmod_free_dense		free a dense matrix
 *
 * cholmod_sparse_to_dense	create a dense matrix copy of a sparse matrix
 * cholmod_dense_to_sparse	create a sparse matrix copy of a dense matrix
 * cholmod_copy_dense		create a copy of a dense matrix
 * cholmod_copy_dense2		copy a dense matrix (pre-allocated)
 *
 * ============================================================================
 * === cholmod_triplet ========================================================
 * ============================================================================
 *
 * A sparse matrix held in triplet form is the simplest one for a user to
 * create.  It consists of a list of nz entries in arbitrary order, held in
 * three arrays: i, j, and x, each of length nk.  The kth entry is in row i[k],
 * column j[k], with value x[k].  There may be duplicate values; if A(i,j)
 * appears more than once, its value is the sum of the entries with those row
 * and column indices.
 *
 * cholmod_allocate_triplet	allocate a triplet matrix
 * cholmod_free_triplet		free a triplet matrix
 *
 * cholmod_reallocate_triplet	change the # of entries in a triplet matrix
 * cholmod_sparse_to_triplet	create a triplet matrix copy of a sparse matrix
 * cholmod_triplet_to_sparse	create a sparse matrix copy of a triplet matrix
 * cholmod_copy_triplet		create a copy of a triplet matrix
 * 
 * ============================================================================
 * === other ==================================================================
 * ============================================================================
 *
 * cholmod_scalar		a real or complex scalar, used as parameters to
 *				CHOLMOD routines.  This is a primary object, but
 *				has no routines associated with it.
 *
 * cholmod_postorder		postorder a tree (a secondary routine)
 *
 * ============================================================================
 * === memory management ======================================================
 * ============================================================================
 *
 * cholmod_malloc		malloc wrapper
 * cholmod_calloc		calloc wrapper
 * cholmod_free			free wrapper
 * cholmod_realloc		realloc wrapper
 * cholmod_realloc_multiple	realloc wrapper for multiple objects
 *
 * ============================================================================
 * === Core CHOLMOD prototypes ================================================
 * ============================================================================
 *
 * All CHOLMOD routines (in all modules) use the following protocol for return
 * values:
 *
 * int			TRUE (1) if successful, or FALSE (0) otherwise.
 * long			a value >= 0 if successful, or -1 otherwise.
 * double		a value >= 0 if successful, or -1 otherwise.
 * size_t		a value > 0 if successful, or 0 otherwise.
 * void *		a non-NULL pointer to newly allocated memory if
 *			successful, or NULL otherwise.
 * cholmod_sparse *	a non-NULL pointer to a newly allocated matrix
 *			if successful, or NULL otherwise.
 * cholmod_factor *	a non-NULL pointer to a newly allocated factor
 *			if successful, or NULL otherwise.
 * cholmod_triplet *	a non-NULL pointer to a newly allocated triplet
 *			matrix if successful, or NULL otherwise.
 * cholmod_dense *	a non-NULL pointer to a newly allocated triplet
 *			matrix if successful, or NULL otherwise.
 *
 * The last parameter to all routines is always a pointer to the CHOLMOD
 * Common object.
 *
 * TRUE and FALSE are not defined here, since they may conflict with the user
 * program.  A routine that described here returning TRUE or FALSE returns 1
 * or 0, respectively.  Any TRUE/FALSE parameter is true if nonzero, false if
 * zero.
 *
 * With two exceptions, all void * parameters are either int * or long *,
 * depending on Common->itype.  The itype of Common and all matrix and factor
 * parameters must match.  You cannot mix different itypes in the same call.
 * You can mix and match itypes within one program, but you must use multiple
 * Common objects to do so.
 *
 * If the itype is CHOLMOD_INT or CHOLMOD_INTLONG, the void * parameters
 * are all of type int * (pointer to an int array of the stated size).
 * Otherwise, they are long * (pointer to to a long array of the stated size).
 *
 * The two exceptions to this rule are cholmod_free and cholmod_realloc, which
 * can take an arbitrary void * pointer.  The type that the void * pointer
 * points to is not defined, just like the ANSI free and realloc.
 */

#ifndef CHOLMOD_CORE_H
#define CHOLMOD_CORE_H

/* required for size_t definition */
#include <stddef.h>

/* ========================================================================== */
/* === CHOLMOD integer type ================================================= */
/* ========================================================================== */

/* itype defines the types of integer used: */

#define CHOLMOD_INT 0		/* all integer arrays are int */
#define CHOLMOD_INTLONG 1	/* most are int, some are long */
#define CHOLMOD_LONG 2		/* all integer arrays are long */

/* The itype of all parameters for all CHOLMOD routines must match.
 *
 * FUTURE WORK: only CHOLMOD_INT is currently supported.
 */


/* ========================================================================== */
/* === CHOLMOD numerical type =============================================== */
/* ========================================================================== */

/* xtype defines the kind of numerical values used */

#define CHOLMOD_PATTERN 0	/* pattern only, no numerical values */
#define CHOLMOD_REAL 1		/* a real matrix */
#define CHOLMOD_COMPLEX 2	/* a complex matrix (ANSI compatible) */
#define CHOLMOD_ZOMPLEX 3	/* a complex matrix (MATLAB compatible) */

/* The xtype of all parameters for all CHOLMOD routines must match.
 * 
 * CHOLMOD_PATTERN: x and z are ignored.
 * CHOLMOD_DOUBLE:  x is non-null of size nzmax, z is ignored.
 * CHOLMOD_COMPLEX: x is non-null of size 2*nzmax, z is ignored.
 * CHOLMOD_ZOMPLEX: x and z are non-null of size nzmax
 *
 * In the real case, z is ignored.  The kth entry in the matrix is x [k].
 * There are two methods for the complex case.  In the ANSI-compatible
 * CHOLMOD_COMPLEX case, the real and imaginary parts of the kth entry
 * are in x [2*k] and x [2*k+1], respectively.  z is ignored.  In the
 * MATLAB-compatible CHOLMOD_ZOMPLEX case , the real and imaginary
 * parts of the kth entry are in x [k] and z [k].
 *
 * FUTURE WORK: the complex (or zomplex) case is not supported yet.
 */


/* ========================================================================== */
/* === CHOLMOD floating-point type ========================================== */
/* ========================================================================== */

/* dtype defines what the numerical type is (double or float) */

#define CHOLMOD_DOUBLE 0	/* all numerical values are double */
#define CHOLMOD_FLOAT 1		/* all numerical values are float */

/* The dtype of all parameters for all CHOLMOD routines must match.
 *
 * FUTURE WORK: the float case is not supported yet.
 */

/* ========================================================================== */
/* === Core/cholmod_scalar ================================================== */
/* ========================================================================== */

/* A real or complex floating-point scalar.  This is used only for passing
 * single scalars as parameters to a CHOLMOD routine, so that any CHOLMOD
 * routine can handle any real or complex type with no change in the calling
 * sequence.  It is used in both the CHOLMOD_DOUBLE and CHOLMOD_FLOAT cases.
 * It is not used as an component of any matrix data structure.
 */

typedef struct cholmod_scalar_struct
{
    double x ;	    /* real part */
    double z ;	    /* imaginary part, ignored when the matrix is real */

} cholmod_scalar ;



/* ========================================================================== */
/* === Core/cholmod_common ================================================== */
/* ========================================================================== */

/* control parameters, workspace, and statistics */

#define CHOLMOD_MAXMETHODS 8	/* maximum number of different methods that
				 * cholmod_analyze can try. Must be >= 8. */

/* Common->status values.  zero means success, negative means a fatal error,
 * positive is a warning. */
#define CHOLMOD_OK 0			/* success */
#define CHOLMOD_NOT_INSTALLED (-1)	/* failure: method not installed */
#define CHOLMOD_OUT_OF_MEMORY (-2)	/* failure: out of memory */
#define CHOLMOD_TOO_LARGE (-3)		/* failure: int overflow occured */
#define CHOLMOD_INVALID (-4)		/* failure: invalid input */
#define CHOLMOD_NOT_POSDEF (1)		/* warning: matrix not pos. def. */

/* ordering method (also used for L->ordering) */
#define CHOLMOD_NATURAL 0	/* use natural ordering */
#define CHOLMOD_GIVEN 1		/* use given permutation */
#define CHOLMOD_AMD 2		/* use minimum degree (AMD or COLAMD) */
#define CHOLMOD_METIS 3		/* use METIS' nested dissection */
#define CHOLMOD_ND 4		/* use CHOLMOD's version of nested dissection:
    * node bisector applied recursively, followed by constrained minimum
    * degree (CSYMAMD or CCOLAMD) */

typedef struct cholmod_common_struct
{
    /* ---------------------------------------------------------------------- */
    /* parameters for symbolic/numeric factorization and update/downdate */
    /* ---------------------------------------------------------------------- */

    double dmin ;	/* Smallest absolute value of diagonal entries of D.
			 * for LDL' factorization and update/downdate/rowadd/
	* rowdel.  Entries in the range 0 to dmin are replaced with dmin.
	* Entries in the range -dmin to 0 are replaced with -dmin.  No changes
	* are made to D if dmin <= 0.  Default: zero */

    double grow0 ;	/* For a dynamic LDL' factorization, L->i and L->x can
			 * grow if necessary.  grow0 is the factor by which
	* it grows.  For the initial space, L is of size MAX (1,grow0) times
	* the required space.  If L runs out of space, the new size of L is
	* MAX(1.2,grow0) times the new required space.   If you do not plan on
	* modifying the LDL' factorization in the Modify module, set grow0 to
	* zero.
	*
	* Default: 1.2 */

    double grow1 ;

    size_t grow2 ;	/* For a dynamic LDL' factorization, each column j of L
			 * is initialized with space equal to
	* grow1*L->ColCount[j] + grow2.  If grow1 is < 1 or grow2 == 0, then the
	* space allocated is exactly equal to L->ColCount[j].  If the column j
	* runs out of space, it increases to grow1*need + grow2 in size, where
	* need is the total # of nonzeros in that column.  If you do not plan on
	* modifying the LDL' factorization in the Modify module, set grow2 to
	* zero.
	*
	* Default: grow1 = 1.2, grow2 = 5. */

    size_t maxrank ;	/* rank of maximum update/downdate.  Valid values:
			 * 2, 4, or 8.  A value < 2 is set to 2, and a
	* value > 8 is set to 8.  It is then rounded up to the next highest
	* power of 2, if not already a power of 2.  Workspace (Xwork, below) of
	* size nrow-by-maxrank double's is allocated for the update/downdate.
	* If an update/downdate of rank-k is requested, with rank > maxrank,
	* it is done in steps of maxrank.  Default: 2, which saves memory but
	* can be slower than the maximum allowed value of 8.
	*
	* maxrank is also used to determine how the workspace in cholmod_solve
	* should be allocated.  If the solver requires space <= nrow*maxrank
	* doubles then it is allocated in Common->Xwork and left there for a
	* subsequent CHOLMOD operation.  If it is larger, then temporary space
	* is allocated instead (by cholmod_malloc) and freed when the solve is
	* done.  See cholmod_solve for a description of how much space it uses.
	*
	* Finally, maxrank is also used to determine how to allocate temporary
	* space in the supernodal factorization.  If L->maxcsize < maxrank*n*
	* sizeof(double) then Xwork is used.  Otherwise, temporary space is
	* allocated, and then freed after the factorization is computed.  See
	* cholmod_super_numeric for a description of how much space it uses.
	*/

    int supernodal ;	/* If TRUE, then cholmod_analyze performs a supernodal
			 * analysis for a subsequent supernodal factorization.
			 * Otherwise, do simplicial analysis.  Default: TRUE. */

    int final_ftype ;	/* If equal to CHOLMOD_LDL_PACKED, CHOLMOD_LDL_UNPACKED,
			 * CHOLMOD_LDL_DYNAMIC, or CHOLMOD_LL_PACKED, then
	* cholmod_factorize converts the factorization into this form when done.
	* If negative, or equal to CHOLMOD_AS_IS, then the factorization is left
	* as is (CHOLMOD_LL_SUPER if supernodal is TRUE, CHOLMOD_LDL_DYNAMIC
	* otherwise).  Default: CHOLMOD_AS_IS. */

    int final_resymbol ;/* if cholmod_factorize performed a supernodal
			 * factorization, final_resymbol is true, and
	* final_ftype is a simplicial numeric factorization, then numerically
	* zero entries that resulted from relaxed supernodal amalgamation are
	* removed.  This does not remove entries that are zero due to exact
	* numeric cancellation, since doing so would break the update/downdate
	* rowadd/rowdel routines.  Default: FALSE. */

    int blas_conform ;	/* TRUE if your BLAS conforms to the reference BLAS,
			 * in which dgemm and dsyrk do not access C on input
			 * when the BLAS beta input argument is zero. */

    int ieee ;		/* TRUE if (double) 0. is an all-zero bit field, which
		         * is true for IEEE floating-point. */

    /* supernodal relaxed amalgamation parameters: */
    size_t nrelax [3] ;

    double zrelax [3] ;

	/* Let ns be the total number of columns in two adjacent supernodes.
	 * Let z be the fraction of zero entries in the two supernodes if they
	 * are merged (z includes zero entries from prior amalgamations).  The
	 * two supernodes are merged if:
	 *
	 *    (ns <= nrelax [0]) || (no new zero entries added) ||
	 *    (ns <= nrelax [1] && z < zrelax [0]) ||
	 *    (ns <= nrelax [2] && z < zrelax [1]) || (z < zrelax [2])
	 *
	 * Default parameters result in the following rule:
	 *
	 *    (ns <= 4) || (no new zero entries added) ||
	 *    (ns <= 16 && z < 0.8) || (ns <= 48 && z < 0.1) || (z < 0.05)
	 */

    /* ---------------------------------------------------------------------- */
    /* printing and error handling options */
    /* ---------------------------------------------------------------------- */

    int print ;		/* print level. Default: 3 */
    int precise ;	/* if TRUE, print 16 digits.  Otherwise print 5 */
    void *file ;	/* output file (FILE *).  A void pointer is used to
			 * avoid including <stdio.h> in this file.
			 * Default: NULL, in which case stdout is used. */

    int try_catch ;	/* if TRUE, then ignore errors; CHOLMOD is in the middle
			 * of a try/catch block.  No error message is printed
			 * and the Common->error_handler function is not
			 * called. */

    void (*error_handler) (int status, char *message) ;

	/* Common->error_handler is the user's error handling routine.  If not
	 * NULL, this routine is called if an error occurs in CHOLMOD. */

    /* ---------------------------------------------------------------------- */
    /* ordering options */
    /* ---------------------------------------------------------------------- */

    /* The cholmod_analyze routine can try many different orderings and select
     * the best one.  It can also try one ordering method multiple times, with
     * different parameter settings.  The default is to use three orderings,
     * the user's permutation (if provided), AMD which is the fastest ordering
     * and generally gives good fill-in, and METIS.  CHOLMOD's nested dissection
     * (METIS with a constrained AMD) usually gives a better ordering than METIS
     * alone (by about 5% to 10%) but it takes more time.
     *
     * If you know the method that is best for your matrix, set Common->nmethods
     * to 1 and set Common->method [0] to the set of parameters for that method.
     * If you set it to 1 and do not provide a permutation, then only AMD will
     * be called.
     *
     * If METIS is not available, the default # of methods tried is 2 (the user
     * permutation, if any, and AMD).
     *
     * To try other methods, set Common->nmethods to the number of methods you
     * want to try.  The suite of default methods and their parameters is
     * described in the cholmod_start routine in cholmod_util.c.
     * You can modify the suite of methods you wish to try by modifying
     * Common.method [...] after calling cholmod_start.
     *
     * If you are going to factorize hundreds or more matrices with the same
     * nonzero pattern, you may wish to spend a great deal of time finding a
     * good permutation.  In this case, try setting Common->nmethods to 16.
     * The time spent in cholmod_analysis will be very high, but you need to
     * call it only once.   You should not do this if you do not have METIS
     * installed for use by CHOLMOD.
     *
     * cholmod_analyze sets Common->method to a value between 0 and nmethods-1.
     * Each ordering method uses the set of options defined by this parameter.
     */

    int nmethods ;	/* The number of ordering methods to try.  Default: 3,
			 * or 2 if METIS is not installed.  The first three
	* methods in the default suite of orderings is (1) use the given
	* permutation (if provided), (2) use AMD, and (3) use METIS.  Thus,
	* if you do not have METIS and do not provide a permutation, only
	* one method (AMD) will be used.  Maximum allowed value is
	* CHOLMOD_MAXMETHODS. */

    int current ;	/* The current method being tried.  Default: 0.  Valid
			 * range is 0 to nmethods-1. */

    int selected ;	/* The best method found. */

    /* The suite of ordering methods and parameters: */

    struct cholmod_method_struct
    {
	/* statistics for this method */
	double lnz ;	    /* nnz(L) excl. zeros from supernodal amalgamation*/

	double fl ;	    /* flop count for a simplicial LL' factorization */

	/* CHOLMOD nested dissection parameters: */
	double nd_prune ;   /* for pruning dense nodes prior to CHOLMOD's nested
			     * dissection.  Rows/columns with more than
	    * MAX (16, nd_prune * sqrt (n)) in an n-by-n matrix are removed
	    * prior to nested dissection.  They appear at the end of the
	    * reordered matrix.  Default is 10 (the same as AMD) */

	int ordering ;	    /* fill-reducing ordering to use */

	int nd_small ;	    /* do not partition graphs with fewer nodes than
			     * nd_small.  Default: 200 (same as METIS) */

	int nd_compress ;   /* If TRUE, compress the graph and subgraphs before
			     * partitioning them.  Default: TRUE */

	int nd_camd ;	    /* If TRUE, follow the nested dissection ordering
			     * with a constrained minimum degree ordering that
	    * respects the partitioning just found.  If you set nd_small very
	    * small, you may not need this ordering, and can save time by
	    * turning it off.  Default: TRUE. */

	/* FUTURE WORK:  AMD and COLAMD options could be added here. */

    } method [CHOLMOD_MAXMETHODS + 1] ;

    /* ---------------------------------------------------------------------- */
    /* memory management routines */
    /* ---------------------------------------------------------------------- */

    void *(*malloc_memory) (size_t) ;		/* pointer to malloc */
    void *(*realloc_memory) (void *, size_t) ;  /* pointer to realloc */
    void (*free_memory) (void *) ;		/* pointer to free */
    void *(*calloc_memory) (size_t, size_t) ;	/* pointer to calloc */

    /* ---------------------------------------------------------------------- */
    /* METIS workarounds */
    /* ---------------------------------------------------------------------- */

    double metis_memory ;   /* This is a parameter for CHOLMOD's interface to
			     * METIS, not a parameter to METIS itself.  METIS
	* uses an amount of memory that is difficult to estimate precisely
	* beforehand.  If it runs out of memory, it terminates your program.
	* All routines in CHOLMOD except for CHOLMOD's interface to METIS
	* return an error status and safely return to your program if they run
	* out of memory.  To mitigate this problem, the CHOLMOD interface
	* allocates a single block of memory equal in size to an empirical
	* upper bound of METIS's memory usage times the Common->metis_memory
	* parameter, and then immediately frees it.  It then calls METIS.  If
	* this pre-allocation fails, it is possible that METIS will fail as
	* well, and so CHOLMOD returns with an out-of-memory condition without
	* calling METIS.
	*
	* METIS_NodeND (used in the CHOLMOD_METIS ordering option) with its
	* default parameter settings typically uses about (4*nz+40n+4096)
	* times sizeof(int) memory, where nz is equal to the number of entries
	* in A for the symmetric case or AA' if an unsymmetric matrix is
	* being ordered (where nz includes both the upper and lower parts
	* of A or AA').  The observed "upper bound" (with 2 exceptions),
	* measured in an instrumented copy of METIS 4.0.1 on thousands of
	* matrices, is (10*nz+50*n+4096) * sizeof(int).  Two large matrices
	* exceeded this bound, one by almost a factor of 2 (Gupta/gupta2).
	*
	* To account for possible memory lost to fragmentation, the default for
	* metis_memory is 2.0.  If you don't mind your program being terminated
	* by METIS, set metis_memory to zero.
	*
	* If a matrix exceeds this predicted memory usage, AMD is attempted
	* instead.  It, too, may run out of memory, but if it does so it will
	* not terminate your program.
	*/

    double metis_dswitch ;	/* METIS_NodeND in METIS 4.0.1 gives a seg */
    size_t metis_nswitch ;	/* fault with one matrix of order n = 3005 and
				 * nz = 6,036,025.  This is a very dense graph.
     * The workaround is to use AMD instead of METIS for matrices of dimension
     * greater than Common->metis_nswitch (default 3000) or more and with
     * density of Common->metis_dswitch (default 0.66) or more.
     * cholmod_nested_dissection has no problems with the same matrix, even
     * though it uses METIS_NodeComputeSeparator on this matrix.  If this
     * seg fault does not affect you, set metis_nswitch to zero or less,
     * and CHOLMOD will not switch to AMD based just on the density of the
     * matrix (it will still switch to AMD if the metis_memory parameter
     * causes the switch).
     */

    /* ---------------------------------------------------------------------- */
    /* workspace */
    /* ---------------------------------------------------------------------- */

    /* CHOLMOD has several routines that take less time than the size of
     * workspace they require.  Allocating and initializing the workspace would
     * dominate the run time, unless workspace is allocated and initialized
     * just once.  CHOLMOD allocates this space when needed, and holds it here
     * between calls to CHOLMOD.  cholmod_start sets these pointers to NULL
     * (which is why it must be the first routine called in CHOLMOD).
     * cholmod_finish frees the workspace (which is why it must be the last
     * call to CHOLMOD).
     */

    size_t nrow ;	/* size of Flag and Head */
    long mark ;		/* mark value for Flag array */
    size_t iworksize ;	/* size of Iwork.  Upper bound: 2*nrow+ncol */
    size_t xworksize ;	/* size of Xwork, in bytes.  Upper bound: the larger of
			 * maxrank*nrow*sizeof(double) and 4*nrow*sizeof(int)
			 */

    /* initialized workspace: contents needed between calls to CHOLMOD */
    void *Flag ;	/* size nrow.  If cleared, Flag [i] < mark for all i */

    void *Head ;	/* size nrow+1. If cleared Head [i] = EMPTY for all i */

    void *Xwork ; 	/* Xwork is size nrow*sizeof(double) for cholmod_rowfac
			 * size 2*nrow*sizeof(double) for cholmod_rowadd,
	* cholmod_rowdel; size nrow*maxrank*sizeof(double) for cholmod_updown.
	* Normally cleared, with Xwork [i] = 0 for all i. 
	* FUTURE WORK: Xwork will be twice the size for the complex case.
	* It type is determined by Common->xtype and Common->dtype. */

    /* uninitialized workspace, contents not needed between calls to CHOLMOD */
    void *Iwork ;	/* size iworksize, typically 2*nrow+ncol */

    int itype ;		/* If CHOLMOD_LONG, Flag, Head, and Iwork are long.
			 * Otherwise all three arrays are int. */
    int xtype ;		/* real, complex, or zomplex */
    int dtype ;		/* double or float */

	/* Common->itype, xtype, and dtype is used to define the types of all
	 * sparse matrices, triplet matrices, dense matrices, and factors
	 * created using this Common struct.  The types of all parameters to all
	 * CHOLMOD routines must match.
	 *
	 * FUTURE WORK: only itype = CHOLMOD_INT, xtype = CHOLMOD_REAL
	 * and dtype = CHOLMOMD_DOUBLE are currently supported. */

    /* ---------------------------------------------------------------------- */
    /* statistics */
    /* ---------------------------------------------------------------------- */

    /* fl and lnz are set only in cholmod_analyze and cholmod_rowcolcounts (and
     * thus in the cholmod_factorize routine, all in the Cholesky module.
     * modfl is set only in the Modify module. */

    int status ;	    /* error code */
    double fl ;		    /* Cholesky flop count from most recent analysis */
    double lnz ;	    /* fundamental nz in L */
    double modfl ;	    /* flop count from most recent update/downdate/
			     * rowadd/rowdel (excluding flops to modify the 
			     * solution to Lx=b, if computed) */
    int malloc_count ;	    /* # of objects malloc'ed minus the # free'ed*/
    size_t memory_usage ;   /* peak memory usage in bytes */
    size_t memory_inuse ;   /* current memory usage in bytes */

} cholmod_common ;

/* -------------------------------------------------------------------------- */
/* cholmod_start */
/* -------------------------------------------------------------------------- */

int cholmod_start
(
    int itype,
    int xtype,
    int dtype,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_finish */
/* -------------------------------------------------------------------------- */

int cholmod_finish
(
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_defaults */
/* -------------------------------------------------------------------------- */

int cholmod_defaults
(
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_maxrank */
/* -------------------------------------------------------------------------- */

size_t cholmod_maxrank
(
    size_t n,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_allocate_work */
/* -------------------------------------------------------------------------- */

long cholmod_allocate_work
(
    size_t nrow,	/* A matrix will be nrow-by-ncol */
    size_t iworksize,	/* Iwork size in sizeof(int)'s, typically 2*nrow+ncol */
    size_t xwork,	/* number of items in Xwork */
    size_t xunits,	/* size of each item in Xwork */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_free_work */
/* -------------------------------------------------------------------------- */

int cholmod_free_work
(
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_clear_flag */
/* -------------------------------------------------------------------------- */

long cholmod_clear_flag
(
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_error */
/* -------------------------------------------------------------------------- */

int cholmod_error
(
    int status,
    char *msg,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_dmin */
/* -------------------------------------------------------------------------- */

double cholmod_dmin
(
    double dj,
    cholmod_common *Common
) ;


/* ========================================================================== */
/* === Core/cholmod_sparse ================================================== */
/* ========================================================================== */

/* A sparse matrix stored in compressed-column form. */

typedef struct cholmod_sparse_struct
{
    size_t nrow ;	/* the matrix is nrow-by-ncol */
    size_t ncol ;
    size_t nzmax ;	/* maximum number of entries in the matrix */

    void *p ;		/* p [0..ncol], the column pointers */
    void *i ;		/* i [0..nzmax-1], the row indices */

    /* for unpacked matrices only: */
    void *nz ;		/* nz [0..ncol-1], the # of nonzeros in each col.  In
			 * packed form, the nonzero pattern of column j is in
	* A->i [A->p [j] ... A->p [j+1]-1].  In unpacked form, column j is in
	* A->i [A->p [j] ... A->p [j]+A->nz[j]-1] instead.  In both cases, the
	* numerical values (if present) are in the corresponding locations in
	* the array x (or z if A->xtype is CHOLMOD_ZOMPLEX). */

    /* double or float: */
    void *x ;		/* size nzmax or 2*nzmax, if present */
    void *z ;		/* size nzmax, if present */

    int stype ;		/* Describes what parts of the matrix are considered:
			 *
	* 0:  matrix is "unsymmetric": use both upper and lower triangular parts
	*     (the matrix may actually be symmetric in pattern and value, but
	*     both parts are explicitly stored and used).  May be square or
	*     rectangular.
	* >0: matrix is square and symmetric, use upper triangular part.
	*     Entries in the lower triangular part are ignored.
	* <0: matrix is square and symmetric, use lower triangular part.
	*     Entries in the upper triangular part are ignored.
	*
	* Note that stype>0 and stype<0 are different for cholmod_sparse and
	* cholmod_triplet.  See the cholmod_triplet data structure for more
	* details.
	*/

    int itype ;		/* CHOLMOD_INT:     p, i, and nz are int.
			 * CHOLMOD_INTLONG: p is long, i and nz are int.
			 * CHOLMOD_LONG:    p, i, and nz are long.  */

    int xtype ;		/* pattern, real, complex, or zomplex */
    int dtype ;		/* x and z are double or float */
    int sorted ;	/* TRUE if columns are sorted, FALSE otherwise */
    int packed ;	/* TRUE if packed (nz must present), FALSE otherwise
			 * (nz is ignored) */

} cholmod_sparse ;

/* -------------------------------------------------------------------------- */
/* cholmod_allocate_sparse */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_allocate_sparse
(
    size_t nrow,    /* # of rows of A */
    size_t ncol,    /* # of columns of A */
    size_t nzmax,   /* max # of nonzeros of A */
    int sorted,	    /* TRUE if columns of A will be sorted, FALSE otherwise */
    int packed,	    /* TRUE if A will be packed, FALSE otherwise */
    int symmetry,   /* symmetry of A */
    int values,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_free_sparse */
/* -------------------------------------------------------------------------- */

int cholmod_free_sparse
(
    cholmod_sparse **A,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_reallocate_sparse */
/* -------------------------------------------------------------------------- */

int cholmod_reallocate_sparse
(
    cholmod_sparse *A,
    size_t newsize,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_nnz */
/* -------------------------------------------------------------------------- */

long cholmod_nnz
(
    cholmod_sparse *A,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_speye */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_speye
(
    size_t nrow,
    size_t ncol,
    size_t nzmax,
    int values,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_zero */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_zero
(
    size_t nrow,
    size_t ncol,
    size_t nzmax,
    int values,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_transpose_unsym */
/* -------------------------------------------------------------------------- */

int cholmod_transpose_unsym
(
    cholmod_sparse *A,	/* nrow-by-ncol, stored in column form */
    int values,		/* TRUE if values transposed, FALSE otherwise */
    void *Perm,		/* size nrow, if present (can be NULL) */
    void *fset,
    size_t fsize,
    cholmod_sparse *F,	/* A (:,f) stored in row form */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_transpose_sym */
/* -------------------------------------------------------------------------- */

int cholmod_transpose_sym
(
    cholmod_sparse *A,
    int values,		/* TRUE if values transposed, FALSE otherwise */
    void *Perm,		/* size nrow, if present (can be NULL) */
    cholmod_sparse *F,	/* A (Perm,Perm) stored in row form */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_transpose */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_transpose
(
    /* inputs, not modified: */
    cholmod_sparse *A,	/* nrow-by-ncol */
    int values,		/* TRUE if values transposed, FALSE otherwise */
    void *Perm,		/* if non-NULL, F = A(p,f) or A(p,p) */
    void *fset,		/* if non-NULL, F = A(*,f) */
    size_t fsize,	/* fset and fsize are used for unsymmetric case only */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_sort */
/* -------------------------------------------------------------------------- */

int cholmod_sort
(
    cholmod_sparse *A,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_copy_sparse */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_copy_sparse
(
    cholmod_sparse *A,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_band */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_band
(
    cholmod_sparse *A,	/* input */
    long k1,		/* ignore entries below the k1-st diagonal */
    long k2,		/* ignore entries above the k2-nd diagonal */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_band_inplace */
/* -------------------------------------------------------------------------- */

int cholmod_band_inplace
(
    cholmod_sparse *A,	/* input/output */
    long k1,		/* ignore entries below the k1-st diagonal */
    long k2,		/* ignore entries above the k2-nd diagonal */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_copy */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_copy
(
    cholmod_sparse *A,	/* input */
    int stype,		/* requested stype of C */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_copy_sym_to_unsym */
/* -------------------------------------------------------------------------- */

/* Same as C = cholmod_copy (A, 0, mode, Common) where A->stype != 0, except
 * that C can be allocated with extra space. */

cholmod_sparse *cholmod_copy_sym_to_unsym
(
    /* inputs, not modified on output */
    cholmod_sparse *A,
    size_t extra,	/* extra space to allocate in C */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_aat */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_aat
(
    cholmod_sparse *A,	/* input */
    void *fset,
    size_t fsize,
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_add		    C = alpha*A + beta*B */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_add
(
    cholmod_sparse *A,
    cholmod_sparse *B,
    cholmod_scalar alpha,
    cholmod_scalar beta,
    int values,		    /* if TRUE compute the numerical values of C */
    int sorted,
    cholmod_common *Common
) ;


/* ========================================================================== */
/* === Core/cholmod_factor ================================================== */
/* ========================================================================== */

/* A symbolic and numeric factorization, either simplicial or supernodal.
 * In all cases, the row indices in the columns of L are kept sorted. */

/* symbolic factorizations (values not present): */
#define CHOLMOD_SYMBOLIC_SUPER (-2) /* supernodal */
#define CHOLMOD_SYMBOLIC (-1)	    /* simplicial (same for LDL', LL', ...) */

/* numeric factorizations (values present): */
#define CHOLMOD_AS_IS (-1)	/* leave factor as is */
#define CHOLMOD_LDL_PACKED 0	/* simplicial LDL' factorization, packed */
#define CHOLMOD_LDL_UNPACKED 1	/* simplicial LDL' factorization, unpacked */
#define CHOLMOD_LDL_DYNAMIC 2	/* simplicial LDL' factorization, dynamic */
#define CHOLMOD_LL_PACKED 3	/* simplicial LL' factorization, packed */
#define CHOLMOD_LL_SUPER 4	/* supernodal LL' factorization */

typedef struct cholmod_factor_struct
{
    /* ---------------------------------------------------------------------- */
    /* for both simplicial and supernodal factorizations */
    /* ---------------------------------------------------------------------- */

    size_t n ;		/* L is n-by-n */

    size_t minor ;	/* If the factorization failed, L->minor is the column
			 * at which it failed (in the range 0 to n-1).  A value
			 * of n means the factorization was successful or
			 * the matrix has not yet been factorized. */

    /* ---------------------------------------------------------------------- */
    /* symbolic ordering and analysis */
    /* ---------------------------------------------------------------------- */

    void *Perm ;	/* permutation used */
    void *ColCount ;	/* column counts for simplicial L */

    /* ---------------------------------------------------------------------- */
    /* simplicial factorization */
    /* ---------------------------------------------------------------------- */

    size_t nzmax ;	/* size of i and x */

    void *p ;		/* p [0..ncol], the column pointers */
    void *i ;		/* i [0..nzmax-1], the row indices */

    void *x ;		/* x [0..nzmax-1], the numerical values */
    void *z ;

    /* used for unpacked and dynamic factors only: */
    void *nz ;		/* nz [0..ncol-1], the # of nonzeros in each col.  In
			 * packed form, the nonzero pattern of column j is in
			 * i [p [j] ... p [j+1]-1].
			 * In unpacked form, the pattern of column j is in
			 * i [p [j] ... p [j]+nz[j]-1] instead.  In both cases,
			 * the numerical values (if present) are in the
			 * corresponding locations in the array x. */

    /* used for dynamic factors only: */
    void *next ;	/* size ncol+2. next [j] is the next column in i/x */
    void *prev ;	/* size ncol+2. prev [j] is the prior column in i/x.
			 * head of the list is ncol+1, and the tail is ncol.
			 * next and prev NULL if the matrix is not dynamic.
			 * A dynamic matrix is also unpacked, so nz must be
			 * present as well. */

    /* ---------------------------------------------------------------------- */
    /* supernodal factorization */
    /* ---------------------------------------------------------------------- */

    /* Note that L->x of size xsize is shared with simplicial data structure */

    size_t nsuper ;	/* number of supernodes */
    size_t ssize ;	/* size of s, integer part of supernodes */
    size_t xsize ;	/* size of x, real part of supernodes */
    size_t maxcsize ;	/* size of largest update matrix */
    size_t maxesize ;	/* max # of rows in supernodes, excl. triangular part */

    void *super ;	/* size nsuper+1, first col in each supernode */
    void *pi ;		/* size nsuper+1, pointers to integer patterns */
    void *px ;		/* size nsuper+1, pointers to real parts */
    void *s ;		/* size ssize, integer part of supernodes */

    /* ---------------------------------------------------------------------- */
    /* factorization type */
    /* ---------------------------------------------------------------------- */

    int ordering ;	/* ordering method used */

    int ftype ;		/* ftype < 0 denotes a symbolic factorization,
			 * ftype >= 0 denotes a numeric factorization.
	*
	* CHOLMOD_SYMBOLIC_SUPER (-2):
	*	a supernodal symbolic factorization.  Nothing is present except
	*	for L->Perm and L->ColCount.
	*
	* CHOLMOD_SYMBOLIC (-1):
	*	a simplicial symbolic factorization.  The simplicial symbolic
	*	information is present (Perm and ColCount), as is all of the
	*	supernodal factorization except for the numerica values (L->x).
	*
	* CHOLMOD_LDL_PACKED (0):
	*	a simplicial LDL' factorization in packed form.  The unit
	*	diagonal of L is not stored.  D is stored at the diagonal
	*	entry of L.  col j is in i/x [p[j]..p[j+1]-1].  
	*
	* CHOLMOD_LDL_UNPACKED (1):
	*	unpacked simplicial LDL' factorization.  col j is in
	*	i/x [p[j]..p[j]+nz[j]-1], columns are in order.  nz must be
	*	present.
	*
	* CHOLMOD_LDL_DYNAMIC (2):
	*	dynamic simplicial LDL' factorization.   col j is in i/x as
	*	unpacked case, columns need not be in order.  Columns appear in
	*	i/x in order given by next/prev list. 	nz, next, and prev
	*	must be present.
	*
	* CHOLMOD_LL_PACKED (3):
	*	a simplicial LL' factorization (packed only).
	*
	* CHOLMOD_LL_SUPER (4):
	*	a supernodal LL' factorization
	*
	* All numerical factors include all diagonal entries of L.  For all
	* simplicial numeric factors, i [p [k]] is k.
	*/

    int itype ;		/* CHOLMOD_INT:     all void arrays are int.
			 * CHOLMOD_INTLONG: p, pi, px are long, all others int.
			 * CHOLMOD_LONG:    all void arrays are long. */
    int xtype ;		/* pattern, real, complex, or zomplex */
    int dtype ;		/* x and z double or float */
    /* FUTURE WORK: only CHOLMOD_INT and CHOLMOD_DOUBLE are supported. */

} cholmod_factor ;

/* -------------------------------------------------------------------------- */
/* cholmod_allocate_factor */
/* -------------------------------------------------------------------------- */

cholmod_factor *cholmod_allocate_factor
(
    size_t n,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_free_factor */
/* -------------------------------------------------------------------------- */

int cholmod_free_factor
(
    cholmod_factor **L,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_reallocate_factor */
/* -------------------------------------------------------------------------- */

int cholmod_reallocate_factor
(
    cholmod_factor *L,
    size_t newsize,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_change_ftype */
/* -------------------------------------------------------------------------- */

int cholmod_change_ftype
(
    cholmod_factor *L,
    int factor,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_pack_factor */
/* -------------------------------------------------------------------------- */

int cholmod_pack_factor
(
    cholmod_factor *L,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_reallocate_column */
/* -------------------------------------------------------------------------- */

int cholmod_reallocate_column
(
    cholmod_factor *L,
    size_t j,		    /* the column to reallocate */
    size_t need,	    /* required size of column j */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_factor_to_sparse */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_factor_to_sparse
(
    cholmod_factor *L,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_copy_factor */
/* -------------------------------------------------------------------------- */

cholmod_factor *cholmod_copy_factor
(
    cholmod_factor *L,
    cholmod_common *Common
) ;


/* ========================================================================== */
/* === Core/cholmod_dense =================================================== */
/* ========================================================================== */

/* A dense matrix in column-oriented form.  It has no itype since it contains
 * no integers.  Entry in row i and column j is located in x [i+j*d].
 */

typedef struct cholmod_dense_struct
{
    size_t nrow ;	/* the matrix is nrow-by-ncol */
    size_t ncol ;
    size_t nzmax ;	/* maximum number of entries in the matrix */
    size_t d ;		/* leading dimension (d >= nrow must hold) */
    void *x ;		/* size nzmax or 2*nzmax, if present */
    void *z ;		/* size nzmax, if present */
    int xtype ;		/* pattern, real, complex, or zomplex */
    int dtype ;		/* x and z double or float */

} cholmod_dense ;

/*
 * FUTURE WORK: Ensure d is always of type size_t when used in X [i+j*d].
 * FUTURE WORK: only xtype of CHOLMOD_REAL is supported.
 */

/* -------------------------------------------------------------------------- */
/* cholmod_allocate_dense */
/* -------------------------------------------------------------------------- */

/* The values parameter determines how the space is initialized: */

#define CHOLMOD_NONE (-1)	/* do not initialize */
#define CHOLMOD_ZERO 0		/* set to zero */
#define CHOLMOD_ONE  1		/* set to one */
#define CHOLMOD_EYE  2		/* set to the identitiy matrix */

cholmod_dense *cholmod_allocate_dense
(
    size_t nrow,
    size_t ncol,
    size_t d,
    int values,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_free_dense */
/* -------------------------------------------------------------------------- */

int cholmod_free_dense
(
    cholmod_dense **X,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_sparse_to_dense */
/* -------------------------------------------------------------------------- */

cholmod_dense *cholmod_sparse_to_dense
(
    cholmod_sparse *A,
    int values,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_dense_to_sparse */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_dense_to_sparse
(
    cholmod_dense *X,
    int stype,
    int values,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_copy_dense */
/* -------------------------------------------------------------------------- */

cholmod_dense *cholmod_copy_dense
(
    cholmod_dense *X,
    size_t d,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_copy_dense2 */
/* -------------------------------------------------------------------------- */

int cholmod_copy_dense2
(
    cholmod_dense *X,
    cholmod_dense *Y,
    cholmod_common *Common
) ;


/* ========================================================================== */
/* === Core/cholmod_triplet ================================================= */
/* ========================================================================== */

/* A sparse matrix stored in triplet form. */

typedef struct cholmod_triplet_struct
{
    size_t nrow ;	/* the matrix is nrow-by-ncol */
    size_t ncol ;
    size_t nzmax ;	/* maximum number of entries in the matrix */
    size_t nnz ;	/* number of nonzeros in the matrix */

    void *i ;		/* i [0..nzmax-1], the row indices */
    void *j ;		/* j [0..nzmax-1], the column indices */
    void *x ;		/* size nzmax or 2*nzmax, if present */
    void *z ;		/* size nzmax, if present */

    int stype ;		/* Describes what parts of the matrix are considered:
			 *
	* 0:  matrix is "unsymmetric": use both upper and lower triangular parts
	*     (the matrix may actually be symmetric in pattern and value, but
	*     both parts are explicitly stored and used).  May be square or
	*     rectangular.
	* >0: matrix is square and symmetric.  Entries in the lower triangular
	*     part are transposed and added to the upper triangular part when
	*     the matrix is converted to cholmod_sparse form.
	* <0: matrix is square and symmetric.  Entries in the upper triangular
	*     part are transposed and added to the lower triangular part when
	*     the matrix is converted to cholmod_sparse form.
	*
	* Note that stype>0 and stype<0 are different for cholmod_sparse and
	* cholmod_triplet.  The reason is simple.  You can permute a symmetric
	* triplet matrix by simply replacing a row and column index with their
	* new row and column indices, via an inverse permutation.  Suppose
	* P = L->Perm is your permutation, and Pinv is an array of size n.
	* Suppose a symmetric matrix A is represent by a triplet matrix T, with
	* entries only in the upper triangular part.  Then the following code:
	*
	*	Ti = T->i ;
	*	Tj = T->j ;
	*	for (k = 0 ; k < n  ; k++) Pinv [P [k]] = k ;
	*	for (k = 0 ; k < nz ; k++) Ti [k] = Pinv [Ti [k]] ;
	*	for (k = 0 ; k < nz ; k++) Tj [k] = Pinv [Tj [k]] ;
	*
	* creates the triplet form of C=P*A*P'.  However, if T initially
	* contains just the upper triangular entries (T->stype = 1), after
	* permutation it has entries in both the upper and lower triangular
	* parts.  These entries should be transposed when constructing the
	* cholmod_sparse form of A, which is what cholmod_triplet_to_sparse
	* does.  Thus:
	*
	*	C = cholmod_triplet_to_sparse (T, &Common) ;
	*
	* will return the matrix C = P*A*P'.
	*
	* Since the triplet matrix T is so simple to generate, it's quite easy
	* to remove entries that you do not want, prior to converting T to the
	* cholmod_sparse form.  So if you include these entries in T, CHOLMOD
	* assumes that there must be a reason (such as the one above).  Thus,
	* no entry in a triplet matrix is ever ignored.
	*/

    int itype ;		/* CHOLMOD_LONG: i and j are long.  Otherwise int. */
    int xtype ;		/* pattern, real, complex, or zomplex */
    int dtype ;		/* x and z are double or float */

    /* FUTURE WORK: only itype of CHOLMOD_INT, xtype of CHOLMOD_PATTERN or
     * CHOLMOD_REAL, and dtype CHOLMOD_DOUBLE are currently supported. */

} cholmod_triplet ;


/* -------------------------------------------------------------------------- */
/* cholmod_allocate_triplet */
/* -------------------------------------------------------------------------- */

cholmod_triplet *cholmod_allocate_triplet
(
    size_t nrow,    /* # of rows of T */
    size_t ncol,    /* # of columns of T */
    size_t nzmax,   /* max # of nonzeros of T */
    int symmetry,   /* symmetry of T */
    int values,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_free_triplet */
/* -------------------------------------------------------------------------- */

int cholmod_free_triplet
(
    cholmod_triplet **T,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_free_triplet */
/* -------------------------------------------------------------------------- */

int cholmod_reallocate_triplet
(
    cholmod_triplet *T,
    size_t nznew,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_sparse_to_triplet */
/* -------------------------------------------------------------------------- */

cholmod_triplet *cholmod_sparse_to_triplet
(
    cholmod_sparse *A,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_triplet_to_sparse */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_triplet_to_sparse
(
    cholmod_triplet *T,
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_copy_triplet */
/* -------------------------------------------------------------------------- */

cholmod_triplet *cholmod_copy_triplet
(
    cholmod_triplet *T,
    cholmod_common *Common
) ;


/* ========================================================================== */
/* === Core/cholmod_postorder =============================================== */
/* ========================================================================== */

/* Compute the postorder of a tree.  Several CHOLMOD unrelated modules depend
 * on this routine, which is why it is in the Core module.
 */

long cholmod_postorder
(
    /* inputs, not modified on output: */
    void *Parent,	/* size nrow.  Parent [j] = p if p is the parent of j */
    size_t nrow,

    /* outputs, not defined on input */
    void *Post,		/* size nrow.  Post [k] = j if j is kth node in
			 * postordered tree */

    cholmod_common *Common
) ;


/* ========================================================================== */
/* === Core memory management =============================================== */
/* ========================================================================== */

/* The void * pointers here are truly void *.  For all other routines in
 * CHOLMOD, void * parameters are used to point to integer arrays of type
 * int or long, depending on Common->itype.
 *
 * The user may make use of these, just like malloc and free.  You can even
 * malloc an object and safely free it with cholmod_free, and visa versa
 * (except that the memory usage statistics will be corrupted).  These routines
 * do differ from malloc and free.  If cholmod_free is given a NULL pointer,
 * for example, it does nothing (unlike the ANSI free).  cholmod_realloc does
 * not return NULL if given a non-NULL pointer and a nonzero size, even if it
 * fails (it sets an error code in Common->status instead).
 *
 * CHOLMOD keeps track of the amount of memory it has allocated, and so the
 * cholmod_free routine also takes the size of the object being freed.  This
 * is only used for statistics.  If you, the user of CHOLMOD, pass the wrong
 * size, the only consequence is that the memory usage statistics will be
 * corrupted.
 */

void *cholmod_malloc
(
    size_t n,	    /* number of items */
    size_t size,    /* size of each item */
    cholmod_common *Common
) ;

void *cholmod_calloc
(
    size_t n,	    /* number of items */
    size_t size,    /* size of each item */
    cholmod_common *Common
) ;

void *cholmod_free
(
    void *p,
    size_t n,	    /* number of items */
    size_t size,    /* size of each item */
    cholmod_common *Common
) ;

void *cholmod_realloc
(
    void *p,
    size_t *n,	    /* on input, n = current size. on output, n = new size */
    size_t new,	    /* requested size of the block of memory */
    size_t size,    /* size of each item */
    cholmod_common *Common
) ;

int cholmod_realloc_multiple
(
    int nint,	    /* number of int/long blocks */
    int nreal,	    /* number of real blocks */
    void **I,	    /* int or long block */
    void **J,	    /* int or long block */
    void **X,	    /* complex, double, or float block */
    void **Z,	    /* double or float block */
    size_t *old,    /* current size of both blocks */
    size_t new,	    /* requested size */
    cholmod_common *Common
) ;

#endif
