#ifndef F2C_INCLUDE
#define F2C_INCLUDE
#define FORTRAN_STRLEN

#include "ml_lapack.h"
#include "string.h"
/* typedef long int integer; */
typedef int integer;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
/* typedef long long longint; */ /* system-dependent */

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnint;
#else
typedef long flag;
typedef long ftnint;
#endif

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

typedef long Long;	/* No longer used; formerly in Namelist */

struct Vardesc {	/* for Namelist */
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;


/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef NOTSTRICT_PROTO
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef real (*R_fp)();
typedef doublereal (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif
#endif
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;	/* complex function */
typedef VOID H_f;	/* character function */
typedef VOID Z_f;	/* double complex function */
typedef doublereal E_f;	/* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif
#endif
extern  double ml_d_sign(doublereal *a, doublereal *b);
extern  int ml_s_cat(char *lp, char *rpp[], int rnp[], int *np, ftnlen ll);
extern  int ml_s_copy(register char *a, register char *b, ftnlen la, ftnlen lb);
extern  double ml_pow_di(doublereal *ap, integer *bp);
extern  integer ml_s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb);

#ifndef ML_DTRSV_FUNC


/* Subroutine */ int MLFORTRAN(dtrsv)(char *uplo, char *trans, char *diag, integer *n, 
	doublereal *a, integer *lda, doublereal *x, integer *incx)
{


    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i, j;
    static integer ix, jx, kx;
    static logical nounit;


/*  Purpose   
    =======   

    DTRSV  solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix.   

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   A'*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular matrix and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular matrix and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             ML_max( 1, n ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*ML_abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! MLFORTRAN(lsame)(uplo, "U") && ! MLFORTRAN(lsame)(uplo, "L")) {
	info = 1;
    } else if (! MLFORTRAN(lsame)(trans, "N") && ! MLFORTRAN(lsame)(trans, "T") &&
	     ! MLFORTRAN(lsame)(trans, "C")) {
	info = 2;
    } else if (! MLFORTRAN(lsame)(diag, "U") && ! MLFORTRAN(lsame)(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < ML_max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	MLFORTRAN(xerbla)("DTRSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = MLFORTRAN(lsame)(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (MLFORTRAN(lsame)(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (MLFORTRAN(lsame)(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (X(j) != 0.) {
			if (nounit) {
			    X(j) /= A(j,j);
			}
			temp = X(j);
			for (i = j - 1; i >= 1; --i) {
			    X(i) -= temp * A(i,j);
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    if (X(jx) != 0.) {
			if (nounit) {
			    X(jx) /= A(j,j);
			}
			temp = X(jx);
			ix = jx;
			for (i = j - 1; i >= 1; --i) {
			    ix -= *incx;
			    X(ix) -= temp * A(i,j);
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = 1; j <= *n; ++j) {
		    if (X(j) != 0.) {
			if (nounit) {
			    X(j) /= A(j,j);
			}
			temp = X(j);
			for (i = j + 1; i <= *n; ++i) {
			    X(i) -= temp * A(i,j);
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		for (j = 1; j <= *n; ++j) {
		    if (X(jx) != 0.) {
			if (nounit) {
			    X(jx) /= A(j,j);
			}
			temp = X(jx);
			ix = jx;
			for (i = j + 1; i <= *n; ++i) {
			    ix += *incx;
			    X(ix) -= temp * A(i,j);
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (MLFORTRAN(lsame)(uplo, "U")) {
	    if (*incx == 1) {
		for (j = 1; j <= *n; ++j) {
		    temp = X(j);
		    for (i = 1; i <= j-1; ++i) {
			temp -= A(i,j) * X(i);
/* L90: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(j) = temp;
/* L100: */
		}
	    } else {
		jx = kx;
		for (j = 1; j <= *n; ++j) {
		    temp = X(jx);
		    ix = kx;
		    for (i = 1; i <= j-1; ++i) {
			temp -= A(i,j) * X(ix);
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(jx) = temp;
		    jx += *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = X(j);
		    for (i = *n; i >= j+1; --i) {
			temp -= A(i,j) * X(i);
/* L130: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(j) = temp;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = X(jx);
		    ix = kx;
		    for (i = *n; i >= j+1; --i) {
			temp -= A(i,j) * X(ix);
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= A(j,j);
		    }
		    X(jx) = temp;
		    jx -= *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRSV . */

} /* dtrsv_ */

#endif

#ifndef ML_DGETRF_FUNC

/* Subroutine */ int MLFORTRAN(dgetrf)(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGETRF computes an LU factorization of a general M-by-N matrix A   
    using partial pivoting with row interchanges.   

    The factorization has the form   
       A = P * L * U   
    where P is a permutation matrix, L is lower triangular with unit   
    diagonal elements (lower trapezoidal if m > n), and U is upper   
    triangular (upper trapezoidal if m < n).   

    This is the right-looking Level 3 BLAS version of the algorithm.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix to be factored.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,M).   

    IPIV    (output) INTEGER array, dimension (ML_min(M,N))   
            The pivot indices; for 1 <= i <= ML_min(M,N), row i of the   
            matrix was interchanged with row IPIV(i).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero. The factorization 
  
                  has been completed, but the factor U is exactly   
                  singular, and division by zero will occur if it is used 
  
                  to solve a system of equations.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static doublereal c_b16 = 1.;
    static doublereal c_b19 = -1.;
    
    /* System generated locals */
    integer   i__1, i__3, i__4, i__5;
    /* Local variables */
    static integer i, j;
    static integer iinfo;
    static integer jb, nb;


#define IPIV(I) ipiv[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ML_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DGETRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = MLFORTRAN(ilaenv)(&c__1, "DGETRF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
    if (nb <= 1 || nb >= ML_min(*m,*n)) {

/*        Use unblocked code. */

	MLFORTRAN(dgetf2)(m, n, &A(1,1), lda, &IPIV(1), info);
    } else {

/*        Use blocked code. */

	i__1 = ML_min(*m,*n);
	for (j = 1; nb < 0 ? j >= ML_min(*m,*n) : j <= ML_min(*m,*n); j += nb) {
/* Computing MIN */
	    i__3 = ML_min(*m,*n) - j + 1;
	    jb = ML_min(i__3,nb);

/*           Factor diagonal and subdiagonal blocks and test for e
xact   
             singularity. */

	    i__3 = *m - j + 1;
	    MLFORTRAN(dgetf2)(&i__3, &jb, &A(j,j), lda, &IPIV(j), &iinfo);

/*           Adjust INFO and the pivot indices. */

	    if (*info == 0 && iinfo > 0) {
		*info = iinfo + j - 1;
	    }
/* Computing MIN */
	    i__4 = *m, i__5 = j + jb - 1;
	    i__3 = ML_min(i__4,i__5);
	    for (i = j; i <= ML_min(*m,j+jb-1); ++i) {
		IPIV(i) = j - 1 + IPIV(i);
/* L10: */
	    }

/*           Apply interchanges to columns 1:J-1. */

	    i__3 = j - 1;
	    i__4 = j + jb - 1;
	    MLFORTRAN(dlaswp)(&i__3, &A(1,1), lda, &j, &i__4, &IPIV(1), &c__1);

	    if (j + jb <= *n) {

/*              Apply interchanges to columns J+JB:N. */

		i__3 = *n - j - jb + 1;
		i__4 = j + jb - 1;
		MLFORTRAN(dlaswp)(&i__3, &A(1,j+jb), lda, &j, &i__4, &
			IPIV(1), &c__1);

/*              Compute block row of U. */

		i__3 = *n - j - jb + 1;
		MLFORTRAN(dtrsm)("Left", "Lower", "No transpose", "Unit", &jb, &i__3, &
			c_b16, &A(j,j), lda, &A(j,j+jb), lda
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Lower"),strlen("No transpose"),strlen("Unit")
#endif
);
		if (j + jb <= *m) {

/*                 Update trailing submatrix. */

		    i__3 = *m - j - jb + 1;
		    i__4 = *n - j - jb + 1;
		    MLFORTRAN(dgemm)("No transpose", "No transpose", &i__3, &i__4, &jb, 
			    &c_b19, &A(j+jb,j), lda, &A(j,j+jb), lda, &c_b16, &A(j+jb,j+jb), lda
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("No transpose")
#endif
);
		}
	    }
/* L20: */
	}
    }
    return 0;

/*     End of DGETRF */

} /* dgetrf_ */
#endif

#ifndef ML_DGETRS_FUNC

/* Subroutine */ int MLFORTRAN(dgetrs)(char *trans, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, integer *info, int dummy1)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGETRS solves a system of linear equations   
       A * X = B  or  A' * X = B   
    with a general N-by-N matrix A using the LU factorization computed   
    by DGETRF.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A'* X = B  (Transpose)   
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The factors L and U from the factorization A = P*L*U   
            as computed by DGETRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,N).   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices from DGETRF; for 1<=i<=N, row i of the   
            matrix was interchanged with row IPIV(i).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= ML_max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b12 = 1.;
    static integer c_n1 = -1;
    
    /* System generated locals */
    integer     i__1;
    /* Local variables */
    static logical notran;



#define IPIV(I) ipiv[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    notran = MLFORTRAN(lsame)(trans, "N");
    if (! notran && ! MLFORTRAN(lsame)(trans, "T") && ! MLFORTRAN(lsame)(trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < ML_max(1,*n)) {
	*info = -5;
    } else if (*ldb < ML_max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DGETRS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    if (notran) {

/*        Solve A * X = B.   

          Apply row interchanges to the right hand sides. */

	MLFORTRAN(dlaswp)(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c__1);

/*        Solve L*X = B, overwriting B with X. */

	MLFORTRAN(dtrsm)("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Lower"),strlen("No transpose"),strlen("Unit")
#endif
);

/*        Solve U*X = B, overwriting B with X. */

	MLFORTRAN(dtrsm)("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b12, &
		A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Upper"),strlen("No transpose"),strlen("Non-unit")
#endif
);
    } else {

/*        Solve A' * X = B.   

          Solve U'*X = B, overwriting B with X. */

	MLFORTRAN(dtrsm)("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Upper"),strlen("Transpose"),strlen("Non-unit")
#endif
);

/*        Solve L'*X = B, overwriting B with X. */

	MLFORTRAN(dtrsm)("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Lower"),strlen("Transpose"),strlen("Unit")
#endif
);

/*        Apply row interchanges to the solution vectors. */

	MLFORTRAN(dlaswp)(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c_n1);
    }

    return 0;

/*     End of DGETRS */

} /* dgetrs_ */
#endif

#ifndef ML_XERBLA_FUNC

/* Subroutine */ int MLFORTRAN(xerbla)(char *srname, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    XERBLA  is an error handler for the LAPACK routines.   
    It is called by an LAPACK routine if an input parameter has an   
    invalid value.  A message is printed and execution stops.   

    Installers may consider modifying the STOP statement in order to   
    call system-specific exception-handling facilities.   

    Arguments   
    =========   

    SRNAME  (input) CHARACTER*6   
            The name of the routine which called XERBLA.   

    INFO    (input) INTEGER   
            The position of the invalid parameter in the parameter list   

            of the calling routine.   

   ===================================================================== 
*/

    printf("** On entry to %6s, parameter number %2ld had an illegal value\n",
		srname, (long int) *info);

/*     End of XERBLA */

    return 0;
} /* xerbla_ */
#endif

#ifndef ML_LSAME_FUNC
logical MLFORTRAN(lsame)(char *ca, char *cb)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    LSAME returns .TRUE. if CA is the same letter as CB regardless of   
    case.   

    Arguments   
    =========   

    CA      (input) CHARACTER*1   
    CB      (input) CHARACTER*1   
            CA and CB specify the single characters to be compared.   

   ===================================================================== 
  


       Test if the characters are equal */
    /* System generated locals */
    logical ret_val;
    /* Local variables */
    static integer inta, intb, zcode;


    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
       machines, on which ICHAR returns a value with bit 8 set.   
       ICHAR('A') on Prime machines returns 193 which is the same as   
       ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower o
r   
          upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
 or   
          upper case 'Z'. */

	if ((inta >= 129 && inta <= 137) || (inta >= 145 && inta <= 153) || (inta 
		>= 162 && inta <= 169)) {
	    inta += 64;
	}
	if ((intb >= 129 && intb <= 137) || (intb >= 145 && intb <= 153) || (intb 
		>= 162 && intb <= 169)) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
e   
          plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN   

       End of LSAME */

    return ret_val;
} /* lsame_ */
#endif

#ifndef ML_DTRSM_FUNC

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Subroutine */ int MLFORTRAN(dtrsm)(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, doublereal *alpha, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, int dummy1, int dummy2, int dummy3,
				      int dummy4)
{


    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i, j, k;
    static logical lside;
    static integer nrowa;
    static logical upper;
    static logical nounit;


/*  Purpose   
    =======   

    DTRSM  solves one of the matrix equations   

       op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,   

    where alpha is a scalar, X and B are m by n matrices, A is a unit, or 
  
    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of 
  

       op( A ) = A   or   op( A ) = A'.   

    The matrix X is overwritten on B.   

    Parameters   
    ==========   

    SIDE   - CHARACTER*1.   
             On entry, SIDE specifies whether op( A ) appears on the left 
  
             or right of X as follows:   

                SIDE = 'L' or 'l'   op( A )*X = alpha*B.   

                SIDE = 'R' or 'r'   X*op( A ) = alpha*B.   

             Unchanged on exit.   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix A is an upper or 
  
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in 
  
             the matrix multiplication as follows:   

                TRANSA = 'N' or 'n'   op( A ) = A.   

                TRANSA = 'T' or 't'   op( A ) = A'.   

                TRANSA = 'C' or 'c'   op( A ) = A'.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit triangular 
  
             as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of B. M must be at 
  
             least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of B.  N must be 
  
             at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry,  ALPHA specifies the scalar  alpha. When  alpha is 
  
             zero then  A is not referenced and  B need not be set before 
  
             entry.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m 
  
             when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. 
  
             Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k 
  
             upper triangular part of the array  A must contain the upper 
  
             triangular matrix  and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k 
  
             lower triangular part of the array  A must contain the lower 
  
             triangular matrix  and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u',  the diagonal elements of 
  
             A  are not referenced either,  but are assumed to be  unity. 
  
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then 
  
             LDA  must be at least  ML_max( 1, m ),  when  SIDE = 'R' or 'r' 
  
             then LDA must be at least ML_max( 1, n ).   
             Unchanged on exit.   

    B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
             Before entry,  the leading  m by n part of the array  B must 
  
             contain  the  right-hand  side  matrix  B,  and  on exit  is 
  
             overwritten by the solution matrix  X.   

    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared 
  
             in  the  calling  (sub)  program.   LDB  must  be  at  least 
  
             ML_max( 1, m ).   
             Unchanged on exit.   


    Level 3 Blas routine.   


    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    lside = MLFORTRAN(lsame)(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = MLFORTRAN(lsame)(diag, "N");
    upper = MLFORTRAN(lsame)(uplo, "U");

    info = 0;
    if (! lside && ! MLFORTRAN(lsame)(side, "R")) {
	info = 1;
    } else if (! upper && ! MLFORTRAN(lsame)(uplo, "L")) {
	info = 2;
    } else if (! MLFORTRAN(lsame)(transa, "N") && ! MLFORTRAN(lsame)(transa, "T") 
	    && ! MLFORTRAN(lsame)(transa, "C")) {
	info = 3;
    } else if (! MLFORTRAN(lsame)(diag, "U") && ! MLFORTRAN(lsame)(diag, "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < ML_max(1,nrowa)) {
	info = 9;
    } else if (*ldb < ML_max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	MLFORTRAN(xerbla)("DTRSM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	for (j = 1; j <= *n; ++j) {
	    for (i = 1; i <= *m; ++i) {
		B(i,j) = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (MLFORTRAN(lsame)(transa, "N")) {

/*           Form  B := alpha*inv( A )*B. */

	    if (upper) {
		for (j = 1; j <= *n; ++j) {
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = *alpha * B(i,j);
/* L30: */
			}
		    }
		    for (k = *m; k >= 1; --k) {
			if (B(k,j) != 0.) {
			    if (nounit) {
				B(k,j) /= A(k,k);
			    }
			    for (i = 1; i <= k-1; ++i) {
				B(i,j) -= B(k,j) * A(i,k);
/* L40: */
			    }
			}
/* L50: */
		    }
/* L60: */
		}
	    } else {
		for (j = 1; j <= *n; ++j) {
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = *alpha * B(i,j);
/* L70: */
			}
		    }
		    for (k = 1; k <= *m; ++k) {
			if (B(k,j) != 0.) {
			    if (nounit) {
				B(k,j) /= A(k,k);
			    }
			    for (i = k + 1; i <= *m; ++i) {
				B(i,j) -= B(k,j) * A(i,k);
/* L80: */
			    }
			}
/* L90: */
		    }
/* L100: */
		}
	    }
	} else {

/*           Form  B := alpha*inv( A' )*B. */

	    if (upper) {
		for (j = 1; j <= *n; ++j) {
		    for (i = 1; i <= *m; ++i) {
			temp = *alpha * B(i,j);
			for (k = 1; k <= i-1; ++k) {
			    temp -= A(k,i) * B(k,j);
/* L110: */
			}
			if (nounit) {
			    temp /= A(i,i);
			}
			B(i,j) = temp;
/* L120: */
		    }
/* L130: */
		}
	    } else {
		for (j = 1; j <= *n; ++j) {
		    for (i = *m; i >= 1; --i) {
			temp = *alpha * B(i,j);
			for (k = i + 1; k <= *m; ++k) {
			    temp -= A(k,i) * B(k,j);
/* L140: */
			}
			if (nounit) {
			    temp /= A(i,i);
			}
			B(i,j) = temp;
/* L150: */
		    }
/* L160: */
		}
	    }
	}
    } else {
	if (MLFORTRAN(lsame)(transa, "N")) {

/*           Form  B := alpha*B*inv( A ). */

	    if (upper) {
		for (j = 1; j <= *n; ++j) {
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = *alpha * B(i,j);
/* L170: */
			}
		    }
		    for (k = 1; k <= j-1; ++k) {
			if (A(k,j) != 0.) {
			    for (i = 1; i <= *m; ++i) {
				B(i,j) -= A(k,j) * B(i,k);
/* L180: */
			    }
			}
/* L190: */
		    }
		    if (nounit) {
			temp = 1. / A(j,j);
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = temp * B(i,j);
/* L200: */
			}
		    }
/* L210: */
		}
	    } else {
		for (j = *n; j >= 1; --j) {
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = *alpha * B(i,j);
/* L220: */
			}
		    }
		    for (k = j + 1; k <= *n; ++k) {
			if (A(k,j) != 0.) {
			    for (i = 1; i <= *m; ++i) {
				B(i,j) -= A(k,j) * B(i,k);
/* L230: */
			    }
			}
/* L240: */
		    }
		    if (nounit) {
			temp = 1. / A(j,j);
			for (i = 1; i <= *m; ++i) {
			    B(i,j) = temp * B(i,j);
/* L250: */
			}
		    }
/* L260: */
		}
	    }
	} else {

/*           Form  B := alpha*B*inv( A' ). */

	    if (upper) {
		for (k = *n; k >= 1; --k) {
		    if (nounit) {
			temp = 1. / A(k,k);
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = temp * B(i,k);
/* L270: */
			}
		    }
		    for (j = 1; j <= k-1; ++j) {
			if (A(j,k) != 0.) {
			    temp = A(j,k);
			    for (i = 1; i <= *m; ++i) {
				B(i,j) -= temp * B(i,k);
/* L280: */
			    }
			}
/* L290: */
		    }
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = *alpha * B(i,k);
/* L300: */
			}
		    }
/* L310: */
		}
	    } else {
		for (k = 1; k <= *n; ++k) {
		    if (nounit) {
			temp = 1. / A(k,k);
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = temp * B(i,k);
/* L320: */
			}
		    }
		    for (j = k + 1; j <= *n; ++j) {
			if (A(j,k) != 0.) {
			    temp = A(j,k);
			    for (i = 1; i <= *m; ++i) {
				B(i,j) -= temp * B(i,k);
/* L330: */
			    }
			}
/* L340: */
		    }
		    if (*alpha != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = *alpha * B(i,k);
/* L350: */
			}
		    }
/* L360: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRSM . */

} /* dtrsm_ */
#endif

#ifndef ML_DLASWP_FUNC

/* Subroutine */ int MLFORTRAN(dlaswp)(integer *n, doublereal *a, integer *lda, integer 
	*k1, integer *k2, integer *ipiv, integer *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASWP performs a series of row interchanges on the matrix A.   
    One row interchange is initiated for each of rows K1 through K2 of A. 
  

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of columns of the matrix A.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the matrix of column dimension N to which the row   
            interchanges will be applied.   
            On exit, the permuted matrix.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   

    K1      (input) INTEGER   
            The first element of IPIV for which a row interchange will   
            be done.   

    K2      (input) INTEGER   
            The last element of IPIV for which a row interchange will   
            be done.   

    IPIV    (input) INTEGER array, dimension (M*ML_abs(INCX))   
            The vector of pivot indices.  Only the elements in positions 
  
            K1 through K2 of IPIV are accessed.   
            IPIV(K) = L implies rows K and L are to be interchanged.   

    INCX    (input) INTEGER   
            The increment between successive values of IPIV.  If IPIV   
            is negative, the pivots are applied in reverse order.   

   ===================================================================== 
  


       Interchange row I with row IPIV(I) for each of rows K1 through K2. 
  

    
   Parameter adjustments   
       Function Body */
    /* Local variables */
    static integer i;
    static integer ip, ix;


#define IPIV(I) ipiv[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (*incx == 0) {
	return 0;
    }
    if (*incx > 0) {
	ix = *k1;
    } else {
	ix = (1 - *k2) * *incx + 1;
    }
    if (*incx == 1) {
	for (i = *k1; i <= *k2; ++i) {
	    ip = IPIV(i);
	    if (ip != i) {
		MLFORTRAN(dswap)(n, &A(i,1), lda, &A(ip,1), lda);
	    }
/* L10: */
	}
    } else if (*incx > 1) {
	for (i = *k1; i <= *k2; ++i) {
	    ip = IPIV(ix);
	    if (ip != i) {
		MLFORTRAN(dswap)(n, &A(i,1), lda, &A(ip,1), lda);
	    }
	    ix += *incx;
/* L20: */
	}
    } else if (*incx < 0) {
	for (i = *k2; i >= *k1; --i) {
	    ip = IPIV(ix);
	    if (ip != i) {
		MLFORTRAN(dswap)(n, &A(i,1), lda, &A(ip,1), lda);
	    }
	    ix += *incx;
/* L30: */
	}
    }

    return 0;

/*     End of DLASWP */

} /* dlaswp_ */
#endif

#ifndef ML_DGETF2_FUNC

/* Subroutine */ int MLFORTRAN(dgetf2)(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1992   


    Purpose   
    =======   

    DGETF2 computes an LU factorization of a general m-by-n matrix A   
    using partial pivoting with row interchanges.   

    The factorization has the form   
       A = P * L * U   
    where P is a permutation matrix, L is lower triangular with unit   
    diagonal elements (lower trapezoidal if m > n), and U is upper   
    triangular (upper trapezoidal if m < n).   

    This is the right-looking Level 2 BLAS version of the algorithm.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix to be factored.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,M).   

    IPIV    (output) INTEGER array, dimension (ML_min(M,N))   
            The pivot indices; for 1 <= i <= ML_min(M,N), row i of the   
            matrix was interchanged with row IPIV(i).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, U(k,k) is exactly zero. The factorization   
                 has been completed, but the factor U is exactly   
                 singular, and division by zero will occur if it is used 
  
                 to solve a system of equations.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b6 = -1.;
    
    /* System generated locals */
    integer   i__1, i__2, i__3;
    doublereal d__1;
    /* Local variables */
    static integer j;
    static integer jp;



#define IPIV(I) ipiv[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ML_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DGETF2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return 0;
    }

    i__1 = ML_min(*m,*n);
    for (j = 1; j <= ML_min(*m,*n); ++j) {

/*        Find pivot and test for singularity. */

	i__2 = *m - j + 1;
	jp = j - 1 + MLFORTRAN(idamax)(&i__2, &A(j,j), &c__1);
	IPIV(j) = jp;
	if (A(jp,j) != 0.) {

/*           Apply the interchange to columns 1:N. */

	    if (jp != j) {
		MLFORTRAN(dswap)(n, &A(j,1), lda, &A(jp,1), lda);
	    }

/*           Compute elements J+1:M of J-th column. */

	    if (j < *m) {
		i__2 = *m - j;
		d__1 = 1. / A(j,j);
		MLFORTRAN(dscal)(&i__2, &d__1, &A(j+1,j), &c__1);
	    }

	} else if (*info == 0) {

	    *info = j;
	}

	if (j < ML_min(*m,*n)) {

/*           Update trailing submatrix. */

	    i__2 = *m - j;
	    i__3 = *n - j;
	    MLFORTRAN(dger)(&i__2, &i__3, &c_b6, &A(j+1,j), &c__1, &A(j,j+1), lda, &A(j+1,j+1), lda);
	}
/* L10: */
    }
    return 0;

/*     End of DGETF2 */

} /* dgetf2_ */
#endif

#ifndef ML_DSWAP_FUNC

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Subroutine */ int MLFORTRAN(dswap)(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{


    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     interchanges two vectors.   
       uses unrolled loops for increments equal one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal   
           to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    for (i = 1; i <= *n; ++i) {
	dtemp = DX(ix);
	DX(ix) = DY(iy);
	DY(iy) = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1   


         clean-up loop */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    for (i = 1; i <= m; ++i) {
	dtemp = DX(i);
	DX(i) = DY(i);
	DY(i) = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 3) {
	dtemp = DX(i);
	DX(i) = DY(i);
	DY(i) = dtemp;
	dtemp = DX(i + 1);
	DX(i + 1) = DY(i + 1);
	DY(i + 1) = dtemp;
	dtemp = DX(i + 2);
	DX(i + 2) = DY(i + 2);
	DY(i + 2) = dtemp;
/* L50: */
    }
    return 0;
} /* dswap_ */
#endif
#ifndef ML_DASUM_FUNC
doublereal MLFORTRAN(dasum)(integer *n, doublereal *dx, integer *incx)
{


    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer nincx, mp1;


/*     takes the sum of the absolute values.
       jack dongarra, linpack, 3/11/78.
       modified 3/93 to return if incx .le. 0.
       modified 12/3/93, array(1) declarations changed to array(*)



   Parameter adjustments
       Function Body */
#define DX(I) dx[(I)-1]


    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0 || *incx <= 0) {
        return ret_val;
    }
    if (*incx == 1) {
        goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
        dtemp += (d__1 = DX(i), ML_abs(d__1));
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for increment equal to 1


          clean-up loop */

L20:
    m = *n % 6;
    if (m == 0) {
        goto L40;
    }
    for (i = 1; i <= m; ++i) {
        dtemp += (d__1 = DX(i), ML_abs(d__1));
/* L30: */
    }
    if (*n < 6) {
        goto L60;
    }
L40:
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 6) {
        dtemp = dtemp + (d__1 = DX(i), ML_abs(d__1)) + (d__2 = DX(i + 1), ML_abs(
                d__2)) + (d__3 = DX(i + 2), ML_abs(d__3)) + (d__4 = DX(i + 3),
                ML_abs(d__4)) + (d__5 = DX(i + 4), ML_abs(d__5)) + (d__6 = DX(i + 5)
                , ML_abs(d__6));
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* dasum_ */
#endif

#ifndef ML_DAXPY_FUNC
/* Subroutine */ int MLFORTRAN(daxpy)(integer *n, doublereal *da, doublereal *dx,
        integer *incx, doublereal *dy, integer *incy)
{



    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a vector plus a vector.
       uses unrolled loops for increments equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*)



   Parameter adjustments
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
        return 0;
    }
    if (*da == 0.) {
        return 0;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

/*        code for unequal increments or equal increments
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    for (i = 1; i <= *n; ++i) {
        DY(iy) += *da * DX(ix);
        ix += *incx;
        iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1


          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
        goto L40;
    }
    for (i = 1; i <= m; ++i) {
        DY(i) += *da * DX(i);
/* L30: */
    }
    if (*n < 4) {
        return 0;
    }
L40:
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 4) {
        DY(i) += *da * DX(i);
        DY(i + 1) += *da * DX(i + 1);
        DY(i + 2) += *da * DX(i + 2);
        DY(i + 3) += *da * DX(i + 3);
/* L50: */
    }
    return 0;
} /* daxpy_ */

#endif

#ifndef ML_DDOT_FUNC
doublereal MLFORTRAN(ddot)(integer *n, doublereal *dx, integer *incx, doublereal *dy,
        integer *incy)
{


    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer i, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     forms the dot product of two vectors.
       uses unrolled loops for increments equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*)



   Parameter adjustments
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
        return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
        goto L20;
    }

/*        code for unequal increments or equal increments
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
        iy = (-(*n) + 1) * *incy + 1;
    }
    for (i = 1; i <= *n; ++i) {
        dtemp += DX(ix) * DY(iy);
        ix += *incx;
        iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
        goto L40;
    }
    for (i = 1; i <= m; ++i) {
        dtemp += DX(i) * DY(i);
/* L30: */
    }
    if (*n < 5) {
        goto L60;
    }
L40:
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 5) {
        dtemp = dtemp + DX(i) * DY(i) + DX(i + 1) * DY(i + 1) + DX(i + 2) *
                DY(i + 2) + DX(i + 3) * DY(i + 3) + DX(i + 4) * DY(i + 4);
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

#endif


/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifndef ML_DSCAL_FUNC

/* Subroutine */ int MLFORTRAN(dscal)(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{



    /* Local variables */
    static integer i, m, nincx, mp1;


/*     scales a vector by a constant.   
       uses unrolled loops for increment equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DX(I) dx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	DX(i) = *da * DX(i);
/* L10: */
    }
    return 0;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    for (i = 1; i <= m; ++i) {
	DX(i) = *da * DX(i);
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 5) {
	DX(i) = *da * DX(i);
	DX(i + 1) = *da * DX(i + 1);
	DX(i + 2) = *da * DX(i + 2);
	DX(i + 3) = *da * DX(i + 3);
	DX(i + 4) = *da * DX(i + 4);
/* L50: */
    }
    return 0;
} /* dscal_ */
#endif


/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


#ifndef ML_DGEMM_FUNC
/* Subroutine */ int MLFORTRAN(dgemm)(char *transa, char *transb, integer *m, integer *
	n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *beta, doublereal *c, integer 
	*ldc,int dummy1,int dummy2)
{


    /* Local variables */
    static integer info;
    static logical nota, notb;
    static doublereal temp;
    static integer i, j, l;
    static integer nrowa, nrowb;


/*  Purpose   
    =======   

    DGEMM  performs one of the matrix-matrix operations   

       C := alpha*op( A )*op( B ) + beta*C,   

    where  op( X ) is one of   

       op( X ) = X   or   op( X ) = X',   

    alpha and beta are scalars, and A, B and C are matrices, with op( A ) 
  
    an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. 
  

    Parameters   
    ==========   

    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in 
  
             the matrix multiplication as follows:   

                TRANSA = 'N' or 'n',  op( A ) = A.   

                TRANSA = 'T' or 't',  op( A ) = A'.   

                TRANSA = 'C' or 'c',  op( A ) = A'.   

             Unchanged on exit.   

    TRANSB - CHARACTER*1.   
             On entry, TRANSB specifies the form of op( B ) to be used in 
  
             the matrix multiplication as follows:   

                TRANSB = 'N' or 'n',  op( B ) = B.   

                TRANSB = 'T' or 't',  op( B ) = B'.   

                TRANSB = 'C' or 'c',  op( B ) = B'.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry,  M  specifies  the number  of rows  of the  matrix 
  
             op( A )  and of the  matrix  C.  M  must  be at least  zero. 
  
             Unchanged on exit.   

    N      - INTEGER.   
             On entry,  N  specifies the number  of columns of the matrix 
  
             op( B ) and the number of columns of the matrix C. N must be 
  
             at least zero.   
             Unchanged on exit.   

    K      - INTEGER.   
             On entry,  K  specifies  the number of columns of the matrix 
  
             op( A ) and the number of rows of the matrix op( B ). K must 
  
             be at least  zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is 
  
             k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.   
             Before entry with  TRANSA = 'N' or 'n',  the leading  m by k 
  
             part of the array  A  must contain the matrix  A,  otherwise 
  
             the leading  k by m  part of the array  A  must contain  the 
  
             matrix A.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. When  TRANSA = 'N' or 'n' then 
  
             LDA must be at least  ML_max( 1, m ), otherwise  LDA must be at 
  
             least  ML_max( 1, k ).   
             Unchanged on exit.   

    B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is 
  
             n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.   
             Before entry with  TRANSB = 'N' or 'n',  the leading  k by n 
  
             part of the array  B  must contain the matrix  B,  otherwise 
  
             the leading  n by k  part of the array  B  must contain  the 
  
             matrix B.   
             Unchanged on exit.   

    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared 
  
             in the calling (sub) program. When  TRANSB = 'N' or 'n' then 
  
             LDB must be at least  ML_max( 1, k ), otherwise  LDB must be at 
  
             least  ML_max( 1, n ).   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry,  BETA  specifies the scalar  beta.  When  BETA  is 
  
             supplied as zero then C need not be set on input.   
             Unchanged on exit.   

    C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry, the leading  m by n  part of the array  C must 
  
             contain the matrix  C,  except when  beta  is zero, in which 
  
             case C need not be set on entry.   
             On exit, the array  C  is overwritten by the  m by n  matrix 
  
             ( alpha*op( A )*op( B ) + beta*C ).   

    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared 
  
             in  the  calling  (sub)  program.   LDC  must  be  at  least 
  
             ML_max( 1, m ).   
             Unchanged on exit.   


    Level 3 Blas routine.   

    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   



       Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not 
  
       transposed and set  NROWA, NCOLA and  NROWB  as the number of rows 
  
       and  columns of  A  and the  number of  rows  of  B  respectively. 
  

    
   Parameter adjustments   
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    nota = MLFORTRAN(lsame)(transa, "N");
    notb = MLFORTRAN(lsame)(transb, "N");
    if (nota) {
	nrowa = *m;
    } else {
	nrowa = *k;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }

/*     Test the input parameters. */

    info = 0;
    if (! nota && ! MLFORTRAN(lsame)(transa, "C") && ! MLFORTRAN(lsame)(transa, "T")) {
	info = 1;
    } else if (! notb && ! MLFORTRAN(lsame)(transb, "C") && ! MLFORTRAN(lsame)(transb, 
	    "T")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < ML_max(1,nrowa)) {
	info = 8;
    } else if (*ldb < ML_max(1,nrowb)) {
	info = 10;
    } else if (*ldc < ML_max(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	MLFORTRAN(xerbla)("DGEMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || ((*alpha == 0. || *k == 0) && *beta == 1.)) {
	return 0;
    }

/*     And if  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (*beta == 0.) {
	    for (j = 1; j <= *n; ++j) {
		for (i = 1; i <= *m; ++i) {
		    C(i,j) = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    for (j = 1; j <= *n; ++j) {
		for (i = 1; i <= *m; ++i) {
		    C(i,j) = *beta * C(i,j);
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (notb) {
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    for (i = 1; i <= *m; ++i) {
			C(i,j) = 0.;
/* L50: */
		    }
		} else if (*beta != 1.) {
		    for (i = 1; i <= *m; ++i) {
			C(i,j) = *beta * C(i,j);
/* L60: */
		    }
		}
		for (l = 1; l <= *k; ++l) {
		    if (B(l,j) != 0.) {
			temp = *alpha * B(l,j);
			for (i = 1; i <= *m; ++i) {
			    C(i,j) += temp * A(i,l);
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else {

/*           Form  C := alpha*A'*B + beta*C */

	    for (j = 1; j <= *n; ++j) {
		for (i = 1; i <= *m; ++i) {
		    temp = 0.;
		    for (l = 1; l <= *k; ++l) {
			temp += A(l,i) * B(l,j);
/* L100: */
		    }
		    if (*beta == 0.) {
			C(i,j) = *alpha * temp;
		    } else {
			C(i,j) = *alpha * temp + *beta * C(i,j);
		    }
/* L110: */
		}
/* L120: */
	    }
	}
    } else {
	if (nota) {

/*           Form  C := alpha*A*B' + beta*C */

	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    for (i = 1; i <= *m; ++i) {
			C(i,j) = 0.;
/* L130: */
		    }
		} else if (*beta != 1.) {
		    for (i = 1; i <= *m; ++i) {
			C(i,j) = *beta * C(i,j);
/* L140: */
		    }
		}
		for (l = 1; l <= *k; ++l) {
		    if (B(j,l) != 0.) {
			temp = *alpha * B(j,l);
			for (i = 1; i <= *m; ++i) {
			    C(i,j) += temp * A(i,l);
/* L150: */
			}
		    }
/* L160: */
		}
/* L170: */
	    }
	} else {

/*           Form  C := alpha*A'*B' + beta*C */

	    for (j = 1; j <= *n; ++j) {
		for (i = 1; i <= *m; ++i) {
		    temp = 0.;
		    for (l = 1; l <= *k; ++l) {
			temp += A(l,i) * B(j,l);
/* L180: */
		    }
		    if (*beta == 0.) {
			C(i,j) = *alpha * temp;
		    } else {
			C(i,j) = *alpha * temp + *beta * C(i,j);
		    }
/* L190: */
		}
/* L200: */
	    }
	}
    }

    return 0;

/*     End of DGEMM . */

} /* dgemm_ */
#endif


/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifndef ML_IDAMAX_FUNC

integer MLFORTRAN(idamax)(integer *n, doublereal *dx, integer *incx)
{


    /* System generated locals */
    integer ret_val;
    doublereal d__1;

    /* Local variables */
    static doublereal dmax__;
    static integer i, ix;


/*     finds the index of element having max. absolute value.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DX(I) dx[(I)-1]


    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    dmax__ = ML_abs(DX(1));
    ix += *incx;
    for (i = 2; i <= *n; ++i) {
	if ((d__1 = DX(ix), ML_abs(d__1)) <= dmax__) {
	    goto L5;
	}
	ret_val = i;
	dmax__ = (d__1 = DX(ix), ML_abs(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    dmax__ = ML_abs(DX(1));
    for (i = 2; i <= *n; ++i) {
	if ((d__1 = DX(i), ML_abs(d__1)) <= dmax__) {
	    goto L30;
	}
	ret_val = i;
	dmax__ = (d__1 = DX(i), ML_abs(d__1));
L30:
	;
    }
    return ret_val;
} /* idamax_ */
#endif


/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


#ifndef ML_DGER_FUNC
/* Subroutine */ int MLFORTRAN(dger)(integer *m, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *a, integer *lda)
{


    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i, j, ix, jy, kx;


/*  Purpose   
    =======   

    DGER   performs the rank 1 operation   

       A := alpha*x*y' + A,   

    where alpha is a scalar, x is an m element vector, y is an n element 
  
    vector and A is an m by n matrix.   

    Parameters   
    ==========   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( m - 1 )*ML_abs( INCX ) ).   
             Before entry, the incremented array X must contain the m   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*ML_abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients. On exit, A is   
             overwritten by the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             ML_max( 1, m ).   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < ML_max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	MLFORTRAN(xerbla)("DGER  ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0.) {
	return 0;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	for (j = 1; j <= *n; ++j) {
	    if (Y(jy) != 0.) {
		temp = *alpha * Y(jy);
		for (i = 1; i <= *m; ++i) {
		    A(i,j) += X(i) * temp;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	for (j = 1; j <= *n; ++j) {
	    if (Y(jy) != 0.) {
		temp = *alpha * Y(jy);
		ix = kx;
		for (i = 1; i <= *m; ++i) {
		    A(i,j) += X(ix) * temp;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }

    return 0;

/*     End of DGER  . */

} /* dger_ */
#endif


#ifndef ML_ILAENV_FUNC
integer MLFORTRAN(ilaenv)(integer *ispec, char *name, char *opts, integer *n1, integer *
	n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ILAENV is called from the LAPACK routines to choose problem-dependent 
  
    parameters for the local environment.  See ISPEC for a description of 
  
    the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

    This routine will not function correctly if it is converted to all   
    lower case.  Converting it to all upper case is allowed.   

    Arguments   
    =========   

    ISPEC   (input) INTEGER   
            Specifies the parameter to be returned as the value of   
            ILAENV.   
            = 1: the optimal blocksize; if this value is 1, an unblocked 
  
                 algorithm will give the best performance.   
            = 2: the minimum block size for which the block routine   
                 should be used; if the usable block size is less than   
                 this value, an unblocked routine should be used.   
            = 3: the crossover point (in a block routine, for N less   
                 than this value, an unblocked routine should be used)   
            = 4: the number of shifts, used in the nonsymmetric   
                 eigenvalue routines   
            = 5: the minimum column dimension for blocking to be used;   
                 rectangular blocks must have dimension at least k by m, 
  
                 where k is given by ILAENV(2,...) and m by ILAENV(5,...) 
  
            = 6: the crossover point for the SVD (when reducing an m by n 
  
                 matrix to bidiagonal form, if ML_max(m,n)/ML_min(m,n) exceeds 
  
                 this value, a QR factorization is used first to reduce   
                 the matrix to a triangular form.)   
            = 7: the number of processors   
            = 8: the crossover point for the multishift QR and QZ methods 
  
                 for nonsymmetric eigenvalue problems.   

    NAME    (input) CHARACTER*(*)   
            The name of the calling subroutine, in either upper case or   
            lower case.   

    OPTS    (input) CHARACTER*(*)   
            The character options to the subroutine NAME, concatenated   
            into a single character string.  For example, UPLO = 'U',   
            TRANS = 'T', and DIAG = 'N' for a triangular routine would   
            be specified as OPTS = 'UTN'.   

    N1      (input) INTEGER   
    N2      (input) INTEGER   
    N3      (input) INTEGER   
    N4      (input) INTEGER   
            Problem dimensions for the subroutine NAME; these may not all 
  
            be required.   

   (ILAENV) (output) INTEGER   
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if ILAENV = -k, the k-th argument had an illegal value. 
  

    Further Details   
    ===============   

    The following conventions have been used when calling ILAENV from the 
  
    LAPACK routines:   
    1)  OPTS is a concatenation of all of the character options to   
        subroutine NAME, in the same order that they appear in the   
        argument list for NAME, even if they are not used in determining 
  
        the value of the parameter specified by ISPEC.   
    2)  The problem dimensions N1, N2, N3, N4 are specified in the order 
  
        that they appear in the argument list for NAME.  N1 is used   
        first, N2 second, and so on, and unused problem dimensions are   
        passed a value of -1.   
    3)  The parameter value returned by ILAENV is checked for validity in 
  
        the calling subroutine.  For example, ILAENV is used to retrieve 
  
        the optimal blocksize for STRTRI as follows:   

        NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )   
        IF( NB.LE.1 ) NB = ML_MAX( 1, N )   

    ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    integer ret_val;
    /* Local variables */
    static integer i;
    static logical cname, sname;
    static integer nbmin;
    static char c1[1], c2[2], c3[3], c4[2];
    static integer ic, nb, iz, nx;
    static char subnam[6];



    switch (*ispec) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L100;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L600;
	case 7:  goto L700;
	case 8:  goto L800;
    }

/*     Invalid value for ISPEC */

    ret_val = -1;
    return ret_val;

L100:

/*     Convert NAME to upper case if the first character is lower case. */

    ret_val = 1;
    ml_s_copy(subnam, name, 6L, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

	if (ic >= 97 && ic <= 122) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i = 2; i <= 6; ++i) {
		ic = *(unsigned char *)&subnam[i - 1];
		if (ic >= 97 && ic <= 122) {
		    *(unsigned char *)&subnam[i - 1] = (char) (ic - 32);
		}
/* L10: */
	    }
	}

    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

	if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && 
		ic <= 169)) {
	    *(unsigned char *)subnam = (char) (ic + 64);
	    for (i = 2; i <= 6; ++i) {
		ic = *(unsigned char *)&subnam[i - 1];
		if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 
			162 && ic <= 169)) {
		    *(unsigned char *)&subnam[i - 1] = (char) (ic + 64);
		}
/* L20: */
	    }
	}

    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

	if (ic >= 225 && ic <= 250) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i = 2; i <= 6; ++i) {
		ic = *(unsigned char *)&subnam[i - 1];
		if (ic >= 225 && ic <= 250) {
		    *(unsigned char *)&subnam[i - 1] = (char) (ic - 32);
		}
/* L30: */
	    }
	}
    }

    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if (! (cname || sname)) {
	return ret_val;
    }
    ml_s_copy(c2, subnam + 1, 2L, 2L);
    ml_s_copy(c3, subnam + 3, 3L, 3L);
    ml_s_copy(c4, c3 + 1, 2L, 2L);

    switch (*ispec) {
	case 1:  goto L110;
	case 2:  goto L200;
	case 3:  goto L300;
    }

L110:

/*     ISPEC = 1:  block size   

       In these examples, separate code is provided for setting NB for   
       real and complex.  We assume that NB will take the same value in   
       single or double precision. */

    nb = 1;

    if (ml_s_cmp(c2, "GE", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (ml_s_cmp(c3, "QRF", 3L, 3L) == 0 || ml_s_cmp(c3, "RQF", 3L, 3L) 
		== 0 || ml_s_cmp(c3, "LQF", 3L, 3L) == 0 || ml_s_cmp(c3, "QLF", 3L, 
		3L) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (ml_s_cmp(c3, "HRD", 3L, 3L) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (ml_s_cmp(c3, "BRD", 3L, 3L) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (ml_s_cmp(c3, "TRI", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (ml_s_cmp(c2, "PO", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (ml_s_cmp(c2, "SY", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (sname && ml_s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nb = 1;
	} else if (sname && ml_s_cmp(c3, "GST", 3L, 3L) == 0) {
	    nb = 64;
	}
    } else if (cname && ml_s_cmp(c2, "HE", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    nb = 64;
	} else if (ml_s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nb = 1;
	} else if (ml_s_cmp(c3, "GST", 3L, 3L) == 0) {
	    nb = 64;
	}
    } else if (sname && ml_s_cmp(c2, "OR", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	}
    } else if (cname && ml_s_cmp(c2, "UN", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nb = 32;
	    }
	}
    } else if (ml_s_cmp(c2, "GB", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (ml_s_cmp(c2, "PB", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (ml_s_cmp(c2, "TR", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRI", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (ml_s_cmp(c2, "LA", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "UUM", 3L, 3L) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (sname && ml_s_cmp(c2, "ST", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "EBZ", 3L, 3L) == 0) {
	    nb = 1;
	}
    }
    ret_val = nb;
    return ret_val;

L200:

/*     ISPEC = 2:  minimum block size */

    nbmin = 2;
    if (ml_s_cmp(c2, "GE", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "QRF", 3L, 3L) == 0 || ml_s_cmp(c3, "RQF", 3L, 3L) == 0 || 
		ml_s_cmp(c3, "LQF", 3L, 3L) == 0 || ml_s_cmp(c3, "QLF", 3L, 3L) == 
		0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (ml_s_cmp(c3, "HRD", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (ml_s_cmp(c3, "BRD", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (ml_s_cmp(c3, "TRI", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	}
    } else if (ml_s_cmp(c2, "SY", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRF", 3L, 3L) == 0) {
	    if (sname) {
		nbmin = 8;
	    } else {
		nbmin = 8;
	    }
	} else if (sname && ml_s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nbmin = 2;
	}
    } else if (cname && ml_s_cmp(c2, "HE", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nbmin = 2;
	}
    } else if (sname && ml_s_cmp(c2, "OR", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	}
    } else if (cname && ml_s_cmp(c2, "UN", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nbmin = 2;
	    }
	}
    }
    ret_val = nbmin;
    return ret_val;

L300:

/*     ISPEC = 3:  crossover point */

    nx = 0;
    if (ml_s_cmp(c2, "GE", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "QRF", 3L, 3L) == 0 || ml_s_cmp(c3, "RQF", 3L, 3L) == 0 || 
		ml_s_cmp(c3, "LQF", 3L, 3L) == 0 || ml_s_cmp(c3, "QLF", 3L, 3L) == 
		0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (ml_s_cmp(c3, "HRD", 3L, 3L) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (ml_s_cmp(c3, "BRD", 3L, 3L) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	}
    } else if (ml_s_cmp(c2, "SY", 2L, 2L) == 0) {
	if (sname && ml_s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nx = 1;
	}
    } else if (cname && ml_s_cmp(c2, "HE", 2L, 2L) == 0) {
	if (ml_s_cmp(c3, "TRD", 3L, 3L) == 0) {
	    nx = 1;
	}
    } else if (sname && ml_s_cmp(c2, "OR", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nx = 128;
	    }
	}
    } else if (cname && ml_s_cmp(c2, "UN", 2L, 2L) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (ml_s_cmp(c4, "QR", 2L, 2L) == 0 || ml_s_cmp(c4, "RQ", 2L, 2L) == 0 
		    || ml_s_cmp(c4, "LQ", 2L, 2L) == 0 || ml_s_cmp(c4, "QL", 2L, 2L)
		     == 0 || ml_s_cmp(c4, "HR", 2L, 2L) == 0 || ml_s_cmp(c4, "TR", 
		    2L, 2L) == 0 || ml_s_cmp(c4, "BR", 2L, 2L) == 0) {
		nx = 128;
	    }
	}
    }
    ret_val = nx;
    return ret_val;

L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

    ret_val = 6;
    return ret_val;

L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

    ret_val = 2;
    return ret_val;

L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

    ret_val = (integer) ((real) ML_min(*n1,*n2) * 1.6f);
    return ret_val;

L700:

/*     ISPEC = 7:  number of processors (not used) */

    ret_val = 1;
    return ret_val;

L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    ret_val = 50;
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */
#endif

#ifndef ML_DGEQRF_FUNC

/* Subroutine */ int MLFORTRAN(dgeqrf)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGEQRF computes a QR factorization of a real M-by-N matrix A:   
    A = Q * R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the elements on and above the diagonal of the array 
  
            contain the ML_min(M,N)-by-N upper trapezoidal matrix R (R is   
            upper triangular if m >= n); the elements below the diagonal, 
  
            with the array TAU, represent the orthogonal matrix Q as a   
            product of ML_min(m,n) elementary reflectors (see Further   
            Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,M).   

    TAU     (output) DOUBLE PRECISION array, dimension (ML_min(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= ML_max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is   
            the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = ML_min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), 
  
    and tau in TAU(i).   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer   i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i, k, nbmin, iinfo;
    static integer ib, nb;
    static integer nx;
    static integer ldwork, iws;



#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ML_max(1,*m)) {
	*info = -4;
    } else if (*lwork < ML_max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DGEQRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    k = ML_min(*m,*n);
    if (k == 0) {
	WORK(1) = 1.;
	return 0;
    }

/*     Determine the block size. */

    nb = MLFORTRAN(ilaenv)(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = MLFORTRAN(ilaenv)(&c__3, "DGEQRF", " ", m, n, &c_n1, &c_n1, 6L,
		 1L);
	nx = ML_max(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduc
e NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = MLFORTRAN(ilaenv)(&c__2, "DGEQRF", " ", m, n, &c_n1, &
			c_n1, 6L, 1L);
		nbmin = ML_max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i = 1; nb < 0 ? i >= k-nx : i <= k-nx; i += nb) {
/* Computing MIN */
	    i__3 = k - i + 1;
	    ib = ML_min(i__3,nb);

/*           Compute the QR factorization of the current block   
             A(i:m,i:i+ib-1) */

	    i__3 = *m - i + 1;
	    MLFORTRAN(dgeqr2)(&i__3, &ib, &A(i,i), lda, &TAU(i), &WORK(1), &
		    iinfo);
	    if (i + ib <= *n) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__3 = *m - i + 1;
		MLFORTRAN(dlarft)("Forward", "Columnwise", &i__3, &ib, &A(i,i), lda, &TAU(i), &WORK(1), &ldwork);

/*              Apply H' to A(i:m,i+ib:n) from the left */

		i__3 = *m - i + 1;
		i__4 = *n - i - ib + 1;
		MLFORTRAN(dlarfb)("Left", "Transpose", "Forward", "Columnwise", &i__3, &
			i__4, &ib, &A(i,i), lda, &WORK(1), &ldwork,
			 &A(i,i+ib), lda, &WORK(ib + 1), &
			ldwork);
	    }
/* L10: */
	}
    } else {
	i = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i <= k) {
	i__2 = *m - i + 1;
	i__1 = *n - i + 1;
	MLFORTRAN(dgeqr2)(&i__2, &i__1, &A(i,i), lda, &TAU(i), &WORK(1), &
		iinfo);
    }

    WORK(1) = (doublereal) iws;
    return 0;

/*     End of DGEQRF */

} /* dgeqrf_ */
#endif


#ifndef ML_DGEQR2_FUNC
/* Subroutine */ int MLFORTRAN(dgeqr2)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DGEQR2 computes a QR factorization of a real m by n matrix A:   
    A = Q * R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix A.   
            On exit, the elements on and above the diagonal of the array 
  
            contain the ML_min(m,n) by n upper trapezoidal matrix R (R is   
            upper triangular if m >= n); the elements below the diagonal, 
  
            with the array TAU, represent the orthogonal matrix Q as a   
            product of elementary reflectors (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,M).   

    TAU     (output) DOUBLE PRECISION array, dimension (ML_min(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = ML_min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), 
  
    and tau in TAU(i).   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer   i__1, i__2, i__3;
    /* Local variables */
    static integer i, k;
    static doublereal aii;



#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ML_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DGEQR2", &i__1);
	return 0;
    }

    k = ML_min(*m,*n);

    i__1 = k;
    for (i = 1; i <= k; ++i) {

/*        Generate elementary reflector H(i) to annihilate A(i+1:m,i) 
*/

	i__2 = *m - i + 1;
/* Computing MIN */
	i__3 = i + 1;
	MLFORTRAN(dlarfg)(&i__2, &A(i,i), &A(ML_min(i+1,*m),i), &
		c__1, &TAU(i));
	if (i < *n) {

/*           Apply H(i) to A(i:m,i+1:n) from the left */

	    aii = A(i,i);
	    A(i,i) = 1.;
	    i__2 = *m - i + 1;
	    i__3 = *n - i;
	    MLFORTRAN(dlarf)("Left", &i__2, &i__3, &A(i,i), &c__1, &TAU(i), &
		    A(i,i+1), lda, &WORK(1));
	    A(i,i) = aii;
	}
/* L10: */
    }
    return 0;

/*     End of DGEQR2 */

} /* dgeqr2_ */
#endif


#ifndef ML_DLARFT_FUNC
/* Subroutine */ int MLFORTRAN(dlarft)(char *direct, char *storev, integer *n, integer *
	k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, 
	integer *ldt)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARFT forms the triangular factor T of a real block reflector H   
    of order n, which is defined as a product of k elementary reflectors. 
  

    If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular; 
  

    If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. 
  

    If STOREV = 'C', the vector which defines the elementary reflector   
    H(i) is stored in the i-th column of the array V, and   

       H  =  I - V * T * V'   

    If STOREV = 'R', the vector which defines the elementary reflector   
    H(i) is stored in the i-th row of the array V, and   

       H  =  I - V' * T * V   

    Arguments   
    =========   

    DIRECT  (input) CHARACTER*1   
            Specifies the order in which the elementary reflectors are   
            multiplied to form the block reflector:   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Specifies how the vectors which define the elementary   
            reflectors are stored (see also Further Details):   
            = 'C': columnwise   
            = 'R': rowwise   

    N       (input) INTEGER   
            The order of the block reflector H. N >= 0.   

    K       (input) INTEGER   
            The order of the triangular factor T (= the number of   
            elementary reflectors). K >= 1.   

    V       (input/output) DOUBLE PRECISION array, dimension   
                                 (LDV,K) if STOREV = 'C'   
                                 (LDV,N) if STOREV = 'R'   
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C', LDV >= ML_max(1,N); if STOREV = 'R', LDV >= K. 
  

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i).   

    T       (output) DOUBLE PRECISION array, dimension (LDT,K)   
            The k by k triangular factor T of the block reflector.   
            If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is 
  
            lower triangular. The rest of the array is not used.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    Further Details   
    ===============   

    The shape of the matrix V and the storage of the vectors which define 
  
    the H(i) is best illustrated by the following example with n = 5 and 
  
    k = 3. The elements equal to 1 are not stored; the corresponding   
    array elements are modified but restored on exit. The rest of the   
    array is not used.   

    DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R': 
  

                 V = (  1       )                 V = (  1 v1 v1 v1 v1 ) 
  
                     ( v1  1    )                     (     1 v2 v2 v2 ) 
  
                     ( v1 v2  1 )                     (        1 v3 v3 ) 
  
                     ( v1 v2 v3 )   
                     ( v1 v2 v3 )   

    DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R': 
  

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1       ) 
  
                     ( v1 v2 v3 )                     ( v2 v2 v2  1    ) 
  
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 ) 
  
                     (     1 v3 )   
                     (        1 )   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b8 = 0.;
    
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    /* Local variables */
    static integer i, j;
    static doublereal vii;



#define TAU(I) tau[(I)-1]

#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]

    if (*n == 0) {
	return 0;
    }

    if (MLFORTRAN(lsame)(direct, "F")) {
	i__1 = *k;
	for (i = 1; i <= *k; ++i) {
	    if (TAU(i) == 0.) {

/*              H(i)  =  I */

		i__2 = i;
		for (j = 1; j <= i; ++j) {
		    T(j,i) = 0.;
/* L10: */
		}
	    } else {

/*              general case */

		vii = V(i,i);
		V(i,i) = 1.;
		if (MLFORTRAN(lsame)(storev, "C")) {

/*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' 
* V(i:n,i) */

		    i__2 = *n - i + 1;
		    i__3 = i - 1;
		    d__1 = -TAU(i);
		    MLFORTRAN(dgemv)("Transpose", &i__2, &i__3, &d__1, &V(i,1), 
			    ldv, &V(i,i), &c__1, &c_b8, &T(1,i), &c__1
#ifdef FORTRAN_STRLEN
				     , strlen("Transpose")
#endif
				     );
		} else {

/*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) *
 V(i,i:n)' */

		    i__2 = i - 1;
		    i__3 = *n - i + 1;
		    d__1 = -TAU(i);
		    MLFORTRAN(dgemv)("No transpose", &i__2, &i__3, &d__1, &V(1,i), ldv, &V(i,i), ldv, &c_b8, &T(1,i), &c__1
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose")
#endif
);
		}
		V(i,i) = vii;

/*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */

		i__2 = i - 1;
		MLFORTRAN(dtrmv)("Upper", "No transpose", "Non-unit", &i__2, &T(1,1), ldt, &T(1,i), &c__1);
		T(i,i) = TAU(i);
	    }
/* L20: */
	}
    } else {
	for (i = *k; i >= 1; --i) {
	    if (TAU(i) == 0.) {

/*              H(i)  =  I */

		i__1 = *k;
		for (j = i; j <= *k; ++j) {
		    T(j,i) = 0.;
/* L30: */
		}
	    } else {

/*              general case */

		if (i < *k) {
		    if (MLFORTRAN(lsame)(storev, "C")) {
			vii = V(*n-*k+i,i);
			V(*n-*k+i,i) = 1.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(1:n-k+i,i+1
:k)' * V(1:n-k+i,i) */

			i__1 = *n - *k + i;
			i__2 = *k - i;
			d__1 = -TAU(i);
			MLFORTRAN(dgemv)("Transpose", &i__1, &i__2, &d__1, &V(1,i+1), ldv, &V(1,i), &c__1, &
				c_b8, &T(i+1,i), &c__1
#ifdef FORTRAN_STRLEN
				     , strlen("Transpose")
#endif
);
			V(*n-*k+i,i) = vii;
		    } else {
			vii = V(i,*n-*k+i);
			V(i,*n-*k+i) = 1.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(i+1:k,1:n-k
+i) * V(i,1:n-k+i)' */

			i__1 = *k - i;
			i__2 = *n - *k + i;
			d__1 = -TAU(i);
			MLFORTRAN(dgemv)("No transpose", &i__1, &i__2, &d__1, &V(i+1,1), ldv, &V(i,1), ldv, &c_b8, &
				T(i+1,i), &c__1
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose")
#endif
);
			V(i,*n-*k+i) = vii;
		    }

/*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,
i) */

		    i__1 = *k - i;
		    MLFORTRAN(dtrmv)("Lower", "No transpose", "Non-unit", &i__1, &T(i+1,i+1), ldt, &T(i+1,i)
			    , &c__1);
		}
		T(i,i) = TAU(i);
	    }
/* L40: */
	}
    }
    return 0;

/*     End of DLARFT */

} /* dlarft_ */
#endif


#ifndef ML_DLARFB_FUNC
/* Subroutine */ int MLFORTRAN(dlarfb)(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublereal *v, integer *
	ldv, doublereal *t, integer *ldt, doublereal *c, integer *ldc, 
	doublereal *work, integer *ldwork)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARFB applies a real block reflector H or its transpose H' to a   
    real m by n matrix C, from either the left or the right.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply H or H' from the Left   
            = 'R': apply H or H' from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply H (No transpose)   
            = 'T': apply H' (Transpose)   

    DIRECT  (input) CHARACTER*1   
            Indicates how H is formed from a product of elementary   
            reflectors   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Indicates how the vectors which define the elementary   
            reflectors are stored:   
            = 'C': Columnwise   
            = 'R': Rowwise   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    K       (input) INTEGER   
            The order of the matrix T (= the number of elementary   
            reflectors whose product defines the block reflector).   

    V       (input) DOUBLE PRECISION array, dimension   
                                  (LDV,K) if STOREV = 'C'   
                                  (LDV,M) if STOREV = 'R' and SIDE = 'L' 
  
                                  (LDV,N) if STOREV = 'R' and SIDE = 'R' 
  
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C' and SIDE = 'L', LDV >= ML_max(1,M);   
            if STOREV = 'C' and SIDE = 'R', LDV >= ML_max(1,N);   
            if STOREV = 'R', LDV >= K.   

    T       (input) DOUBLE PRECISION array, dimension (LDT,K)   
            The triangular k by k matrix T in the representation of the   
            block reflector.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by H*C or H'*C or C*H or C*H'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDA >= ML_max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)   

    LDWORK  (input) INTEGER   
            The leading dimension of the array WORK.   
            If SIDE = 'L', LDWORK >= ML_max(1,N);   
            if SIDE = 'R', LDWORK >= ML_max(1,M).   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b14 = 1.;
    static doublereal c_b25 = -1.;
    
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i, j;
    static char transt[1];




#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define T(I,J) t[(I)-1 + ((J)-1)* ( *ldt)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]
#undef WORK
#define WORK(I,J) work[(I)-1 + ((J)-1)* ( *ldwork)]

    if (*m <= 0 || *n <= 0) {
	return 0;
    }

    if (MLFORTRAN(lsame)(trans, "N")) {
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }

    if (MLFORTRAN(lsame)(storev, "C")) {

	if (MLFORTRAN(lsame)(direct, "F")) {

/*           Let  V =  ( V1 )    (first K rows)   
                       ( V2 )   
             where  V1  is unit lower triangular. */

	    if (MLFORTRAN(lsame)(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in 
WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    MLFORTRAN(dcopy)(n, &C(j,1), ldc, &WORK(1,j), &
			    c__1);
/* L10: */
		}

/*              W := W * V1 */

		MLFORTRAN(dtrmm)("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &V(1,1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen("No transpose"),strlen("Unit")
#endif
				 );
		if (*m > *k) {

/*                 W := W + C2'*V2 */

		    i__1 = *m - *k;
		    MLFORTRAN(dgemm)("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    C(*k+1,1), ldc, &V(*k+1,1), ldv,
			     &c_b14, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				     , strlen("Transpose"),strlen("No transpose")
#endif
);
		}

/*              W := W * T'  or  W * T */

		MLFORTRAN(dtrmm)("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen(transt),strlen("Non-unit")
#endif
				 );

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2 * W' */

		    i__1 = *m - *k;
		    MLFORTRAN(dgemm)("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    V(*k+1,1), ldv, &WORK(1,1), 
			    ldwork, &c_b14, &C(*k+1,1), ldc
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("Transpose")
#endif
)
			    ;
		}

/*              W := W * V1' */

		MLFORTRAN(dtrmm)("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			V(1,1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen("Transpose"),strlen("Unit")
#endif
				 );

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    for (i = 1; i <= *n; ++i) {
			C(j,i) -= WORK(i,j);
/* L20: */
		    }
/* L30: */
		}

	    } else if (MLFORTRAN(lsame)(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WOR
K)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    MLFORTRAN(dcopy)(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
/* L40: */
		}

/*              W := W * V1 */

		MLFORTRAN(dtrmm)("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &V(1,1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen("No transpose"),strlen("Unit")
#endif
				 );
		if (*n > *k) {

/*                 W := W + C2 * V2 */

		    i__1 = *n - *k;
		    MLFORTRAN(dgemm)("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &C(1,*k+1), ldc, &V(*k+1,1), ldv, &c_b14, &WORK(1,1), 
			    ldwork
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("No transpose")
#endif
);
		}

/*              W := W * T  or  W * T' */

		MLFORTRAN(dtrmm)("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen(trans),strlen("Non-unit")
#endif
);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C2 := C2 - W * V2' */

		    i__1 = *n - *k;
		    MLFORTRAN(dgemm)("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    WORK(1,1), ldwork, &V(*k+1,1), 
			    ldv, &c_b14, &C(1,*k+1), ldc
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("Transpose")
#endif
);
		}

/*              W := W * V1' */

		MLFORTRAN(dtrmm)("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			V(1,1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen("Transpose"),strlen("Unit")
#endif
);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    for (i = 1; i <= *m; ++i) {
			C(i,j) -= WORK(i,j);
/* L50: */
		    }
/* L60: */
		}
	    }

	} else {

/*           Let  V =  ( V1 )   
                       ( V2 )    (last K rows)   
             where  V2  is unit upper triangular. */

	    if (MLFORTRAN(lsame)(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in 
WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    MLFORTRAN(dcopy)(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
/* L70: */
		}

/*              W := W * V2 */

		MLFORTRAN(dtrmm)("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &V(*m-*k+1,1), ldv, &WORK(1,1), 
			ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen("No transpose"),strlen("Unit")
#endif
);
		if (*m > *k) {

/*                 W := W + C1'*V1 */

		    i__1 = *m - *k;
		    MLFORTRAN(dgemm)("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    C(1,1), ldc, &V(1,1), ldv, &c_b14, &
			    WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				     , strlen("Transpose"),strlen("No transpose")
#endif
);
		}

/*              W := W * T'  or  W * T */

		MLFORTRAN(dtrmm)("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen(transt),strlen("Non-unit")
#endif
);

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1 * W' */

		    i__1 = *m - *k;
		    MLFORTRAN(dgemm)("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    V(1,1), ldv, &WORK(1,1), ldwork, &
			    c_b14, &C(1,1), ldc
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("Transpose")
#endif
);
		}

/*              W := W * V2' */

		MLFORTRAN(dtrmm)("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			V(*m-*k+1,1), ldv, &WORK(1,1), 
			ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen("Transpose"),strlen("Unit")
#endif
				 );

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    for (i = 1; i <= *n; ++i) {
			C(*m-*k+j,i) -= WORK(i,j)
				;
/* L80: */
		    }
/* L90: */
		}

	    } else if (MLFORTRAN(lsame)(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WOR
K)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    MLFORTRAN(dcopy)(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
/* L100: */
		}

/*              W := W * V2 */

		MLFORTRAN(dtrmm)("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &V(*n-*k+1,1), ldv, &WORK(1,1), 
			ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen("No transpose"),strlen("Unit")
#endif
				 );
		if (*n > *k) {

/*                 W := W + C1 * V1 */

		    i__1 = *n - *k;
		    MLFORTRAN(dgemm)("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &C(1,1), ldc, &V(1,1), ldv, &
			    c_b14, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("No transpose")
#endif
);
		}

/*              W := W * T  or  W * T' */

		MLFORTRAN(dtrmm)("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen(trans),strlen("Non-unit")
#endif
);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C1 := C1 - W * V1' */

		    i__1 = *n - *k;
		    MLFORTRAN(dgemm)("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    WORK(1,1), ldwork, &V(1,1), ldv, &
			    c_b14, &C(1,1), ldc
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("Transpose")
#endif
				     );
		}

/*              W := W * V2' */

		MLFORTRAN(dtrmm)("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			V(*n-*k+1,1), ldv, &WORK(1,1), 
			ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen("Transpose"),strlen("Unit")
#endif
				 );

/*              C2 := C2 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    for (i = 1; i <= *m; ++i) {
			C(i,*n-*k+j) -= WORK(i,j);
/* L110: */
		    }
/* L120: */
		}
	    }
	}

    } else if (MLFORTRAN(lsame)(storev, "R")) {

	if (MLFORTRAN(lsame)(direct, "F")) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns)   
             where  V1  is unit upper triangular. */

	    if (MLFORTRAN(lsame)(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored i
n WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    MLFORTRAN(dcopy)(n, &C(j,1), ldc, &WORK(1,j), &
			    c__1);
/* L130: */
		}

/*              W := W * V1' */

		MLFORTRAN(dtrmm)("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			V(1,1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen("Transpose"),strlen("Unit")
#endif
);
		if (*m > *k) {

/*                 W := W + C2'*V2' */

		    i__1 = *m - *k;
		    MLFORTRAN(dgemm)("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(*k+1,1), ldc, &V(1,*k+1), 
			    ldv, &c_b14, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				     , strlen("Transpose"),strlen("Transpose")
#endif
);
		}

/*              W := W * T'  or  W * T */

		MLFORTRAN(dtrmm)("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen(transt),strlen("Non-unit")
#endif
);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2' * W' */

		    i__1 = *m - *k;
		    MLFORTRAN(dgemm)("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,*k+1), ldv, &WORK(1,1), 
			    ldwork, &c_b14, &C(*k+1,1), ldc
#ifdef FORTRAN_STRLEN
				     , strlen("Transpose"),strlen("Transpose")
#endif
);
		}

/*              W := W * V1 */

		MLFORTRAN(dtrmm)("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &V(1,1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen("No transpose"),strlen("Unit")
#endif
);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    for (i = 1; i <= *n; ++i) {
			C(j,i) -= WORK(i,j);
/* L140: */
		    }
/* L150: */
		}

	    } else if (MLFORTRAN(lsame)(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in 
WORK)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    MLFORTRAN(dcopy)(m, &C(1,j), &c__1, &WORK(1,j), &c__1);
/* L160: */
		}

/*              W := W * V1' */

		MLFORTRAN(dtrmm)("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			V(1,1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen("Transpose"),strlen("Unit")
#endif
);
		if (*n > *k) {

/*                 W := W + C2 * V2' */

		    i__1 = *n - *k;
		    MLFORTRAN(dgemm)("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    C(1,*k+1), ldc, &V(1,*k+1), ldv, &c_b14, &WORK(1,1), 
			    ldwork
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("Transpose")
#endif
);
		}

/*              W := W * T  or  W * T' */

		MLFORTRAN(dtrmm)("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen(trans),strlen("Non-unit")
#endif
);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

		    i__1 = *n - *k;
		    MLFORTRAN(dgemm)("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &WORK(1,1), ldwork, &V(1,*k+1), ldv, &c_b14, &C(1,*k+1), ldc
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("No transpose")
#endif
);
		}

/*              W := W * V1 */

		MLFORTRAN(dtrmm)("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &V(1,1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Upper"),strlen("No transpose"),strlen("Unit")
#endif
);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    for (i = 1; i <= *m; ++i) {
			C(i,j) -= WORK(i,j);
/* L170: */
		    }
/* L180: */
		}

	    }

	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns)   
             where  V2  is unit lower triangular. */

	    if (MLFORTRAN(lsame)(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored i
n WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    MLFORTRAN(dcopy)(n, &C(*m-*k+j,1), ldc, &WORK(1,j), &c__1);
/* L190: */
		}

/*              W := W * V2' */

		MLFORTRAN(dtrmm)("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			V(1,*m-*k+1), ldv, &WORK(1,1)
			, ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen("Transpose"),strlen("Unit")
#endif
);
		if (*m > *k) {

/*                 W := W + C1'*V1' */

		    i__1 = *m - *k;
		    MLFORTRAN(dgemm)("Transpose", "Transpose", n, k, &i__1, &c_b14, &C(1,1), ldc, &V(1,1), ldv, &c_b14, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				     , strlen("Transpose"),strlen("Transpose")
#endif
);
		}

/*              W := W * T'  or  W * T */

		MLFORTRAN(dtrmm)("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen(transt),strlen("Non-unit")
#endif
);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1' * W' */

		    i__1 = *m - *k;
		    MLFORTRAN(dgemm)("Transpose", "Transpose", &i__1, n, k, &c_b25, &V(1,1), ldv, &WORK(1,1), ldwork, &
			    c_b14, &C(1,1), ldc
#ifdef FORTRAN_STRLEN
				     , strlen("Transpose"),strlen("Transpose")
#endif
);
		}

/*              W := W * V2 */

		MLFORTRAN(dtrmm)("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &V(1,*m-*k+1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen("No transpose"),strlen("Unit")
#endif
);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    for (i = 1; i <= *n; ++i) {
			C(*m-*k+j,i) -= WORK(i,j)
				;
/* L200: */
		    }
/* L210: */
		}

	    } else if (MLFORTRAN(lsame)(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in 
WORK)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    MLFORTRAN(dcopy)(m, &C(1,*n-*k+j), &c__1, &WORK(1,j), &c__1);
/* L220: */
		}

/*              W := W * V2' */

		MLFORTRAN(dtrmm)("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			V(1,*n-*k+1), ldv, &WORK(1,1)
			, ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen("Transpose"),strlen("Unit")
#endif
);
		if (*n > *k) {

/*                 W := W + C1 * V1' */

		    i__1 = *n - *k;
		    MLFORTRAN(dgemm)("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    C(1,1), ldc, &V(1,1), ldv, &c_b14, &
			    WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("Transpose")
#endif
);
		}

/*              W := W * T  or  W * T' */

		MLFORTRAN(dtrmm)("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &T(1,1), ldt, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen(trans),strlen("Non-unit")
#endif
);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

		    i__1 = *n - *k;
		    MLFORTRAN(dgemm)("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &WORK(1,1), ldwork, &V(1,1), 
			    ldv, &c_b14, &C(1,1), ldc
#ifdef FORTRAN_STRLEN
				     , strlen("No transpose"),strlen("No transpose")
#endif
);
		}

/*              W := W * V2 */

		MLFORTRAN(dtrmm)("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &V(1,*n-*k+1), ldv, &WORK(1,1), ldwork
#ifdef FORTRAN_STRLEN
				 , strlen("Right"),strlen("Lower"),strlen("No transpose"),strlen("Unit")
#endif
);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= *k; ++j) {
		    for (i = 1; i <= *m; ++i) {
			C(i,*n-*k+j) -= WORK(i,j);
/* L230: */
		    }
/* L240: */
		}

	    }

	}
    }

    return 0;

/*     End of DLARFB */

} /* dlarfb_ */
#endif


/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifndef ML_DCOPY_FUNC
/* Subroutine */ int MLFORTRAN(dcopy)(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{


    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    for (i = 1; i <= *n; ++i) {
	DY(iy) = DX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    for (i = 1; i <= m; ++i) {
	DY(i) = DX(i);
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    for (i = mp1; i <= *n; i += 7) {
	DY(i) = DX(i);
	DY(i + 1) = DX(i + 1);
	DY(i + 2) = DX(i + 2);
	DY(i + 3) = DX(i + 3);
	DY(i + 4) = DX(i + 4);
	DY(i + 5) = DX(i + 5);
	DY(i + 6) = DX(i + 6);
/* L50: */
    }
    return 0;
} /* dcopy_ */
#endif


/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifndef ML_DTRMM_FUNC
/* Subroutine */ int MLFORTRAN(dtrmm)(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, doublereal *alpha, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, int dummy1, int dummy2, 
int dummy3, int dummy4)
{


    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i, j, k;
    static logical lside;
    static integer nrowa;
    static logical upper;
    static logical nounit;


/*  Purpose   
    =======   

    DTRMM  performs one of the matrix-matrix operations   

       B := alpha*op( A )*B,   or   B := alpha*B*op( A ),   

    where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or 
  
    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of 
  

       op( A ) = A   or   op( A ) = A'.   

    Parameters   
    ==========   

    SIDE   - CHARACTER*1.   
             On entry,  SIDE specifies whether  op( A ) multiplies B from 
  
             the left or right as follows:   

                SIDE = 'L' or 'l'   B := alpha*op( A )*B.   

                SIDE = 'R' or 'r'   B := alpha*B*op( A ).   

             Unchanged on exit.   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix A is an upper or 
  
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in 
  
             the matrix multiplication as follows:   

                TRANSA = 'N' or 'n'   op( A ) = A.   

                TRANSA = 'T' or 't'   op( A ) = A'.   

                TRANSA = 'C' or 'c'   op( A ) = A'.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit triangular 
  
             as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of B. M must be at 
  
             least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of B.  N must be 
  
             at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry,  ALPHA specifies the scalar  alpha. When  alpha is 
  
             zero then  A is not referenced and  B need not be set before 
  
             entry.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m 
  
             when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. 
  
             Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k 
  
             upper triangular part of the array  A must contain the upper 
  
             triangular matrix  and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k 
  
             lower triangular part of the array  A must contain the lower 
  
             triangular matrix  and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u',  the diagonal elements of 
  
             A  are not referenced either,  but are assumed to be  unity. 
  
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then 
  
             LDA  must be at least  ML_max( 1, m ),  when  SIDE = 'R' or 'r' 
  
             then LDA must be at least ML_max( 1, n ).   
             Unchanged on exit.   

    B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
             Before entry,  the leading  m by n part of the array  B must 
  
             contain the matrix  B,  and  on exit  is overwritten  by the 
  
             transformed matrix.   

    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared 
  
             in  the  calling  (sub)  program.   LDB  must  be  at  least 
  
             ML_max( 1, m ).   
             Unchanged on exit.   


    Level 3 Blas routine.   

    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    lside = MLFORTRAN(lsame)(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = MLFORTRAN(lsame)(diag, "N");
    upper = MLFORTRAN(lsame)(uplo, "U");

    info = 0;
    if (! lside && ! MLFORTRAN(lsame)(side, "R")) {
	info = 1;
    } else if (! upper && ! MLFORTRAN(lsame)(uplo, "L")) {
	info = 2;
    } else if (! MLFORTRAN(lsame)(transa, "N") && ! MLFORTRAN(lsame)(transa, "T") 
	    && ! MLFORTRAN(lsame)(transa, "C")) {
	info = 3;
    } else if (! MLFORTRAN(lsame)(diag, "U") && ! MLFORTRAN(lsame)(diag, "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < ML_max(1,nrowa)) {
	info = 9;
    } else if (*ldb < ML_max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	MLFORTRAN(xerbla)("DTRMM ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	for (j = 1; j <= *n; ++j) {
	    for (i = 1; i <= *m; ++i) {
		B(i,j) = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }

/*     Start the operations. */

    if (lside) {
	if (MLFORTRAN(lsame)(transa, "N")) {

/*           Form  B := alpha*A*B. */

	    if (upper) {
		for (j = 1; j <= *n; ++j) {
		    for (k = 1; k <= *m; ++k) {
			if (B(k,j) != 0.) {
			    temp = *alpha * B(k,j);
			    for (i = 1; i <= k-1; ++i) {
				B(i,j) += temp * A(i,k);
/* L30: */
			    }
			    if (nounit) {
				temp *= A(k,k);
			    }
			    B(k,j) = temp;
			}
/* L40: */
		    }
/* L50: */
		}
	    } else {
		for (j = 1; j <= *n; ++j) {
		    for (k = *m; k >= 1; --k) {
			if (B(k,j) != 0.) {
			    temp = *alpha * B(k,j);
			    B(k,j) = temp;
			    if (nounit) {
				B(k,j) *= A(k,k);
			    }
			    for (i = k + 1; i <= *m; ++i) {
				B(i,j) += temp * A(i,k);
/* L60: */
			    }
			}
/* L70: */
		    }
/* L80: */
		}
	    }
	} else {

/*           Form  B := alpha*B*A'. */

	    if (upper) {
		for (j = 1; j <= *n; ++j) {
		    for (i = *m; i >= 1; --i) {
			temp = B(i,j);
			if (nounit) {
			    temp *= A(i,i);
			}
			for (k = 1; k <= i-1; ++k) {
			    temp += A(k,i) * B(k,j);
/* L90: */
			}
			B(i,j) = *alpha * temp;
/* L100: */
		    }
/* L110: */
		}
	    } else {
		for (j = 1; j <= *n; ++j) {
		    for (i = 1; i <= *m; ++i) {
			temp = B(i,j);
			if (nounit) {
			    temp *= A(i,i);
			}
			for (k = i + 1; k <= *m; ++k) {
			    temp += A(k,i) * B(k,j);
/* L120: */
			}
			B(i,j) = *alpha * temp;
/* L130: */
		    }
/* L140: */
		}
	    }
	}
    } else {
	if (MLFORTRAN(lsame)(transa, "N")) {

/*           Form  B := alpha*B*A. */

	    if (upper) {
		for (j = *n; j >= 1; --j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= A(j,j);
		    }
		    for (i = 1; i <= *m; ++i) {
			B(i,j) = temp * B(i,j);
/* L150: */
		    }
		    for (k = 1; k <= j-1; ++k) {
			if (A(k,j) != 0.) {
			    temp = *alpha * A(k,j);
			    for (i = 1; i <= *m; ++i) {
				B(i,j) += temp * B(i,k);
/* L160: */
			    }
			}
/* L170: */
		    }
/* L180: */
		}
	    } else {
		for (j = 1; j <= *n; ++j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= A(j,j);
		    }
		    for (i = 1; i <= *m; ++i) {
			B(i,j) = temp * B(i,j);
/* L190: */
		    }
		    for (k = j + 1; k <= *n; ++k) {
			if (A(k,j) != 0.) {
			    temp = *alpha * A(k,j);
			    for (i = 1; i <= *m; ++i) {
				B(i,j) += temp * B(i,k);
/* L200: */
			    }
			}
/* L210: */
		    }
/* L220: */
		}
	    }
	} else {

/*           Form  B := alpha*B*A'. */

	    if (upper) {
		for (k = 1; k <= *n; ++k) {
		    for (j = 1; j <= k-1; ++j) {
			if (A(j,k) != 0.) {
			    temp = *alpha * A(j,k);
			    for (i = 1; i <= *m; ++i) {
				B(i,j) += temp * B(i,k);
/* L230: */
			    }
			}
/* L240: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= A(k,k);
		    }
		    if (temp != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = temp * B(i,k);
/* L250: */
			}
		    }
/* L260: */
		}
	    } else {
		for (k = *n; k >= 1; --k) {
		    for (j = k + 1; j <= *n; ++j) {
			if (A(j,k) != 0.) {
			    temp = *alpha * A(j,k);
			    for (i = 1; i <= *m; ++i) {
				B(i,j) += temp * B(i,k);
/* L270: */
			    }
			}
/* L280: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= A(k,k);
		    }
		    if (temp != 1.) {
			for (i = 1; i <= *m; ++i) {
			    B(i,k) = temp * B(i,k);
/* L290: */
			}
		    }
/* L300: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRMM . */

} /* dtrmm_ */
#endif


/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifndef ML_DGEMV_FUNC
/* Subroutine */ int MLFORTRAN(dgemv)(char *trans, integer *m, integer *n, doublereal *
	alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy, int dummy1)
{


    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer lenx, leny, i, j;
    static integer ix, iy, jx, jy, kx, ky;


/*  Purpose   
    =======   

    DGEMV  performs one of the matrix-vector operations   

       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   

    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   

    Parameters   
    ==========   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             ML_max( 1, m ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( n - 1 )*ML_abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*ML_abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( m - 1 )*ML_abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*ML_abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the 
  
             updated vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! MLFORTRAN(lsame)(trans, "N") && ! MLFORTRAN(lsame)(trans, "T") && ! 
	    MLFORTRAN(lsame)(trans, "C")) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < ML_max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	MLFORTRAN(xerbla)("DGEMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0. && *beta == 1.)) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  
       up the start points in  X  and  Y. */

    if (MLFORTRAN(lsame)(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   

       First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		for (i = 1; i <= leny; ++i) {
		    Y(i) = 0.;
/* L10: */
		}
	    } else {
		for (i = 1; i <= leny; ++i) {
		    Y(i) = *beta * Y(i);
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		for (i = 1; i <= leny; ++i) {
		    Y(iy) = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		for (i = 1; i <= leny; ++i) {
		    Y(iy) = *beta * Y(iy);
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (MLFORTRAN(lsame)(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.) {
		    temp = *alpha * X(jx);
		    for (i = 1; i <= *m; ++i) {
			Y(i) += temp * A(i,j);
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.) {
		    temp = *alpha * X(jx);
		    iy = ky;
		    for (i = 1; i <= *m; ++i) {
			Y(iy) += temp * A(i,j);
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    for (j = 1; j <= *n; ++j) {
		temp = 0.;
		for (i = 1; i <= *m; ++i) {
		    temp += A(i,j) * X(i);
/* L90: */
		}
		Y(jy) += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    for (j = 1; j <= *n; ++j) {
		temp = 0.;
		ix = kx;
		for (i = 1; i <= *m; ++i) {
		    temp += A(i,j) * X(ix);
		    ix += *incx;
/* L110: */
		}
		Y(jy) += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return 0;

/*     End of DGEMV . */

} /* dgemv_ */
#endif


/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifndef ML_DTRMV_FUNC
/* Subroutine */ int MLFORTRAN(dtrmv)(char *uplo, char *trans, char *diag, integer *n, 
	doublereal *a, integer *lda, doublereal *x, integer *incx)
{


    /* System generated locals */

    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i, j;
    static integer ix, jx, kx;
    static logical nounit;


/*  Purpose   
    =======   

    DTRMV  performs one of the matrix-vector operations   

       x := A*x,   or   x := A'*x,   

    where x is an n element vector and  A is an n by n unit, or non-unit, 
  
    upper or lower triangular matrix.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   x := A*x.   

                TRANS = 'T' or 't'   x := A'*x.   

                TRANS = 'C' or 'c'   x := A'*x.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular matrix and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular matrix and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             ML_max( 1, n ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*ML_abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x. On exit, X is overwritten with the   
             tranformed vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! MLFORTRAN(lsame)(uplo, "U") && ! MLFORTRAN(lsame)(uplo, "L")) {
	info = 1;
    } else if (! MLFORTRAN(lsame)(trans, "N") && ! MLFORTRAN(lsame)(trans, "T") &&
	     ! MLFORTRAN(lsame)(trans, "C")) {
	info = 2;
    } else if (! MLFORTRAN(lsame)(diag, "U") && ! MLFORTRAN(lsame)(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < ML_max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	MLFORTRAN(xerbla)("DTRMV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = MLFORTRAN(lsame)(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (MLFORTRAN(lsame)(trans, "N")) {

/*        Form  x := A*x. */

	if (MLFORTRAN(lsame)(uplo, "U")) {
	    if (*incx == 1) {
		for (j = 1; j <= *n; ++j) {
		    if (X(j) != 0.) {
			temp = X(j);
			for (i = 1; i <= j-1; ++i) {
			    X(i) += temp * A(i,j);
/* L10: */
			}
			if (nounit) {
			    X(j) *= A(j,j);
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx;
		for (j = 1; j <= *n; ++j) {
		    if (X(jx) != 0.) {
			temp = X(jx);
			ix = kx;
			for (i = 1; i <= j-1; ++i) {
			    X(ix) += temp * A(i,j);
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    X(jx) *= A(j,j);
			}
		    }
		    jx += *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (X(j) != 0.) {
			temp = X(j);
			for (i = *n; i >= j+1; --i) {
			    X(i) += temp * A(i,j);
/* L50: */
			}
			if (nounit) {
			    X(j) *= A(j,j);
			}
		    }
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (X(jx) != 0.) {
			temp = X(jx);
			ix = kx;
			for (i = *n; i >= j+1; --i) {
			    X(ix) += temp * A(i,j);
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    X(jx) *= A(j,j);
			}
		    }
		    jx -= *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := A'*x. */

	if (MLFORTRAN(lsame)(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = X(j);
		    if (nounit) {
			temp *= A(j,j);
		    }
		    for (i = j - 1; i >= 1; --i) {
			temp += A(i,j) * X(i);
/* L90: */
		    }
		    X(j) = temp;
/* L100: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    temp = X(jx);
		    ix = jx;
		    if (nounit) {
			temp *= A(j,j);
		    }
		    for (i = j - 1; i >= 1; --i) {
			ix -= *incx;
			temp += A(i,j) * X(ix);
/* L110: */
		    }
		    X(jx) = temp;
		    jx -= *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = 1; j <= *n; ++j) {
		    temp = X(j);
		    if (nounit) {
			temp *= A(j,j);
		    }
		    for (i = j + 1; i <= *n; ++i) {
			temp += A(i,j) * X(i);
/* L130: */
		    }
		    X(j) = temp;
/* L140: */
		}
	    } else {
		jx = kx;
		for (j = 1; j <= *n; ++j) {
		    temp = X(jx);
		    ix = jx;
		    if (nounit) {
			temp *= A(j,j);
		    }
		    for (i = j + 1; i <= *n; ++i) {
			ix += *incx;
			temp += A(i,j) * X(ix);
/* L150: */
		    }
		    X(jx) = temp;
		    jx += *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of DTRMV . */

} /* dtrmv_ */
#endif


#ifndef ML_DLARF_FUNC
/* Subroutine */ int MLFORTRAN(dlarf)(char *side, integer *m, integer *n, doublereal *v,
	 integer *incv, doublereal *tau, doublereal *c, integer *ldc, 
	doublereal *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARF applies a real elementary reflector H to a real m by n matrix   
    C, from either the left or the right. H is represented in the form   

          H = I - tau * v * v'   

    where tau is a real scalar and v is a real vector.   

    If tau = 0, then H is taken to be the unit matrix.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': form  H * C   
            = 'R': form  C * H   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    V       (input) DOUBLE PRECISION array, dimension   
                       (1 + (M-1)*ML_abs(INCV)) if SIDE = 'L'   
                    or (1 + (N-1)*ML_abs(INCV)) if SIDE = 'R'   
            The vector v in the representation of H. V is not used if   
            TAU = 0.   

    INCV    (input) INTEGER   
            The increment between elements of v. INCV <> 0.   

    TAU     (input) DOUBLE PRECISION   
            The value tau in the representation of H.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by the matrix H * C if SIDE = 'L', 
  
            or C * H if SIDE = 'R'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= ML_max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension   
                           (N) if SIDE = 'L'   
                        or (M) if SIDE = 'R'   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b4 = 1.;
    static doublereal c_b5 = 0.;
    static integer c__1 = 1;
    
    /* System generated locals */
    doublereal d__1;
    /* Local variables */


#undef V
#undef WORK
#define V(I) v[(I)-1]
#define WORK(I) work[(I)-1]

#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (MLFORTRAN(lsame)(side, "L")) {

/*        Form  H * C */

	if (*tau != 0.) {

/*           w := C' * v */

	    MLFORTRAN(dgemv)("Transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), incv, &
		    c_b5, &WORK(1), &c__1
#ifdef FORTRAN_STRLEN
			     , strlen("Transpose")
#endif
			     );

/*           C := C - v * w' */

	    d__1 = -(*tau);
	    MLFORTRAN(dger)(m, n, &d__1, &V(1), incv, &WORK(1), &c__1, &C(1,1), 
		    ldc);
	}
    } else {

/*        Form  C * H */

	if (*tau != 0.) {

/*           w := C * v */

	    MLFORTRAN(dgemv)("No transpose", m, n, &c_b4, &C(1,1), ldc, &V(1), 
		    incv, &c_b5, &WORK(1), &c__1
#ifdef FORTRAN_STRLEN
			     , strlen("No transpose")
#endif
			     );

/*           C := C - w * v' */

	    d__1 = -(*tau);
	    MLFORTRAN(dger)(m, n, &d__1, &WORK(1), &c__1, &V(1), incv, &C(1,1), 
		    ldc);
	}
    }
    return 0;

/*     End of DLARF */

} /* dlarf_ */
#endif


#ifndef ML_DLARFG_FUNC
/* Subroutine */ int MLFORTRAN(dlarfg)(integer *n, doublereal *alpha, doublereal *x, 
	integer *incx, doublereal *tau)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARFG generates a real elementary reflector H of order n, such   
    that   

          H * ( alpha ) = ( beta ),   H' * H = I.   
              (   x   )   (   0  )   

    where alpha and beta are scalars, and x is an (n-1)-element real   
    vector. H is represented in the form   

          H = I - tau * ( 1 ) * ( 1 v' ) ,   
                        ( v )   

    where tau is a real scalar and v is a real (n-1)-element   
    vector.   

    If the elements of x are all zero, then tau = 0 and H is taken to be 
  
    the unit matrix.   

    Otherwise  1 <= tau <= 2.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the elementary reflector.   

    ALPHA   (input/output) DOUBLE PRECISION   
            On entry, the value alpha.   
            On exit, it is overwritten with the value beta.   

    X       (input/output) DOUBLE PRECISION array, dimension   
                           (1+(N-2)*ML_abs(INCX))   
            On entry, the vector x.   
            On exit, it is overwritten with the vector v.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    TAU     (output) DOUBLE PRECISION   
            The value tau.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    /* Builtin functions */
    /* Local variables */
    static doublereal beta;
    static integer j;
    static doublereal xnorm;
    static doublereal safmin, rsafmn;
    static integer knt;


#define X(I) x[(I)-1]


    if (*n <= 1) {
	*tau = 0.;
	return 0;
    }

    i__1 = *n - 1;
    xnorm = MLFORTRAN(dnrm2)(&i__1, &X(1), incx);

    if (xnorm == 0.) {

/*        H  =  I */

	*tau = 0.;
    } else {

/*        general case */

	d__1 = MLFORTRAN(dlapy2)(alpha, &xnorm);
	beta = -ml_d_sign(&d__1, alpha);
	safmin = MLFORTRAN(dlamch)("S") / MLFORTRAN(dlamch)("E");
	if (ML_abs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute 
them */

	    rsafmn = 1. / safmin;
	    knt = 0;
L10:
	    ++knt;
	    i__1 = *n - 1;
	    MLFORTRAN(dscal)(&i__1, &rsafmn, &X(1), incx);
	    beta *= rsafmn;
	    *alpha *= rsafmn;
	    if (ML_abs(beta) < safmin) {
		goto L10;
	    }

/*           New BETA is at most 1, at least SAFMIN */

	    i__1 = *n - 1;
	    xnorm = MLFORTRAN(dnrm2)(&i__1, &X(1), incx);
	    d__1 = MLFORTRAN(dlapy2)(alpha, &xnorm);
	    beta = -ml_d_sign(&d__1, alpha);
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);
	    MLFORTRAN(dscal)(&i__1, &d__1, &X(1), incx);

/*           If ALPHA is subnormal, it may lose relative accuracy 
*/

	    *alpha = beta;
	    i__1 = knt;
	    for (j = 1; j <= knt; ++j) {
		*alpha *= safmin;
/* L20: */
	    }
	} else {
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);
	    MLFORTRAN(dscal)(&i__1, &d__1, &X(1), incx);
	    *alpha = beta;
	}
    }

    return 0;

/*     End of DLARFG */

} /* dlarfg_ */
#endif


/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


#ifndef ML_DNRM2_FUNC
doublereal MLFORTRAN(dnrm2)(integer *n, doublereal *x, integer *incx)
{


    doublereal ret_val, d__1;


    /* Local variables */
    static doublereal norm, scale, absxi;
    static integer ix;
    static doublereal ssq;


/*  DNRM2 returns the euclidean norm of a vector via the function   
    name, so that   

       DNRM2 := sqrt( x'*x )   



    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to DLASSQ.   
       Sven Hammarling, Nag Ltd.   


    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]


    if (*n < 1 || *incx < 1) {
	norm = 0.;
    } else if (*n == 1) {
	norm = ML_abs(X(1));
    } else {
	scale = 0.;
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK 
  
          auxiliary routine:   
          CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */

	for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
	    if (X(ix) != 0.) {
		absxi = (d__1 = X(ix), ML_abs(d__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    d__1 = scale / absxi;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of DNRM2. */

} /* dnrm2_ */
#endif


#ifndef ML_DLAMCH_FUNC
doublereal MLFORTRAN(dlamch)(char *cmach)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMCH determines double precision machine parameters.   

    Arguments   
    =========   

    CMACH   (input) CHARACTER*1   
            Specifies the value to be returned by DLAMCH:   
            = 'E' or 'e',   DLAMCH := eps   
            = 'S' or 's ,   DLAMCH := sfmin   
            = 'B' or 'b',   DLAMCH := base   
            = 'P' or 'p',   DLAMCH := eps*base   
            = 'N' or 'n',   DLAMCH := t   
            = 'R' or 'r',   DLAMCH := rnd   
            = 'M' or 'm',   DLAMCH := emin   
            = 'U' or 'u',   DLAMCH := rmin   
            = 'L' or 'l',   DLAMCH := emax   
            = 'O' or 'o',   DLAMCH := rmax   

            where   

            eps   = relative machine precision   
            sfmin = safe minimum, such that 1/sfmin does not overflow   
            base  = base of the machine   
            prec  = eps*base   
            t     = number of (base) digits in the mantissa   
            rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
            emin  = minimum exponent before (gradual) underflow   
            rmin  = underflow threshold - base**(emin-1)   
            emax  = largest exponent before overflow   
            rmax  = overflow threshold  - (base**emax)*(1-eps)   

   ===================================================================== 
*/
/* >>Start of File<<   
       Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    doublereal ret_val;
    /* Local variables */
    static doublereal base;
    static integer beta;
    static doublereal emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static doublereal rmin, rmax, t, rmach;
    static doublereal small, sfmin;
    static integer it;
    static doublereal rnd, eps;



    if (first) {
	first = FALSE_;
	MLFORTRAN(dlamc2)(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (doublereal) beta;
	t = (doublereal) it;
	if (lrnd) {
	    rnd = 1.;
	    i__1 = 1 - it;
	    eps = ml_pow_di(&base, &i__1) / 2;
	} else {
	    rnd = 0.;
	    i__1 = 1 - it;
	    eps = ml_pow_di(&base, &i__1);
	}
	prec = eps * base;
	emin = (doublereal) imin;
	emax = (doublereal) imax;
	sfmin = rmin;
	small = 1. / rmax;
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rou
nding   
             causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.);
	}
    }

    if (MLFORTRAN(lsame)(cmach, "E")) {
	rmach = eps;
    } else if (MLFORTRAN(lsame)(cmach, "S")) {
	rmach = sfmin;
    } else if (MLFORTRAN(lsame)(cmach, "B")) {
	rmach = base;
    } else if (MLFORTRAN(lsame)(cmach, "P")) {
	rmach = prec;
    } else if (MLFORTRAN(lsame)(cmach, "N")) {
	rmach = t;
    } else if (MLFORTRAN(lsame)(cmach, "R")) {
	rmach = rnd;
    } else if (MLFORTRAN(lsame)(cmach, "M")) {
	rmach = emin;
    } else if (MLFORTRAN(lsame)(cmach, "U")) {
	rmach = rmin;
    } else if (MLFORTRAN(lsame)(cmach, "L")) {
	rmach = emax;
    } else if (MLFORTRAN(lsame)(cmach, "O")) {
	rmach = rmax;
    }

    ret_val = rmach;
    return ret_val;

/*     End of DLAMCH */

} /* dlamch_ */
#endif


#ifndef ML_DLAMC1_FUNC
/* Subroutine */ int MLFORTRAN(dlamc1)(integer *beta, integer *t, logical *rnd, logical 
	*ieee1)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC1 determines the machine parameters given by BETA, T, RND, and   
    IEEE1.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    IEEE1   (output) LOGICAL   
            Specifies whether rounding appears to be done in the IEEE   
            'round to nearest' style.   

    Further Details   
    ===============   

    The routine is based on the routine  ENVRON  by Malcolm and   
    incorporates suggestions by Gentleman and Marovich. See   

       Malcolm M. A. (1972) Algorithms to reveal properties of   
          floating-point arithmetic. Comms. of the ACM, 15, 949-951.   

       Gentleman W. M. and Marovich S. B. (1974) More on algorithms   
          that reveal properties of floating point arithmetic units.   
          Comms. of the ACM, 17, 276-277.   

   ===================================================================== 
*/
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    doublereal d__1, d__2;
    /* Local variables */
    static logical lrnd;
    static doublereal a, b, c, f;
    static integer lbeta;
    static doublereal savec;
    static logical lieee1;
    static doublereal t1, t2;
    static integer lt;
    static doublereal one, qtr;



    if (first) {
	first = FALSE_;
	one = 1.;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BE
TA,   
          IEEE1, T and RND.   

          Throughout this routine  we use the function  DLAMC3  to ens
ure   
          that relevant values are  stored and not held in registers, 
 or   
          are not affected by optimizers.   

          Compute  a = 2.0**m  with the  smallest positive integer m s
uch   
          that   

             fl( a + 1.0 ) = a. */

	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L10:
	if (c == one) {
	    a *= 2;
	    c = MLFORTRAN(dlamc3)(&a, &one);
	    d__1 = -a;
	    c = MLFORTRAN(dlamc3)(&c, &d__1);
	    goto L10;
	}
/* +       END WHILE   

          Now compute  b = 2.0**m  with the smallest positive integer 
m   
          such that   

             fl( a + b ) .gt. a. */

	b = 1.;
	c = MLFORTRAN(dlamc3)(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
L20:
	qtr =  MLFORTRAN(dlamc3)(&c, &a); /* needed for tflop compiler */
	if (c == a) {
	    b *= 2;
	    c = MLFORTRAN(dlamc3)(&a, &b);
	    goto L20;
	}
/* +       END WHILE   

          Now compute the base.  a and c  are neighbouring floating po
int   
          numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and
 so   
          their difference is beta. Adding 0.25 to c is to ensure that
 it   
          is truncated to beta and not ( beta - 1 ). */

	qtr = one / 4;
	savec = c;
	d__1 = -a;
	c = MLFORTRAN(dlamc3)(&c, &d__1);
	lbeta = (integer) (c + qtr);

/*        Now determine whether rounding or chopping occurs,  by addin
g a   
          bit  less  than  beta/2  and a  bit  more  than  beta/2  to 
 a. */

	b = (doublereal) lbeta;
	d__1 = b / 2;
	d__2 = -b / 100;
	f = MLFORTRAN(dlamc3)(&d__1, &d__2);
	c = MLFORTRAN(dlamc3)(&f, &a);
	if (c == a) {
	    lrnd = TRUE_;
	} else {
	    lrnd = FALSE_;
	}
	d__1 = b / 2;
	d__2 = b / 100;
	f = MLFORTRAN(dlamc3)(&d__1, &d__2);
	c = MLFORTRAN(dlamc3)(&f, &a);
	if (lrnd && c == a) {
	    lrnd = FALSE_;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round
 to   
          nearest' style. B/2 is half a unit in the last place of the 
two   
          numbers A and SAVEC. Furthermore, A is even, i.e. has last  
bit   
          zero, and SAVEC is odd. Thus adding B/2 to A should not  cha
nge   
          A, but adding B/2 to SAVEC should change SAVEC. */

	d__1 = b / 2;
	t1 = MLFORTRAN(dlamc3)(&d__1, &a);
	d__1 = b / 2;
	t2 = MLFORTRAN(dlamc3)(&d__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  integer part
 of   
          log to the base beta of a,  however it is safer to determine
  t   
          by powering.  So we find t as the smallest positive integer 
for   
          which   

             fl( beta**t + 1.0 ) = 1.0. */

	lt = 0;
	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L30:
	if (c == one) {
	    ++lt;
	    a *= lbeta;
	    c = MLFORTRAN(dlamc3)(&a, &one);
	    d__1 = -a;
	    c = MLFORTRAN(dlamc3)(&c, &d__1);
	    goto L30;
	}
/* +       END WHILE */

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return 0;

/*     End of DLAMC1 */

} /* dlamc1_ */
#endif


#ifndef ML_DLAMC2_FUNC
/* Subroutine */ int MLFORTRAN(dlamc2)(integer *beta, integer *t, logical *rnd, 
	doublereal *eps, integer *emin, doublereal *rmin, integer *emax, 
	doublereal *rmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC2 determines the machine parameters specified in its argument   
    list.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    EPS     (output) DOUBLE PRECISION   
            The smallest positive number such that   

               fl( 1.0 - EPS ) .LT. 1.0,   

            where fl denotes the computed value.   

    EMIN    (output) INTEGER   
            The minimum exponent before (gradual) underflow occurs.   

    RMIN    (output) DOUBLE PRECISION   
            The smallest normalized number for the machine, given by   
            BASE**( EMIN - 1 ), where  BASE  is the floating point value 
  
            of BETA.   

    EMAX    (output) INTEGER   
            The maximum exponent before overflow occurs.   

    RMAX    (output) DOUBLE PRECISION   
            The largest positive number for the machine, given by   
            BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point 
  
            value of BETA.   

    Further Details   
    ===============   

    The computation of  EPS  is based on a routine PARANOIA by   
    W. Kahan of the University of California at Berkeley.   

   ===================================================================== 
*/
    /* Table of constant values */
    
    /* Initialized data */
    static logical first = TRUE_;
    static logical iwarn = FALSE_;
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;
    /* Local variables */
    static logical ieee;
    static doublereal half;
    static logical lrnd;
    static doublereal leps, zero, a, b, c;
    static integer i, lbeta;
    static doublereal rbase;
    static integer lemin, lemax, gnmin;
    static doublereal small;
    static integer gpmin;
    static doublereal third, lrmin, lrmax, sixth;
    static logical lieee1;
    static integer lt, ngnmin, ngpmin;
    static doublereal one, two;



    if (first) {
	first = FALSE_;
	zero = 0.;
	one = 1.;
	two = 2.;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values
 of   
          BETA, T, RND, EPS, EMIN and RMIN.   

          Throughout this routine  we use the function  DLAMC3  to ens
ure   
          that relevant values are stored  and not held in registers, 
 or   
          are not affected by optimizers.   

          DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. 
*/

	MLFORTRAN(dlamc1)(&lbeta, &lt, &lrnd, &lieee1);

/*        Start to find EPS. */

	b = (doublereal) lbeta;
	i__1 = -lt;
	a = ml_pow_di(&b, &i__1);
	leps = a;

/*        Try some tricks to see whether or not this is the correct  E
PS. */

	b = two / 3;
	half = one / 2;
	d__1 = -half;
	sixth = MLFORTRAN(dlamc3)(&b, &d__1);
	third = MLFORTRAN(dlamc3)(&sixth, &sixth);
	d__1 = -half;
	b = MLFORTRAN(dlamc3)(&third, &d__1);
	b = MLFORTRAN(dlamc3)(&b, &sixth);
	b = ML_abs(b);
	if (b < leps) {
	    b = leps;
	}

	leps = 1.;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    d__1 = half * leps;
/* Computing 5th power */
	    d__3 = two, d__4 = d__3, d__3 *= d__3;
/* Computing 2nd power */
	    d__5 = leps;
	    d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
	    c = MLFORTRAN(dlamc3)(&d__1, &d__2);
	    d__1 = -c;
	    c = MLFORTRAN(dlamc3)(&half, &d__1);
	    b = MLFORTRAN(dlamc3)(&half, &c);
	    d__1 = -b;
	    c = MLFORTRAN(dlamc3)(&half, &d__1);
	    b = MLFORTRAN(dlamc3)(&half, &c);
	    goto L10;
	}
/* +       END WHILE */

	if (a < leps) {
	    leps = a;
	}

/*        Computation of EPS complete.   

          Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3
)).   
          Keep dividing  A by BETA until (gradual) underflow occurs. T
his   
          is detected when we cannot recover the previous A. */

	rbase = one / lbeta;
	small = one;
	for (i = 1; i <= 3; ++i) {
	    d__1 = small * rbase;
	    small = MLFORTRAN(dlamc3)(&d__1, &zero);
/* L20: */
	}
	a = MLFORTRAN(dlamc3)(&one, &small);
	MLFORTRAN(dlamc4)(&ngpmin, &one, &lbeta);
	d__1 = -one;
	MLFORTRAN(dlamc4)(&ngnmin, &d__1, &lbeta);
	MLFORTRAN(dlamc4)(&gpmin, &a, &lbeta);
	d__1 = -a;
	MLFORTRAN(dlamc4)(&gnmin, &d__1, &lbeta);
	ieee = FALSE_;

	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual under
flow;   
                e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual und
erflow;   
                e.g., IEEE standard followers ) */
	    } else {
		lemin = ML_min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, ML_abs(i__1)) == 1) {
		lemin = ML_max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow
;   
                e.g., CYBER 205 ) */
	    } else {
		lemin = ML_min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if ((i__1 = ngpmin - ngnmin, ML_abs(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - ML_min(ngpmin,ngnmin) == 3) {
		lemin = ML_max(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflo
w;   
                no known machine ) */
	    } else {
		lemin = ML_min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else {
/* Computing MIN */
	    i__1 = ML_min(ngpmin,ngnmin), i__1 = ML_min(i__1,gpmin);
	    lemin = ML_min(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = TRUE_;
	}
/* **   
   Comment out this if block if EMIN is ok */
	if (iwarn) {
	    first = TRUE_;
	    printf("\n\n WARNING. The value EMIN may be incorrect:- ");
	    printf("EMIN = %8d\n",lemin);
	    printf("If, after inspection, the value EMIN looks acceptable");
            printf("please comment out \n the IF block as marked within the"); 
            printf("code of routine DLAMC2, \n otherwise supply EMIN"); 
            printf("explicitly.\n");
	}
/* **   

          Assume IEEE arithmetic if we found denormalised  numbers abo
ve,   
          or if arithmetic seems to round in the  IEEE style,  determi
ned   
          in routine DLAMC1. A true IEEE machine should have both  thi
ngs   
          true; however, faulty machines may have one or the other. */

	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could comp
ute   
          RMIN as BASE**( EMIN - 1 ),  but some machines underflow dur
ing   
          this computation. */

	lrmin = 1.;
	i__1 = 1 - lemin;
	for (i = 1; i <= 1-lemin; ++i) {
	    d__1 = lrmin * rbase;
	    lrmin = MLFORTRAN(dlamc3)(&d__1, &zero);
/* L30: */
	}

/*        Finally, call DLAMC5 to compute EMAX and RMAX. */

	MLFORTRAN(dlamc5)(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;

    return 0;


/*     End of DLAMC2 */

} /* dlamc2_ */
#endif


#ifndef ML_DLAMC3_FUNC
doublereal MLFORTRAN(dlamc3)(doublereal *a, doublereal *b)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC3  is intended to force  A  and  B  to be stored prior to doing 
  
    the addition of  A  and  B ,  for use in situations where optimizers 
  
    might hold one of these in a register.   

    Arguments   
    =========   

    A, B    (input) DOUBLE PRECISION   
            The values A and B.   

   ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    doublereal ret_val;



    ret_val = *a + *b;

    return ret_val;

/*     End of DLAMC3 */

} /* dlamc3_ */
#endif


#ifndef ML_DLAMC4_FUNC
/* Subroutine */ int MLFORTRAN(dlamc4)(integer *emin, doublereal *start, integer *base)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC4 is a service routine for DLAMC2.   

    Arguments   
    =========   

    EMIN    (output) EMIN   
            The minimum exponent before (gradual) underflow, computed by 
  
            setting A = START and dividing by BASE until the previous A   
            can not be recovered.   

    START   (input) DOUBLE PRECISION   
            The starting point for determining EMIN.   

    BASE    (input) INTEGER   
            The base of the machine.   

   ===================================================================== 
*/
    /* System generated locals */
    doublereal d__1;
    /* Local variables */
    static doublereal zero, a;
    static integer i;
    static doublereal rbase, b1, b2, c1, c2, d1, d2;
    static doublereal one;



    a = *start;
    one = 1.;
    rbase = one / *base;
    zero = 0.;
    *emin = 1;
    d__1 = a * rbase;
    b1 = MLFORTRAN(dlamc3)(&d__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.   
      $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	d__1 = a / *base;
	b1 = MLFORTRAN(dlamc3)(&d__1, &zero);
	d__1 = b1 * *base;
	c1 = MLFORTRAN(dlamc3)(&d__1, &zero);
	d1 = zero;
	for (i = 1; i <= *base; ++i) {
	    d1 += b1;
/* L20: */
	}
	d__1 = a * rbase;
	b2 = MLFORTRAN(dlamc3)(&d__1, &zero);
	d__1 = b2 / rbase;
	c2 = MLFORTRAN(dlamc3)(&d__1, &zero);
	d2 = zero;
	for (i = 1; i <= *base; ++i) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
/* +    END WHILE */

    return 0;

/*     End of DLAMC4 */

} /* dlamc4_ */
#endif


#ifndef ML_DLAMC5_FUNC
/* Subroutine */ int MLFORTRAN(dlamc5)(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, doublereal *rmax)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC5 attempts to compute RMAX, the largest machine floating-point   
    number, without overflow.  It assumes that EMAX + ML_abs(EMIN) sum   
    approximately to a power of 2.  It will fail on machines where this   
    assumption does not hold, for example, the Cyber 205 (EMIN = -28625, 
  
    EMAX = 28718).  It will also fail if the value supplied for EMIN is   
    too large (i.e. too close to zero), probably with overflow.   

    Arguments   
    =========   

    BETA    (input) INTEGER   
            The base of floating-point arithmetic.   

    P       (input) INTEGER   
            The number of base BETA digits in the mantissa of a   
            floating-point value.   

    EMIN    (input) INTEGER   
            The minimum exponent before (gradual) underflow.   

    IEEE    (input) LOGICAL   
            A logical flag specifying whether or not the arithmetic   
            system is thought to comply with the IEEE standard.   

    EMAX    (output) INTEGER   
            The largest exponent before overflow   

    RMAX    (output) DOUBLE PRECISION   
            The largest machine floating-point number.   

   ===================================================================== 
  


       First compute LEXP and UEXP, two powers of 2 that bound   
       ML_abs(EMIN). We then assume that EMAX + ML_abs(EMIN) will sum   
       approximately to the bound that is closest to ML_abs(EMIN).   
       (EMAX is the exponent of the required number RMAX). */
    /* Table of constant values */
    static doublereal c_b5 = 0.;
    
    /* System generated locals */
    doublereal d__1;
    /* Local variables */
    static integer lexp;
    static doublereal oldy;
    static integer uexp, i;
    static doublereal y, z;
    static integer nbits;
    static doublereal recbas;
    static integer exbits, expsum, try__;



    lexp = 1;
    exbits = 1;
L10:
    try__ = lexp << 1;
    if (try__ <= -(*emin)) {
	lexp = try__;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try__;
	++exbits;
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater   
       than or equal to EMIN. EXBITS is the number of bits needed to   
       store the exponent. */

    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }

/*     EXPSUM is the exponent range, approximately equal to   
       EMAX - EMIN + 1 . */

    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a   
       floating-point number. */

    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a   
          floating-point number, which is unlikely, or some bits are 
  
          not used in the representation of numbers, which is possible
,   
          (e.g. Cray machines) or the mantissa has an implicit bit,   
          (e.g. IEEE machines, Dec Vax machines), which is perhaps the
   
          most likely. We have to assume the last alternative.   
          If this is true, then we need to reduce EMAX by one because 
  
          there must be some way of representing zero in an implicit-b
it   
          system. On machines like Cray, we are reducing EMAX by one 
  
          unnecessarily. */

	--(*emax);
    }

    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent
   
          for infinity and NaN. */

	--(*emax);
    }

/*     Now create RMAX, the largest machine number, which should   
       be equal to (1.0 - BETA**(-P)) * BETA**EMAX .   

       First compute 1.0 - BETA**(-P), being careful that the   
       result is less than 1.0 . */

    recbas = 1. / *beta;
    z = *beta - 1.;
    y = 0.;
    for (i = 1; i <= *p; ++i) {
	z *= recbas;
	if (y < 1.) {
	    oldy = y;
	}
	y = MLFORTRAN(dlamc3)(&y, &z);
/* L20: */
    }
    if (y >= 1.) {
	y = oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

    for (i = 1; i <= *emax; ++i) {
	d__1 = y * *beta;
	y = MLFORTRAN(dlamc3)(&d__1, &c_b5);
/* L30: */
    }

    *rmax = y;
    return 0;

/*     End of DLAMC5 */

} /* dlamc5_ */
#endif

#ifndef ML_DLAPY2_FUNC
doublereal MLFORTRAN(dlapy2)(doublereal *x, doublereal *y)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary 
  
    overflow.   

    Arguments   
    =========   

    X       (input) DOUBLE PRECISION   
    Y       (input) DOUBLE PRECISION   
            X and Y specify the values x and y.   

    ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    doublereal ret_val, d__1;
    /* Local variables */
    static doublereal xabs, yabs, w, z;



    xabs = ML_abs(*x);
    yabs = ML_abs(*y);
    w = ML_max(xabs,yabs);
    z = ML_min(xabs,yabs);
    if (z == 0.) {
	ret_val = w;
    } else {
/* Computing 2nd power */
	d__1 = z / w;
	ret_val = w * sqrt(d__1 * d__1 + 1.);
    }
    return ret_val;

/*     End of DLAPY2 */

} /* dlapy2_ */
#endif

#ifndef ML_FORTRAN_SUPPORT
/* compare two strings */

integer ml_s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
register unsigned char *a, *aend, *b, *bend;
a = (unsigned char *)a0;
b = (unsigned char *)b0;
aend = a + la;
bend = b + lb;

if(la <= lb)
	{
	while(a < aend)
		if(*a != *b)
			return( *a - *b );
		else
			{ ++a; ++b; }

	while(b < bend)
		if(*b != ' ')
			return( ' ' - *b );
		else	++b;
	}

else
	{
	while(b < bend)
		if(*a == *b)
			{ ++a; ++b; }
		else
			return( *a - *b );
	while(a < aend)
		if(*a != ' ')
			return(*a - ' ');
		else	++a;
	}
return(0);
}

/* assign strings:  a = b */

int ml_s_copy(register char *a, register char *b, ftnlen la, ftnlen lb)
{
register char *aend, *bend;

aend = a + la;

if(la <= lb)
	while(a < aend)
		*a++ = *b++;

else
	{
	bend = b + lb;
	while(b < bend)
		*a++ = *b++;
	while(a < aend)
		*a++ = ' ';
	}
return 0;
}

double ml_d_sign(doublereal *a, doublereal *b)
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}

int ml_s_cat(char *lp, char *rpp[], int rnp[], int *np, ftnlen ll)
{
ftnlen i, n, nc;
char *f__rp;

n = (int)*np;
for(i = 0 ; i < n ; ++i)
        {
        nc = ll;
        if(rnp[i] < nc)
                nc = rnp[i];
        ll -= nc;
        f__rp = rpp[i];
        while(--nc >= 0)
                *lp++ = *f__rp++;
        }
while(--ll >= 0)
        *lp++ = ' ';
        return 0;
}


double ml_pow_di(doublereal *ap, integer *bp)
{
double pow, x;
integer n;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for( ; ; )
		{
		if(n & 01)
			pow *= x;
		if(n >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}
#endif

#ifndef ML_DORGQR_FUNC

/* Subroutine */ int MLFORTRAN(dorgqr)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORGQR generates an M-by-N real matrix Q with orthonormal columns,   
    which is defined as the first N columns of a product of K elementary 
  
    reflectors of order M   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by DGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by DGEQRF in the first k columns of its array   
            argument A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= ML_max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= ML_max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer   i__1, i__2, i__3;
    /* Local variables */
    static integer i, j, l, nbmin, iinfo;
    static integer ib, nb, ki, kk;
    static integer nx;
    static integer ldwork, iws;



#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < ML_max(1,*m)) {
	*info = -5;
    } else if (*lwork < ML_max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DORGQR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	WORK(1) = 1.;
	return 0;
    }

/*     Determine the block size. */

    nb = MLFORTRAN(ilaenv)(&c__1, "DORGQR", " ", m, n, k, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = MLFORTRAN(ilaenv)(&c__3, "DORGQR", " ", m, n, k, &c_n1, 6L, 1L)
		;
	nx = ML_max(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduc
e NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = MLFORTRAN(ilaenv)(&c__2, "DORGQR", " ", m, n, k, &c_n1,
			 6L, 1L);
		nbmin = ML_max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block.   
          The first kk columns are handled by the block method. */

	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
	i__1 = *k, i__2 = ki + nb;
	kk = ML_min(i__1,i__2);

/*        Set A(1:kk,kk+1:n) to zero. */

	i__1 = *n;
	for (j = kk + 1; j <= *n; ++j) {
	    i__2 = kk;
	    for (i = 1; i <= kk; ++i) {
		A(i,j) = 0.;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the last or only block. */

    if (kk < *n) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	MLFORTRAN(dorg2r)(&i__1, &i__2, &i__3, &A(kk+1,kk+1), lda, &
		TAU(kk + 1), &WORK(1), &iinfo);
    }

    if (kk > 0) {

/*        Use blocked code */

	i__1 = -nb;
	for (i = ki + 1; -nb < 0 ? i >= 1 : i <= 1; i += -nb) {
/* Computing MIN */
	    i__2 = nb, i__3 = *k - i + 1;
	    ib = ML_min(i__2,i__3);
	    if (i + ib <= *n) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__2 = *m - i + 1;
		MLFORTRAN(dlarft)("Forward", "Columnwise", &i__2, &ib, &A(i,i), lda, &TAU(i), &WORK(1), &ldwork);

/*              Apply H to A(i:m,i+ib:n) from the left */

		i__2 = *m - i + 1;
		i__3 = *n - i - ib + 1;
		MLFORTRAN(dlarfb)("Left", "No transpose", "Forward", "Columnwise", &
			i__2, &i__3, &ib, &A(i,i), lda, &WORK(1), &
			ldwork, &A(i,i+ib), lda, &WORK(ib + 1),
			 &ldwork);
	    }

/*           Apply H to rows i:m of current block */

	    i__2 = *m - i + 1;
	    MLFORTRAN(dorg2r)(&i__2, &ib, &ib, &A(i,i), lda, &TAU(i), &WORK(
		    1), &iinfo);

/*           Set rows 1:i-1 of current block to zero */

	    i__2 = i + ib - 1;
	    for (j = i; j <= i+ib-1; ++j) {
		i__3 = i - 1;
		for (l = 1; l <= i-1; ++l) {
		    A(l,j) = 0.;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    WORK(1) = (doublereal) iws;
    return 0;

/*     End of DORGQR */

} /* dorgqr_ */
#endif


#ifndef ML_DORG2R_FUNC
/* Subroutine */ int MLFORTRAN(dorg2r)(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORG2R generates an m by n real matrix Q with orthonormal columns,   
    which is defined as the first n columns of a product of k elementary 
  
    reflectors of order m   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by DGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by DGEQRF in the first k columns of its array   
            argument A.   
            On exit, the m-by-n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= ML_max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer  i__1, i__2;
    doublereal d__1;
    /* Local variables */
    static integer i, j, l;


#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < ML_max(1,*m)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DORG2R", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

/*     Initialise columns k+1:n to columns of the unit matrix */

    i__1 = *n;
    for (j = *k + 1; j <= *n; ++j) {
	i__2 = *m;
	for (l = 1; l <= *m; ++l) {
	    A(l,j) = 0.;
/* L10: */
	}
	A(j,j) = 1.;
/* L20: */
    }

    for (i = *k; i >= 1; --i) {

/*        Apply H(i) to A(i:m,i:n) from the left */

	if (i < *n) {
	    A(i,i) = 1.;
	    i__1 = *m - i + 1;
	    i__2 = *n - i;
	    MLFORTRAN(dlarf)("Left", &i__1, &i__2, &A(i,i), &c__1, &TAU(i), &
		    A(i,i+1), lda, &WORK(1));
	}
	if (i < *m) {
	    i__1 = *m - i;
	    d__1 = -TAU(i);
	    MLFORTRAN(dscal)(&i__1, &d__1, &A(i+1,i), &c__1);
	}
	A(i,i) = 1. - TAU(i);

/*        Set A(1:i-1,i) to zero */

	i__1 = i - 1;
	for (l = 1; l <= i-1; ++l) {
	    A(l,i) = 0.;
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     End of DORG2R */

} /* dorg2r_ */
#endif


#ifndef ML_DPOTRS_FUNC

/* Subroutine */ int MLFORTRAN(dpotrs)(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPOTRS solves a system of linear equations A*X = B with a symmetric   
    positive definite matrix A using the Cholesky factorization   
    A = U**T*U or A = L*L**T computed by DPOTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The triangular factor U or L from the Cholesky factorization 
  
            A = U**T*U or A = L*L**T, as computed by DPOTRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,N).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= ML_max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b9 = 1.;
    
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static logical upper;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    upper = MLFORTRAN(lsame)(uplo, "U");
    if (! upper && ! MLFORTRAN(lsame)(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < ML_max(1,*n)) {
	*info = -5;
    } else if (*ldb < ML_max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DPOTRS", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    if (upper) {

/*        Solve A*X = B where A = U'*U.   

          Solve U'*X = B, overwriting B with X. */

	MLFORTRAN(dtrsm)("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b9, &A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Upper"),strlen("Transpose"),strlen("Non-unit")
#endif
);

/*        Solve U*X = B, overwriting B with X. */

	MLFORTRAN(dtrsm)("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b9, &
		A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Upper"),strlen("No transpose"),strlen("Non-unit")
#endif
);
    } else {

/*        Solve A*X = B where A = L*L'.   

          Solve L*X = B, overwriting B with X. */

	MLFORTRAN(dtrsm)("Left", "Lower", "No transpose", "Non-unit", n, nrhs, &c_b9, &
		A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Lower"),strlen("No transpose"),strlen("Non-unit")
#endif
);

/*        Solve L'*X = B, overwriting B with X. */

	MLFORTRAN(dtrsm)("Left", "Lower", "Transpose", "Non-unit", n, nrhs, &c_b9, &A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Lower"),strlen("Transpose"),strlen("Non-unit")
#endif
);
    }

    return 0;

/*     End of DPOTRS */

} /* dpotrs_ */
#endif


#ifndef ML_DGELQ2_FUNC
/* Subroutine */ int MLFORTRAN(dgelq2)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DGELQ2 computes an LQ factorization of a real m by n matrix A:   
    A = L * Q.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix A.   
            On exit, the elements on and below the diagonal of the array 
  
            contain the m by ML_min(m,n) lower trapezoidal matrix L (L is   
            lower triangular if m <= n); the elements above the diagonal, 
  
            with the array TAU, represent the orthogonal matrix Q as a   
            product of elementary reflectors (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,M).   

    TAU     (output) DOUBLE PRECISION array, dimension (ML_min(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (M)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(k) . . . H(2) H(1), where k = ML_min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n), 
  
    and tau in TAU(i).   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1, i__2, i__3;
    /* Local variables */
    static integer i, k;
    static doublereal aii;


#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ML_max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DGELQ2", &i__1);
	return 0;
    }

    k = ML_min(*m,*n);

    i__1 = k;
    for (i = 1; i <= k; ++i) {

/*        Generate elementary reflector H(i) to annihilate A(i,i+1:n) 
*/

	i__2 = *n - i + 1;
/* Computing MIN */
	i__3 = i + 1;
	MLFORTRAN(dlarfg)(&i__2, &A(i,i), &A(i,ML_min(i+1,*n)), lda,
		 &TAU(i));
	if (i < *m) {

/*           Apply H(i) to A(i+1:m,i:n) from the right */

	    aii = A(i,i);
	    A(i,i) = 1.;
	    i__2 = *m - i;
	    i__3 = *n - i + 1;
	    MLFORTRAN(dlarf)("Right", &i__2, &i__3, &A(i,i), lda, &TAU(i), &
		    A(i+1,i), lda, &WORK(1));
	    A(i,i) = aii;
	}
/* L10: */
    }
    return 0;

/*     End of DGELQ2 */

} /* dgelq2_ */
#endif


#ifndef ML_DGELQF_FUNC
/* Subroutine */ int MLFORTRAN(dgelqf)(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGELQF computes an LQ factorization of a real M-by-N matrix A:   
    A = L * Q.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the elements on and below the diagonal of the array 
  
            contain the m-by-ML_min(m,n) lower trapezoidal matrix L (L is   
            lower triangular if m <= n); the elements above the diagonal, 
  
            with the array TAU, represent the orthogonal matrix Q as a   
            product of elementary reflectors (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,M).   

    TAU     (output) DOUBLE PRECISION array, dimension (ML_min(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= ML_max(1,M).   
            For optimum performance LWORK >= M*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(k) . . . H(2) H(1), where k = ML_min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n), 
  
    and tau in TAU(i).   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i, k, nbmin, iinfo;
    static integer ib, nb;
    static integer nx;
    static integer ldwork, iws;



#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ML_max(1,*m)) {
	*info = -4;
    } else if (*lwork < ML_max(1,*m)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DGELQF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    k = ML_min(*m,*n);
    if (k == 0) {
	WORK(1) = 1.;
	return 0;
    }

/*     Determine the block size. */

    nb = MLFORTRAN(ilaenv)(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = MLFORTRAN(ilaenv)(&c__3, "DGELQF", " ", m, n, &c_n1, &c_n1, 6L,
		 1L);
	nx = ML_max(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduc
e NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = MLFORTRAN(ilaenv)(&c__2, "DGELQF", " ", m, n, &c_n1, &
			c_n1, 6L, 1L);
		nbmin = ML_max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i = 1; nb < 0 ? i >= k-nx : i <= k-nx; i += nb) {
/* Computing MIN */
	    i__3 = k - i + 1;
	    ib = ML_min(i__3,nb);

/*           Compute the LQ factorization of the current block   
             A(i:i+ib-1,i:n) */

	    i__3 = *n - i + 1;
	    MLFORTRAN(dgelq2)(&ib, &i__3, &A(i,i), lda, &TAU(i), &WORK(1), &
		    iinfo);
	    if (i + ib <= *m) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__3 = *n - i + 1;
		MLFORTRAN(dlarft)("Forward", "Rowwise", &i__3, &ib, &A(i,i), 
			lda, &TAU(i), &WORK(1), &ldwork);

/*              Apply H to A(i+ib:m,i:n) from the right */

		i__3 = *m - i - ib + 1;
		i__4 = *n - i + 1;
		MLFORTRAN(dlarfb)("Right", "No transpose", "Forward", "Rowwise", &i__3, 
			&i__4, &ib, &A(i,i), lda, &WORK(1), &
			ldwork, &A(i+ib,i), lda, &WORK(ib + 1), &
			ldwork);
	    }
/* L10: */
	}
    } else {
	i = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i <= k) {
	i__2 = *m - i + 1;
	i__1 = *n - i + 1;
	MLFORTRAN(dgelq2)(&i__2, &i__1, &A(i,i), lda, &TAU(i), &WORK(1), &
		iinfo);
    }

    WORK(1) = (doublereal) iws;
    return 0;

/*     End of DGELQF */

} /* dgelqf_ */
#endif

#ifndef ML_DGELS_FUNC

/* Subroutine */ int MLFORTRAN(dgels)(char *trans, integer *m, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGELS solves overdetermined or underdetermined real linear systems   
    involving an M-by-N matrix A, or its transpose, using a QR or LQ   
    factorization of A.  It is assumed that A has full rank.   

    The following options are provided:   

    1. If TRANS = 'N' and m >= n:  find the least squares solution of   
       an overdetermined system, i.e., solve the least squares problem   
                    minimize || B - A*X ||.   

    2. If TRANS = 'N' and m < n:  find the minimum norm solution of   
       an underdetermined system A * X = B.   

    3. If TRANS = 'T' and m >= n:  find the minimum norm solution of   
       an undetermined system A**T * X = B.   

    4. If TRANS = 'T' and m < n:  find the least squares solution of   
       an overdetermined system, i.e., solve the least squares problem   
                    minimize || B - A**T * X ||.   

    Several right hand side vectors b and solution vectors x can be   
    handled in a single call; they are stored as the columns of the   
    M-by-NRHS right hand side matrix B and the N-by-NRHS solution   
    matrix X.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER   
            = 'N': the linear system involves A;   
            = 'T': the linear system involves A**T.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of   
            columns of the matrices B and X. NRHS >=0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit,   
              if M >= N, A is overwritten by details of its QR   
                         factorization as returned by DGEQRF;   
              if M <  N, A is overwritten by details of its LQ   
                         factorization as returned by DGELQF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,M).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the matrix B of right hand side vectors, stored   
            columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS   
            if TRANS = 'T'.   
            On exit, B is overwritten by the solution vectors, stored   
            columnwise:   
            if TRANS = 'N' and m >= n, rows 1 to n of B contain the least 
  
            squares solution vectors; the residual sum of squares for the 
  
            solution in each column is given by the sum of squares of   
            elements N+1 to M in that column;   
            if TRANS = 'N' and m < n, rows 1 to N of B contain the   
            minimum norm solution vectors;   
            if TRANS = 'T' and m >= n, rows 1 to M of B contain the   
            minimum norm solution vectors;   
            if TRANS = 'T' and m < n, rows 1 to M of B contain the   
            least squares solution vectors; the residual sum of squares   
            for the solution in each column is given by the sum of   
            squares of elements M+1 to N in that column.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= ML_MAX(1,M,N).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            LWORK >= ML_min(M,N) + ML_MAX(1,M,N,NRHS).   
            For optimal performance,   
            LWORK >= ML_min(M,N) + ML_MAX(1,M,N,NRHS) * NB   
            where NB is the optimum block size.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static doublereal c_b33 = 0.;
    static integer c__0 = 0;
    static doublereal c_b61 = 1.;
    
    /* System generated locals */
    integer i__1, i__2, i__3;
    /* Local variables */
    static doublereal anrm, bnrm;
    static integer brow;
    static logical tpsd;
    static integer i, j, iascl, ibscl;
    static integer wsize;
    static doublereal rwork[1];
    static integer nb;
    static integer mn;
    static integer scllen;
    static doublereal bignum;
    static doublereal smlnum;



#define RWORK(I) rwork[(I)]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    mn = ML_min(*m,*n);
    if (! (MLFORTRAN(lsame)(trans, "N") || MLFORTRAN(lsame)(trans, "T"))) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*lda < ML_max(1,*m)) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = ML_max(1,*m);
	if (*ldb < ML_max(i__1,*n)) {
	    *info = -8;
	} else /* if(complicated condition) */ {
/* Computing MAX   
   Computing MAX */
	    i__3 = ML_max(*m,*n);
	    i__1 = 1, i__2 = mn + ML_max(i__3,*nrhs);
	    if (*lwork < ML_max(i__1,i__2)) {
		*info = -10;
	    }
	}
    }

/*     Figure out optimal block size */

    if (*info == 0 || *info == -10) {

	tpsd = TRUE_;
	if (MLFORTRAN(lsame)(trans, "N")) {
	    tpsd = FALSE_;
	}

	if (*m >= *n) {
	    nb = MLFORTRAN(ilaenv)(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
	    if (tpsd) {
/* Computing MAX */
		i__1 = nb, i__2 = MLFORTRAN(ilaenv)(&c__1, "DORMQR", "LN", m, nrhs, n, &
			c_n1, 6L, 2L);
		nb = ML_max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = nb, i__2 = MLFORTRAN(ilaenv)(&c__1, "DORMQR", "LT", m, nrhs, n, &
			c_n1, 6L, 2L);
		nb = ML_max(i__1,i__2);
	    }
	} else {
	    nb = MLFORTRAN(ilaenv)(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
	    if (tpsd) {
/* Computing MAX */
		i__1 = nb, i__2 = MLFORTRAN(ilaenv)(&c__1, "DORMLQ", "LT", n, nrhs, m, &
			c_n1, 6L, 2L);
		nb = ML_max(i__1,i__2);
	    } else {
/* Computing MAX */
		i__1 = nb, i__2 = MLFORTRAN(ilaenv)(&c__1, "DORMLQ", "LN", n, nrhs, m, &
			c_n1, 6L, 2L);
		nb = ML_max(i__1,i__2);
	    }
	}

/* Computing MAX */
	i__1 = ML_max(*m,*n);
	wsize = mn + ML_max(i__1,*nrhs) * nb;
	WORK(1) = (doublereal) wsize;

    }

    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DGELS ", &i__1);
	return 0;
    }

/*     Quick return if possible   

   Computing MIN */
    i__1 = ML_min(*m,*n);
    if (ML_min(i__1,*nrhs) == 0) {
	i__1 = ML_max(*m,*n);
	MLFORTRAN(dlaset)("Full", &i__1, nrhs, &c_b33, &c_b33, &B(1,1), ldb);
	return 0;
    }

/*     Get machine parameters */

    smlnum = MLFORTRAN(dlamch)("S") / MLFORTRAN(dlamch)("P");
    bignum = 1. / smlnum;
    MLFORTRAN(dlabad)(&smlnum, &bignum);

/*     Scale A, B if max element outside range [SMLNUM,BIGNUM] */

    anrm = MLFORTRAN(dlange)("M", m, n, &A(1,1), lda, rwork);
    iascl = 0;
    if (anrm > 0. && anrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

	MLFORTRAN(dlascl)("G", &c__0, &c__0, &anrm, &smlnum, m, n, &A(1,1), lda, 
		info);
	iascl = 1;
    } else if (anrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

	MLFORTRAN(dlascl)("G", &c__0, &c__0, &anrm, &bignum, m, n, &A(1,1), lda, 
		info);
	iascl = 2;
    } else if (anrm == 0.) {

/*        Matrix all zero. Return zero solution. */

	i__1 = ML_max(*m,*n);
	MLFORTRAN(dlaset)("F", &i__1, nrhs, &c_b33, &c_b33, &B(1,1), ldb);
	goto L50;
    }

    brow = *m;
    if (tpsd) {
	brow = *n;
    }
    bnrm = MLFORTRAN(dlange)("M", &brow, nrhs, &B(1,1), ldb, rwork);
    ibscl = 0;
    if (bnrm > 0. && bnrm < smlnum) {

/*        Scale matrix norm up to SMLNUM */

	MLFORTRAN(dlascl)("G", &c__0, &c__0, &bnrm, &smlnum, &brow, nrhs, &B(1,1), 
		ldb, info);
	ibscl = 1;
    } else if (bnrm > bignum) {

/*        Scale matrix norm down to BIGNUM */

	MLFORTRAN(dlascl)("G", &c__0, &c__0, &bnrm, &bignum, &brow, nrhs, &B(1,1), 
		ldb, info);
	ibscl = 2;
    }

    if (*m >= *n) {

/*        compute QR factorization of A */

	i__1 = *lwork - mn;
	MLFORTRAN(dgeqrf)(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
		;

/*        workspace at least N, optimally N*NB */

	if (! tpsd) {

/*           Least-Squares Problem min || A * X - B ||   

             B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS) */

	    i__1 = *lwork - mn;
	    MLFORTRAN(dormqr)("Left", "Transpose", m, nrhs, n, &A(1,1), lda, &WORK(
		    1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info)
		    ;

/*           workspace at least NRHS, optimally NRHS*NB   

             B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS) */

	    MLFORTRAN(dtrsm)("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &
		    c_b61, &A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Upper"),strlen("No transpose"),strlen("Non-unit")
#endif
);

	    scllen = *n;

	} else {

/*           Overdetermined system of equations A' * X = B   

             B(1:N,1:NRHS) := inv(R') * B(1:N,1:NRHS) */

	    MLFORTRAN(dtrsm)("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b61, 
		    &A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Upper"),strlen("Transpose"),strlen("Non-unit")
#endif
);

/*           B(N+1:M,1:NRHS) = ZERO */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		i__2 = *m;
		for (i = *n + 1; i <= *m; ++i) {
		    B(i,j) = 0.;
/* L10: */
		}
/* L20: */
	    }

/*           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS) */

	    i__1 = *lwork - mn;
	    MLFORTRAN(dormqr)("Left", "No transpose", m, nrhs, n, &A(1,1), lda, &
		    WORK(1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB */

	    scllen = *m;

	}

    } else {

/*        Compute LQ factorization of A */

	i__1 = *lwork - mn;
	MLFORTRAN(dgelqf)(m, n, &A(1,1), lda, &WORK(1), &WORK(mn + 1), &i__1, info)
		;

/*        workspace at least M, optimally M*NB. */

	if (! tpsd) {

/*           underdetermined system of equations A * X = B   

             B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS) */

	    MLFORTRAN(dtrsm)("Left", "Lower", "No transpose", "Non-unit", m, nrhs, &
		    c_b61, &A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Lower"),strlen("No transpose"),strlen("Non-unit")
#endif
);

/*           B(M+1:N,1:NRHS) = 0 */

	    i__1 = *nrhs;
	    for (j = 1; j <= *nrhs; ++j) {
		i__2 = *n;
		for (i = *m + 1; i <= *n; ++i) {
		    B(i,j) = 0.;
/* L30: */
		}
/* L40: */
	    }

/*           B(1:N,1:NRHS) := Q(1:N,:)' * B(1:M,1:NRHS) */

	    i__1 = *lwork - mn;
	    MLFORTRAN(dormlq)("Left", "Transpose", n, nrhs, m, &A(1,1), lda, &WORK(
		    1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info)
		    ;

/*           workspace at least NRHS, optimally NRHS*NB */

	    scllen = *n;

	} else {

/*           overdetermined system min || A' * X - B ||   

             B(1:N,1:NRHS) := Q * B(1:N,1:NRHS) */

	    i__1 = *lwork - mn;
	    MLFORTRAN(dormlq)("Left", "No transpose", n, nrhs, m, &A(1,1), lda, &
		    WORK(1), &B(1,1), ldb, &WORK(mn + 1), &i__1, info);

/*           workspace at least NRHS, optimally NRHS*NB   

             B(1:M,1:NRHS) := inv(L') * B(1:M,1:NRHS) */

	    MLFORTRAN(dtrsm)("Left", "Lower", "Transpose", "Non-unit", m, nrhs, &c_b61, 
		    &A(1,1), lda, &B(1,1), ldb
#ifdef FORTRAN_STRLEN
				 , strlen("Left"),strlen("Lower"),strlen("Transpose"),strlen("Non-unit")
#endif
);

	    scllen = *m;

	}

    }

/*     Undo scaling */

    if (iascl == 1) {
	MLFORTRAN(dlascl)("G", &c__0, &c__0, &anrm, &smlnum, &scllen, nrhs, &B(1,1)
		, ldb, info);
    } else if (iascl == 2) {
	MLFORTRAN(dlascl)("G", &c__0, &c__0, &anrm, &bignum, &scllen, nrhs, &B(1,1)
		, ldb, info);
    }
    if (ibscl == 1) {
	MLFORTRAN(dlascl)("G", &c__0, &c__0, &smlnum, &bnrm, &scllen, nrhs, &B(1,1)
		, ldb, info);
    } else if (ibscl == 2) {
	MLFORTRAN(dlascl)("G", &c__0, &c__0, &bignum, &bnrm, &scllen, nrhs, &B(1,1)
		, ldb, info);
    }

L50:
    WORK(1) = (doublereal) wsize;

    return 0;

/*     End of DGELS */

} /* dgels_ */
#endif


#ifndef ML_DLASCL_FUNC

/* Subroutine */ int MLFORTRAN(dlascl)(char *type, integer *kl, integer *ku, doublereal 
	*cfrom, doublereal *cto, integer *m, integer *n, doublereal *a, 
	integer *lda, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLASCL multiplies the M by N real matrix A by the real scalar   
    CTO/CFROM.  This is done without over/underflow as long as the final 
  
    result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that 
  
    A may be full, upper triangular, lower triangular, upper Hessenberg, 
  
    or banded.   

    Arguments   
    =========   

    TYPE    (input) CHARACTER*1   
            TYPE indices the storage type of the input matrix.   
            = 'G':  A is a full matrix.   
            = 'L':  A is a lower triangular matrix.   
            = 'U':  A is an upper triangular matrix.   
            = 'H':  A is an upper Hessenberg matrix.   
            = 'B':  A is a symmetric band matrix with lower bandwidth KL 
  
                    and upper bandwidth KU and with the only the lower   
                    half stored.   
            = 'Q':  A is a symmetric band matrix with lower bandwidth KL 
  
                    and upper bandwidth KU and with the only the upper   
                    half stored.   
            = 'Z':  A is a band matrix with lower bandwidth KL and upper 
  
                    bandwidth KU.   

    KL      (input) INTEGER   
            The lower bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    KU      (input) INTEGER   
            The upper bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    CFROM   (input) DOUBLE PRECISION   
    CTO     (input) DOUBLE PRECISION   
            The matrix A is multiplied by CTO/CFROM. A(I,J) is computed   
            without over/underflow if the final result CTO*A(I,J)/CFROM   
            can be represented without over/underflow.  CFROM must be   
            nonzero.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)   
            The matrix to be multiplied by CTO/CFROM.  See TYPE for the   
            storage type.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,M).   

    INFO    (output) INTEGER   
            0  - successful exit   
            <0 - if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static logical done;
    static doublereal ctoc;
    static integer i, j;
    static integer itype, k1, k2, k3, k4;
    static doublereal cfrom1;
    static doublereal bignum, smlnum, mul, cto1;
    static doublereal cfromc;



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;

    if (MLFORTRAN(lsame)(type, "G")) {
	itype = 0;
    } else if (MLFORTRAN(lsame)(type, "L")) {
	itype = 1;
    } else if (MLFORTRAN(lsame)(type, "U")) {
	itype = 2;
    } else if (MLFORTRAN(lsame)(type, "H")) {
	itype = 3;
    } else if (MLFORTRAN(lsame)(type, "B")) {
	itype = 4;
    } else if (MLFORTRAN(lsame)(type, "Q")) {
	itype = 5;
    } else if (MLFORTRAN(lsame)(type, "Z")) {
	itype = 6;
    } else {
	itype = -1;
    }

    if (itype == -1) {
	*info = -1;
    } else if (*cfrom == 0.) {
	*info = -4;
    } else if (*m < 0) {
	*info = -6;
    } else if (*n < 0 || (itype == 4 && *n != *m) || (itype == 5 && *n != *m)){
	*info = -7;
    } else if (itype <= 3 && *lda < ML_max(1,*m)) {
	*info = -9;
    } else if (itype >= 4) {
/* Computing MAX */
	i__1 = *m - 1;
	if (*kl < 0 || *kl > ML_max(i__1,0)) {
	    *info = -2;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = *n - 1;
	    if (*ku < 0 || *ku > ML_max(i__1,0) || ((itype == 4 || itype == 5) && 
		    *kl != *ku)) {
		*info = -3;
	    } else if ((itype == 4 && *lda < *kl + 1) ||(itype == 5 && *lda < *
		    ku + 1)||(itype == 6 && *lda < (*kl << 1) + *ku + 1)) {
		*info = -9;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DLASCL", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return 0;
    }

/*     Get machine parameters */

    smlnum = MLFORTRAN(dlamch)("S");
    bignum = 1. / smlnum;

    cfromc = *cfrom;
    ctoc = *cto;

L10:
    cfrom1 = cfromc * smlnum;
    cto1 = ctoc / bignum;
    if (ML_abs(cfrom1) > ML_abs(ctoc) && ctoc != 0.) {
	mul = smlnum;
	done = FALSE_;
	cfromc = cfrom1;
    } else if (ML_abs(cto1) > ML_abs(cfromc)) {
	mul = bignum;
	done = FALSE_;
	ctoc = cto1;
    } else {
	mul = ctoc / cfromc;
	done = TRUE_;
    }

    if (itype == 0) {

/*        Full matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    for (i = 1; i <= *m; ++i) {
		A(i,j) *= mul;
/* L20: */
	    }
/* L30: */
	}

    } else if (itype == 1) {

/*        Lower triangular matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    for (i = j; i <= *m; ++i) {
		A(i,j) *= mul;
/* L40: */
	    }
/* L50: */
	}

    } else if (itype == 2) {

/*        Upper triangular matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    for (i = 1; i <= ML_min(j,*m); ++i) {
		A(i,j) *= mul;
/* L60: */
	    }
/* L70: */
	}

    } else if (itype == 3) {

/*        Upper Hessenberg matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    for (i = 1; i <= ML_min(j+1,*m); ++i) {
		A(i,j) *= mul;
/* L80: */
	    }
/* L90: */
	}

    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

	k3 = *kl + 1;
	k4 = *n + 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    for (i = 1; i <= ML_min(k3,k4-j); ++i) {
		A(i,j) *= mul;
/* L100: */
	    }
/* L110: */
	}

    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

	k1 = *ku + 2;
	k3 = *ku + 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    for (i = ML_max(k1-j,1); i <= k3; ++i) {
		A(i,j) *= mul;
/* L120: */
	    }
/* L130: */
	}

    } else if (itype == 6) {

/*        Band matrix */

	k1 = *kl + *ku + 2;
	k2 = *kl + 1;
	k3 = (*kl << 1) + *ku + 1;
	k4 = *kl + *ku + 1 + *m;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
/* Computing MIN */
	    for (i = ML_max(k1-j,k2); i <= ML_min(k3,k4-j); ++i) {
		A(i,j) *= mul;
/* L140: */
	    }
/* L150: */
	}

    }

    if (! done) {
	goto L10;
    }

    return 0;

/*     End of DLASCL */

} /* dlascl_ */
#endif

#ifndef ML_DLASET_FUNC
/* Subroutine */ int MLFORTRAN(dlaset)(char *uplo, integer *m, integer *n, doublereal *
	alpha, doublereal *beta, doublereal *a, integer *lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASET initializes an m-by-n matrix A to BETA on the diagonal and   
    ALPHA on the offdiagonals.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be set.   
            = 'U':      Upper triangular part is set; the strictly lower 
  
                        triangular part of A is not changed.   
            = 'L':      Lower triangular part is set; the strictly upper 
  
                        triangular part of A is not changed.   
            Otherwise:  All of the matrix A is set.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    ALPHA   (input) DOUBLE PRECISION   
            The constant to which the offdiagonal elements are to be set. 
  

    BETA    (input) DOUBLE PRECISION   
            The constant to which the diagonal elements are to be set.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On exit, the leading m-by-n submatrix of A is set as follows: 
  

            if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,   
            if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,   
            otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,   

            and, for all UPLO, A(i,i) = BETA, 1<=i<=ML_min(m,n).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,M).   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    /* Local variables */
    static integer i, j;



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (MLFORTRAN(lsame)(uplo, "U")) {

/*        Set the strictly upper triangular or trapezoidal part of the
   
          array to ALPHA. */

	for (j = 2; j <= *n; ++j) {
/* Computing MIN */
	    for (i = 1; i <= ML_min(j-1,*m); ++i) {
		A(i,j) = *alpha;
/* L10: */
	    }
/* L20: */
	}

    } else if (MLFORTRAN(lsame)(uplo, "L")) {

/*        Set the strictly lower triangular or trapezoidal part of the
   
          array to ALPHA. */

	for (j = 1; j <= ML_min(*m,*n); ++j) {
	    for (i = j + 1; i <= *m; ++i) {
		A(i,j) = *alpha;
/* L30: */
	    }
/* L40: */
	}

    } else {

/*        Set the leading m-by-n submatrix to ALPHA. */

	for (j = 1; j <= *n; ++j) {
	    for (i = 1; i <= *m; ++i) {
		A(i,j) = *alpha;
/* L50: */
	    }
/* L60: */
	}
    }

/*     Set the first ML_min(M,N) diagonal elements to BETA. */

    for (i = 1; i <= ML_min(*m,*n); ++i) {
	A(i,i) = *beta;
/* L70: */
    }

    return 0;

/*     End of DLASET */

} /* dlaset_ */
#endif


#ifndef ML_DLASSQ_FUNC
/* Subroutine */ int MLFORTRAN(dlassq)(integer *n, doublereal *x, integer *incx, 
	doublereal *scale, doublereal *sumsq)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASSQ  returns the values  scl  and  smsq  such that   

       ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, 
  

    where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is   
    assumed to be non-negative and  scl  returns the value   

       scl = ML_max( scale, ML_abs( x( i ) ) ).   

    scale and sumsq must be supplied in SCALE and SUMSQ and   
    scl and smsq are overwritten on SCALE and SUMSQ respectively.   

    The routine makes only one pass through the vector x.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of elements to be used from the vector X.   

    X       (input) DOUBLE PRECISION   
            The vector for which a scaled sum of squares is computed.   
               x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.   

    INCX    (input) INTEGER   
            The increment between successive values of the vector X.   
            INCX > 0.   

    SCALE   (input/output) DOUBLE PRECISION   
            On entry, the value  scale  in the equation above.   
            On exit, SCALE is overwritten with  scl , the scaling factor 
  
            for the sum of squares.   

    SUMSQ   (input/output) DOUBLE PRECISION   
            On entry, the value  sumsq  in the equation above.   
            On exit, SUMSQ is overwritten with  smsq , the basic sum of   
            squares from which  scl  has been factored out.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    doublereal d__1;
    /* Local variables */
    static doublereal absxi;
    static integer ix;


#define X(I) x[(I)-1]


    if (*n > 0) {
	for (ix = 1; *incx < 0 ? ix >= (*n-1)**incx+1 : ix <= (*n-1)**incx+1; ix += *incx) {
	    if (X(ix) != 0.) {
		absxi = (d__1 = X(ix), ML_abs(d__1));
		if (*scale < absxi) {
/* Computing 2nd power */
		    d__1 = *scale / absxi;
		    *sumsq = *sumsq * (d__1 * d__1) + 1;
		    *scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / *scale;
		    *sumsq += d__1 * d__1;
		}
	    }
/* L10: */
	}
    }
    return 0;

/*     End of DLASSQ */

} /* dlassq_ */
#endif

#ifndef ML_DLANGE_FUNC

doublereal MLFORTRAN(dlange)(char *norm, integer *m, integer *n, doublereal *a, integer 
	*lda, doublereal *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANGE  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    real matrix A.   

    Description   
    ===========   

    DLANGE returns the value   

       DLANGE = ( ML_max(ML_abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
  
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
  
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
  
    squares).  Note that  ML_max(ML_abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANGE as described   
            above.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.  When M = 0,   
            DLANGE is set to zero.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.  When N = 0, 
  
            DLANGE is set to zero.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The m by n matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(M,1).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= M when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;
    /* Local variables */
    static integer i, j;
    static doublereal scale;
    static doublereal value;
    static doublereal sum;



#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (ML_min(*m,*n) == 0) {
	value = 0.;
    } else if (MLFORTRAN(lsame)(norm, "M")) {

/*        Find ML_max(ML_abs(A(i,j))). */

	value = 0.;
	for (j = 1; j <= *n; ++j) {
	    for (i = 1; i <= *m; ++i) {
/* Computing MAX */
		d__2 = value, d__3 = (d__1 = A(i,j), ML_abs(d__1));
		value = ML_max(d__2,d__3);
/* L10: */
	    }
/* L20: */
	}
    } else if (MLFORTRAN(lsame)(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	for (j = 1; j <= *n; ++j) {
	    sum = 0.;
	    for (i = 1; i <= *m; ++i) {
		sum += (d__1 = A(i,j), ML_abs(d__1));
/* L30: */
	    }
	    value = ML_max(value,sum);
/* L40: */
	}
    } else if (MLFORTRAN(lsame)(norm, "I")) {

/*        Find normI(A). */

	for (i = 1; i <= *m; ++i) {
	    WORK(i) = 0.;
/* L50: */
	}
	for (j = 1; j <= *n; ++j) {
	    for (i = 1; i <= *m; ++i) {
		WORK(i) += (d__1 = A(i,j), ML_abs(d__1));
/* L60: */
	    }
/* L70: */
	}
	value = 0.;
	for (i = 1; i <= *m; ++i) {
/* Computing MAX */
	    d__1 = value, d__2 = WORK(i);
	    value = ML_max(d__1,d__2);
/* L80: */
	}
    } else if (MLFORTRAN(lsame)(norm, "F") || MLFORTRAN(lsame)(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	for (j = 1; j <= *n; ++j) {
	    MLFORTRAN(dlassq)(m, &A(1,j), &c__1, &scale, &sum);
/* L90: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANGE */

} /* dlange_ */
#endif

#ifndef ML_DLABAD_FUNC

/* Subroutine */ int MLFORTRAN(dlabad)(doublereal *small, doublereal *large)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLABAD takes as input the values computed by SLAMCH for underflow and 
  
    overflow, and returns the square root of each of these values if the 
  
    log of LARGE is sufficiently large.  This subroutine is intended to   
    identify machines with a large exponent range, such as the Crays, and 
  
    redefine the underflow and overflow limits to be the square roots of 
  
    the values computed by DLAMCH.  This subroutine is needed because   
    DLAMCH does not compensate for poor arithmetic in the upper half of   
    the exponent range, as is found on a Cray.   

    Arguments   
    =========   

    SMALL   (input/output) DOUBLE PRECISION   
            On entry, the underflow threshold as computed by DLAMCH.   
            On exit, if LOG10(LARGE) is sufficiently large, the square   
            root of SMALL, otherwise unchanged.   

    LARGE   (input/output) DOUBLE PRECISION   
            On entry, the overflow threshold as computed by DLAMCH.   
            On exit, if LOG10(LARGE) is sufficiently large, the square   
            root of LARGE, otherwise unchanged.   

    ===================================================================== 
  


       If it looks like we're on a Cray, take the square root of   
       SMALL and LARGE to avoid overflow and underflow problems. */
#ifdef CRAY
    /* Builtin functions */
  double d_lg10(doublereal *);


    if (d_lg10(large) > 2e3) {
	*small = sqrt(*small);
	*large = sqrt(*large);
    }
#endif

    return 0;

/*     End of DLABAD */

} /* dlabad_ */
#endif

#ifndef ML_DORMQR_FUNC

/* Subroutine */ int MLFORTRAN(dormqr)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c, integer *ldc, doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORMQR overwrites the general real M-by-N matrix C with   

                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q**T * C       C * Q**T   

    where Q is a real orthogonal matrix defined as the product of k   
    elementary reflectors   

          Q = H(1) H(2) . . . H(k)   

    as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N   
    if SIDE = 'R'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q**T from the Left;   
            = 'R': apply Q or Q**T from the Right.   

    TRANS   (input) CHARACTER*1   
            = 'N':  No transpose, apply Q;   
            = 'T':  Transpose, apply Q**T.   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines   
            the matrix Q.   
            If SIDE = 'L', M >= K >= 0;   
            if SIDE = 'R', N >= K >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,K)   
            The i-th column must contain the vector which defines the   
            elementary reflector H(i), for i = 1,2,...,k, as returned by 
  
            DGEQRF in the first k columns of its array argument A.   
            A is modified by the routine but restored on exit.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   
            If SIDE = 'L', LDA >= ML_max(1,M);   
            if SIDE = 'R', LDA >= ML_max(1,N).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the M-by-N matrix C.   
            On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. 
  

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= ML_max(1,M).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If SIDE = 'L', LWORK >= ML_max(1,N);   
            if SIDE = 'R', LWORK >= ML_max(1,M).   
            For optimum performance LWORK >= N*NB if SIDE = 'L', and   
            LWORK >= M*NB if SIDE = 'R', where NB is the optimal   
            blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__2 = 2;
    static integer c__65 = 65;
    
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2], i__4, i__5;
    char ch__1[2];
    /* Local variables */
    static logical left;
    static integer i;
    static doublereal t[4160]	/* was [65][64] */;
    static integer nbmin, iinfo, i1, i2, i3;
    static integer ib, ic, jc, nb, mi, ni;
    static integer nq, nw;
    static logical notran;
    static integer ldwork, iws;



#undef T
#define T(I) t[(I)]
#define WAS(I) was[(I)]
#undef TAU
#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    left = MLFORTRAN(lsame)(side, "L");
    notran = MLFORTRAN(lsame)(trans, "N");

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }
    if (! left && ! MLFORTRAN(lsame)(side, "R")) {
	*info = -1;
    } else if (! notran && ! MLFORTRAN(lsame)(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < ML_max(1,nq)) {
	*info = -7;
    } else if (*ldc < ML_max(1,*m)) {
	*info = -10;
    } else if (*lwork < ML_max(1,nw)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DORMQR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	WORK(1) = 1.;
	return 0;
    }

/*     Determine the block size.  NB may be at most NBMAX, where NBMAX   
       is used to define the local array T.   

   Computing MIN   
   Writing concatenation */
    i__3[0] = 1, a__1[0] = side;
    i__3[1] = 1, a__1[1] = trans;
    ml_s_cat(ch__1, a__1, i__3, &c__2, 2L);
    i__1 = 64, i__2 = MLFORTRAN(ilaenv)(&c__1, "DORMQR", ch__1, m, n, k, &c_n1, 6L, 2L);
    nb = ML_min(i__1,i__2);
    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
/* Computing MAX   
   Writing concatenation */
	    i__3[0] = 1, a__1[0] = side;
	    i__3[1] = 1, a__1[1] = trans;
	    ml_s_cat(ch__1, a__1, i__3, &c__2, 2L);
	    i__1 = 2, i__2 = MLFORTRAN(ilaenv)(&c__2, "DORMQR", ch__1, m, n, k, &c_n1, 
		    6L, 2L);
	    nbmin = ML_max(i__1,i__2);
	}
    } else {
	iws = nw;
    }

    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

	MLFORTRAN(dorm2r)(side, trans, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1)
		, ldc, &WORK(1), &iinfo);
    } else {

/*        Use blocked code */

	if ((left && ! notran) || ( ! left && notran) ) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	i__1 = i2;
	i__2 = i3;
	for (i = i1; i3 < 0 ? i >= i2 : i <= i2; i += i3) {
/* Computing MIN */
	    i__4 = nb, i__5 = *k - i + 1;
	    ib = ML_min(i__4,i__5);

/*           Form the triangular factor of the block reflector   
             H = H(i) H(i+1) . . . H(i+ib-1) */

	    i__4 = nq - i + 1;
	    MLFORTRAN(dlarft)("Forward", "Columnwise", &i__4, &ib, &A(i,i), 
		    lda, &TAU(i), t, &c__65);
	    if (left) {

/*              H or H' is applied to C(i:m,1:n) */

		mi = *m - i + 1;
		ic = i;
	    } else {

/*              H or H' is applied to C(1:m,i:n) */

		ni = *n - i + 1;
		jc = i;
	    }

/*           Apply H or H' */

	    MLFORTRAN(dlarfb)(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &A(i,i), lda, t, &c__65, &C(ic,jc), ldc, 
		    &WORK(1), &ldwork);
/* L10: */
	}
    }
    WORK(1) = (doublereal) iws;
    return 0;

/*     End of DORMQR */

} /* dormqr_ */
#endif


#ifndef ML_DORMLQ_FUNC

/* Subroutine */ int MLFORTRAN(dormlq)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c, integer *ldc, doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORMLQ overwrites the general real M-by-N matrix C with   

                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q**T * C       C * Q**T   

    where Q is a real orthogonal matrix defined as the product of k   
    elementary reflectors   

          Q = H(k) . . . H(2) H(1)   

    as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N   
    if SIDE = 'R'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q**T from the Left;   
            = 'R': apply Q or Q**T from the Right.   

    TRANS   (input) CHARACTER*1   
            = 'N':  No transpose, apply Q;   
            = 'T':  Transpose, apply Q**T.   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines   
            the matrix Q.   
            If SIDE = 'L', M >= K >= 0;   
            if SIDE = 'R', N >= K >= 0.   

    A       (input) DOUBLE PRECISION array, dimension   
                                 (LDA,M) if SIDE = 'L',   
                                 (LDA,N) if SIDE = 'R'   
            The i-th row must contain the vector which defines the   
            elementary reflector H(i), for i = 1,2,...,k, as returned by 
  
            DGELQF in the first k rows of its array argument A.   
            A is modified by the routine but restored on exit.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= ML_max(1,K).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGELQF.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the M-by-N matrix C.   
            On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. 
  

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= ML_max(1,M).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If SIDE = 'L', LWORK >= ML_max(1,N);   
            if SIDE = 'R', LWORK >= ML_max(1,M).   
            For optimum performance LWORK >= N*NB if SIDE = 'L', and   
            LWORK >= M*NB if SIDE = 'R', where NB is the optimal   
            blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__2 = 2;
    static integer c__65 = 65;
    
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2], i__4, i__5;
    char ch__1[2];
    /* Local variables */
    static logical left;
    static integer i;
    static doublereal t[4160]	/* was [65][64] */;
    static integer nbmin, iinfo, i1, i2, i3;
    static integer ib, ic, jc, nb, mi, ni;
    static integer nq, nw;
    static logical notran;
    static integer ldwork;
    static char transt[1];
    static integer iws;


#undef T

#define T(I) t[(I)]
#define WAS(I) was[(I)]
#undef TAU
#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    left = MLFORTRAN(lsame)(side, "L");
    notran = MLFORTRAN(lsame)(trans, "N");

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }
    if (! left && ! MLFORTRAN(lsame)(side, "R")) {
	*info = -1;
    } else if (! notran && ! MLFORTRAN(lsame)(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < ML_max(1,*k)) {
	*info = -7;
    } else if (*ldc < ML_max(1,*m)) {
	*info = -10;
    } else if (*lwork < ML_max(1,nw)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DORMLQ", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	WORK(1) = 1.;
	return 0;
    }

/*     Determine the block size.  NB may be at most NBMAX, where NBMAX   
       is used to define the local array T.   

   Computing MIN   
   Writing concatenation */
    i__3[0] = 1, a__1[0] = side;
    i__3[1] = 1, a__1[1] = trans;
    ml_s_cat(ch__1, a__1, i__3, &c__2, 2L);
    i__1 = 64, i__2 = MLFORTRAN(ilaenv)(&c__1, "DORMLQ", ch__1, m, n, k, &c_n1, 6L, 2L);
    nb = ML_min(i__1,i__2);
    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
/* Computing MAX   
   Writing concatenation */
	    i__3[0] = 1, a__1[0] = side;
	    i__3[1] = 1, a__1[1] = trans;
	    ml_s_cat(ch__1, a__1, i__3, &c__2, 2L);
	    i__1 = 2, i__2 = MLFORTRAN(ilaenv)(&c__2, "DORMLQ", ch__1, m, n, k, &c_n1, 
		    6L, 2L);
	    nbmin = ML_max(i__1,i__2);
	}
    } else {
	iws = nw;
    }

    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

	MLFORTRAN(dorml2)(side, trans, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1)
		, ldc, &WORK(1), &iinfo);
    } else {

/*        Use blocked code */

	if ((left && notran) || (! left && ! notran)) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	if (notran) {
	    *(unsigned char *)transt = 'T';
	} else {
	    *(unsigned char *)transt = 'N';
	}

	i__1 = i2;
	i__2 = i3;
	for (i = i1; i3 < 0 ? i >= i2 : i <= i2; i += i3) {
/* Computing MIN */
	    i__4 = nb, i__5 = *k - i + 1;
	    ib = ML_min(i__4,i__5);

/*           Form the triangular factor of the block reflector   
             H = H(i) H(i+1) . . . H(i+ib-1) */

	    i__4 = nq - i + 1;
	    MLFORTRAN(dlarft)("Forward", "Rowwise", &i__4, &ib, &A(i,i), lda,
		     &TAU(i), t, &c__65);
	    if (left) {

/*              H or H' is applied to C(i:m,1:n) */

		mi = *m - i + 1;
		ic = i;
	    } else {

/*              H or H' is applied to C(1:m,i:n) */

		ni = *n - i + 1;
		jc = i;
	    }

/*           Apply H or H' */

	    MLFORTRAN(dlarfb)(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &A(i,i), lda, t, &c__65, &C(ic,jc), ldc, &
		    WORK(1), &ldwork);
/* L10: */
	}
    }
    WORK(1) = (doublereal) iws;
    return 0;

/*     End of DORMLQ */

} /* dormlq_ */
#endif


#ifndef ML_DORM2R_FUNC


/* Subroutine */ int MLFORTRAN(dorm2r)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c, integer *ldc, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORM2R overwrites the general real m by n matrix C with   

          Q * C  if SIDE = 'L' and TRANS = 'N', or   

          Q'* C  if SIDE = 'L' and TRANS = 'T', or   

          C * Q  if SIDE = 'R' and TRANS = 'N', or   

          C * Q' if SIDE = 'R' and TRANS = 'T',   

    where Q is a real orthogonal matrix defined as the product of k   
    elementary reflectors   

          Q = H(1) H(2) . . . H(k)   

    as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n   
    if SIDE = 'R'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q' from the Left   
            = 'R': apply Q or Q' from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply Q  (No transpose)   
            = 'T': apply Q' (Transpose)   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines   
            the matrix Q.   
            If SIDE = 'L', M >= K >= 0;   
            if SIDE = 'R', N >= K >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,K)   
            The i-th column must contain the vector which defines the   
            elementary reflector H(i), for i = 1,2,...,k, as returned by 
  
            DGEQRF in the first k columns of its array argument A.   
            A is modified by the routine but restored on exit.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   
            If SIDE = 'L', LDA >= ML_max(1,M);   
            if SIDE = 'R', LDA >= ML_max(1,N).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= ML_max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension   
                                     (N) if SIDE = 'L',   
                                     (M) if SIDE = 'R'   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static logical left;
    static integer i;
    static integer i1, i2, i3, ic, jc, mi, ni, nq;
    static logical notran;
    static doublereal aii;


#undef WORK

#undef TAU
#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    left = MLFORTRAN(lsame)(side, "L");
    notran = MLFORTRAN(lsame)(trans, "N");

/*     NQ is the order of Q */

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (! left && ! MLFORTRAN(lsame)(side, "R")) {
	*info = -1;
    } else if (! notran && ! MLFORTRAN(lsame)(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < ML_max(1,nq)) {
	*info = -7;
    } else if (*ldc < ML_max(1,*m)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DORM2R", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if ((left && ! notran) || (! left && notran)) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
	jc = 1;
    } else {
	mi = *m;
	ic = 1;
    }

    i__1 = i2;
    for (i = i1; i3 < 0 ? i >= i2 : i <= i2; i += i3) {
	if (left) {

/*           H(i) is applied to C(i:m,1:n) */

	    mi = *m - i + 1;
	    ic = i;
	} else {

/*           H(i) is applied to C(1:m,i:n) */

	    ni = *n - i + 1;
	    jc = i;
	}

/*        Apply H(i) */

	aii = A(i,i);
	A(i,i) = 1.;
	MLFORTRAN(dlarf)(side, &mi, &ni, &A(i,i), &c__1, &TAU(i), &C(ic,jc), ldc, &WORK(1));
	A(i,i) = aii;
/* L10: */
    }
    return 0;

/*     End of DORM2R */

} /* dorm2r_ */
#endif

#ifndef ML_DORML2_FUNC

/* Subroutine */ int MLFORTRAN(dorml2)(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c, integer *ldc, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORML2 overwrites the general real m by n matrix C with   

          Q * C  if SIDE = 'L' and TRANS = 'N', or   

          Q'* C  if SIDE = 'L' and TRANS = 'T', or   

          C * Q  if SIDE = 'R' and TRANS = 'N', or   

          C * Q' if SIDE = 'R' and TRANS = 'T',   

    where Q is a real orthogonal matrix defined as the product of k   
    elementary reflectors   

          Q = H(k) . . . H(2) H(1)   

    as returned by DGELQF. Q is of order m if SIDE = 'L' and of order n   
    if SIDE = 'R'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q' from the Left   
            = 'R': apply Q or Q' from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply Q  (No transpose)   
            = 'T': apply Q' (Transpose)   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines   
            the matrix Q.   
            If SIDE = 'L', M >= K >= 0;   
            if SIDE = 'R', N >= K >= 0.   

    A       (input) DOUBLE PRECISION array, dimension   
                                 (LDA,M) if SIDE = 'L',   
                                 (LDA,N) if SIDE = 'R'   
            The i-th row must contain the vector which defines the   
            elementary reflector H(i), for i = 1,2,...,k, as returned by 
  
            DGELQF in the first k rows of its array argument A.   
            A is modified by the routine but restored on exit.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= ML_max(1,K).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGELQF.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= ML_max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension   
                                     (N) if SIDE = 'L',   
                                     (M) if SIDE = 'R'   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static logical left;
    static integer i;

    static integer i1, i2, i3, ic, jc, mi, ni, nq;
    static logical notran;
    static doublereal aii;


#undef TAU
#define TAU(I) tau[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    left = MLFORTRAN(lsame)(side, "L");
    notran = MLFORTRAN(lsame)(trans, "N");

/*     NQ is the order of Q */

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (! left && ! MLFORTRAN(lsame)(side, "R")) {
	*info = -1;
    } else if (! notran && ! MLFORTRAN(lsame)(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < ML_max(1,*k)) {
	*info = -7;
    } else if (*ldc < ML_max(1,*m)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DORML2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if ((left && notran) || (! left && ! notran)) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
	jc = 1;
    } else {
	mi = *m;
	ic = 1;
    }

    i__1 = i2;
    for (i = i1; i3 < 0 ? i >= i2 : i <= i2; i += i3) {
	if (left) {

/*           H(i) is applied to C(i:m,1:n) */

	    mi = *m - i + 1;
	    ic = i;
	} else {

/*           H(i) is applied to C(1:m,i:n) */

	    ni = *n - i + 1;
	    jc = i;
	}

/*        Apply H(i) */

	aii = A(i,i);
	A(i,i) = 1.;
	MLFORTRAN(dlarf)(side, &mi, &ni, &A(i,i), lda, &TAU(i), &C(ic,jc), ldc, &WORK(1));
	A(i,i) = aii;
/* L10: */
    }
    return 0;

/*     End of DORML2 */

} /* dorml2_ */
#endif



/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifndef ML_DSYRK_FUNC
/* Subroutine */ int MLFORTRAN(dsyrk)(char *uplo, char *trans, integer *n, integer *k, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *beta, 
	doublereal *c, integer *ldc, int dummy1 , int dummy2)
{


    /* System generated locals */


    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i, j, l;
    static integer nrowa;
    static logical upper;


/*  Purpose   
    =======   

    DSYRK  performs one of the symmetric rank k operations   

       C := alpha*A*A' + beta*C,   

    or   

       C := alpha*A'*A + beta*C,   

    where  alpha and beta  are scalars, C is an  n by n  symmetric matrix 
  
    and  A  is an  n by k  matrix in the first case and a  k by n  matrix 
  
    in the second case.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On  entry,   UPLO  specifies  whether  the  upper  or  lower 
  
             triangular  part  of the  array  C  is to be  referenced  as 
  
             follows:   

                UPLO = 'U' or 'u'   Only the  upper triangular part of  C 
  
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the  lower triangular part of  C 
  
                                    is to be referenced.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry,  TRANS  specifies the operation to be performed as 
  
             follows:   

                TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.   

                TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.   

                TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry,  N specifies the order of the matrix C.  N must be 
  
             at least zero.   
             Unchanged on exit.   

    K      - INTEGER.   
             On entry with  TRANS = 'N' or 'n',  K  specifies  the number 
  
             of  columns   of  the   matrix   A,   and  on   entry   with 
  
             TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number 
  
             of rows of the matrix  A.  K must be at least zero.   
             Unchanged on exit.   

    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is 
  
             k  when  TRANS = 'N' or 'n',  and is  n  otherwise.   
             Before entry with  TRANS = 'N' or 'n',  the  leading  n by k 
  
             part of the array  A  must contain the matrix  A,  otherwise 
  
             the leading  k by n  part of the array  A  must contain  the 
  
             matrix A.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' 
  
             then  LDA must be at least  ML_max( 1, n ), otherwise  LDA must 
  
             be at least  ML_max( 1, k ).   
             Unchanged on exit.   

    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta.   
             Unchanged on exit.   

    C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry  with  UPLO = 'U' or 'u',  the leading  n by n 
  
             upper triangular part of the array C must contain the upper 
  
             triangular part  of the  symmetric matrix  and the strictly 
  
             lower triangular part of C is not referenced.  On exit, the 
  
             upper triangular part of the array  C is overwritten by the 
  
             upper triangular part of the updated matrix.   
             Before entry  with  UPLO = 'L' or 'l',  the leading  n by n 
  
             lower triangular part of the array C must contain the lower 
  
             triangular part  of the  symmetric matrix  and the strictly 
  
             upper triangular part of C is not referenced.  On exit, the 
  
             lower triangular part of the array  C is overwritten by the 
  
             lower triangular part of the updated matrix.   

    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared 
  
             in  the  calling  (sub)  program.   LDC  must  be  at  least 
  
             ML_max( 1, n ).   
             Unchanged on exit.   


    Level 3 Blas routine.   

    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    if (MLFORTRAN(lsame)(trans, "N")) {
	nrowa = *n;
    } else {
	nrowa = *k;
    }
    upper = MLFORTRAN(lsame)(uplo, "U");

    info = 0;
    if (! upper && ! MLFORTRAN(lsame)(uplo, "L")) {
	info = 1;
    } else if (! MLFORTRAN(lsame)(trans, "N") && ! MLFORTRAN(lsame)(trans, "T") &&
	     ! MLFORTRAN(lsame)(trans, "C")) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*k < 0) {
	info = 4;
    } else if (*lda < ML_max(1,nrowa)) {
	info = 7;
    } else if (*ldc < ML_max(1,*n)) {
	info = 10;
    }
    if (info != 0) {
	MLFORTRAN(xerbla)("DSYRK ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0 || ((*alpha == 0. || *k == 0) && *beta == 1.)) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (*alpha == 0.) {
	if (upper) {
	    if (*beta == 0.) {
		for (j = 1; j <= *n; ++j) {
		    for (i = 1; i <= j; ++i) {
			C(i,j) = 0.;
/* L10: */
		    }
/* L20: */
		}
	    } else {
		for (j = 1; j <= *n; ++j) {
		    for (i = 1; i <= j; ++i) {
			C(i,j) = *beta * C(i,j);
/* L30: */
		    }
/* L40: */
		}
	    }
	} else {
	    if (*beta == 0.) {
		for (j = 1; j <= *n; ++j) {
		    for (i = j; i <= *n; ++i) {
			C(i,j) = 0.;
/* L50: */
		    }
/* L60: */
		}
	    } else {
		for (j = 1; j <= *n; ++j) {
		    for (i = j; i <= *n; ++i) {
			C(i,j) = *beta * C(i,j);
/* L70: */
		    }
/* L80: */
		}
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (MLFORTRAN(lsame)(trans, "N")) {

/*        Form  C := alpha*A*A' + beta*C. */

	if (upper) {
	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    for (i = 1; i <= j; ++i) {
			C(i,j) = 0.;
/* L90: */
		    }
		} else if (*beta != 1.) {
		    for (i = 1; i <= j; ++i) {
			C(i,j) = *beta * C(i,j);
/* L100: */
		    }
		}
		for (l = 1; l <= *k; ++l) {
		    if (A(j,l) != 0.) {
			temp = *alpha * A(j,l);
			for (i = 1; i <= j; ++i) {
			    C(i,j) += temp * A(i,l);
/* L110: */
			}
		    }
/* L120: */
		}
/* L130: */
	    }
	} else {
	    for (j = 1; j <= *n; ++j) {
		if (*beta == 0.) {
		    for (i = j; i <= *n; ++i) {
			C(i,j) = 0.;
/* L140: */
		    }
		} else if (*beta != 1.) {
		    for (i = j; i <= *n; ++i) {
			C(i,j) = *beta * C(i,j);
/* L150: */
		    }
		}
		for (l = 1; l <= *k; ++l) {
		    if (A(j,l) != 0.) {
			temp = *alpha * A(j,l);
			for (i = j; i <= *n; ++i) {
			    C(i,j) += temp * A(i,l);
/* L160: */
			}
		    }
/* L170: */
		}
/* L180: */
	    }
	}
    } else {

/*        Form  C := alpha*A'*A + beta*C. */

	if (upper) {
	    for (j = 1; j <= *n; ++j) {
		for (i = 1; i <= j; ++i) {
		    temp = 0.;
		    for (l = 1; l <= *k; ++l) {
			temp += A(l,i) * A(l,j);
/* L190: */
		    }
		    if (*beta == 0.) {
			C(i,j) = *alpha * temp;
		    } else {
			C(i,j) = *alpha * temp + *beta * C(i,j);
		    }
/* L200: */
		}
/* L210: */
	    }
	} else {
	    for (j = 1; j <= *n; ++j) {
		for (i = j; i <= *n; ++i) {
		    temp = 0.;
		    for (l = 1; l <= *k; ++l) {
			temp += A(l,i) * A(l,j);
/* L220: */
		    }
		    if (*beta == 0.) {
			C(i,j) = *alpha * temp;
		    } else {
			C(i,j) = *alpha * temp + *beta * C(i,j);
		    }
/* L230: */
		}
/* L240: */
	    }
	}
    }

    return 0;

/*     End of DSYRK . */

} /* dsyrk_ */
#endif

#ifndef ML_DLAIC1_FUNC

/* Subroutine */ int MLFORTRAN(dlaic1)(integer *job, integer *j, doublereal *x, 
	doublereal *sest, doublereal *w, doublereal *gamma, doublereal *
	sestpr, doublereal *s, doublereal *c)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAIC1 applies one step of incremental condition estimation in   
    its simplest version:   

    Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j 
  
    lower triangular matrix L, such that   
             twonorm(L*x) = sest   
    Then DLAIC1 computes sestpr, s, c such that   
    the vector   
                    [ s*x ]   
             xhat = [  c  ]   
    is an approximate singular vector of   
                    [ L     0  ]   
             Lhat = [ w' gamma ]   
    in the sense that   
             twonorm(Lhat*xhat) = sestpr.   

    Depending on JOB, an estimate for the largest or smallest singular   
    value is computed.   

    Note that [s c]' and sestpr**2 is an eigenpair of the system   

        diag(sest*sest, 0) + [alpha  gamma] * [ alpha ]   
                                              [ gamma ]   

    where  alpha =  x'*w.   

    Arguments   
    =========   

    JOB     (input) INTEGER   
            = 1: an estimate for the largest singular value is computed. 
  
            = 2: an estimate for the smallest singular value is computed. 
  

    J       (input) INTEGER   
            Length of X and W   

    X       (input) DOUBLE PRECISION array, dimension (J)   
            The j-vector x.   

    SEST    (input) DOUBLE PRECISION   
            Estimated singular value of j by j matrix L   

    W       (input) DOUBLE PRECISION array, dimension (J)   
            The j-vector w.   

    GAMMA   (input) DOUBLE PRECISION   
            The diagonal element gamma.   

    SEDTPR  (output) DOUBLE PRECISION   
            Estimated singular value of (j+1) by (j+1) matrix Lhat.   

    S       (output) DOUBLE PRECISION   
            Sine needed in forming xhat.   

    C       (output) DOUBLE PRECISION   
            Cosine needed in forming xhat.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b5 = 1.;
    
    /* System generated locals */
    doublereal d__1, d__2, d__3, d__4;
    /* Builtin functions */
    /* Local variables */
    static doublereal sine, test, zeta1, zeta2, b, t, alpha, norma, s1, s2;
    static doublereal absgam, absalp, cosine, absest, eps, tmp;



#define W(I) w[(I)-1]
#define X(I) x[(I)-1]


    eps = MLFORTRAN(dlamch)("Epsilon");
    alpha = MLFORTRAN(ddot)(j, &X(1), &c__1, &W(1), &c__1);

    absalp = abs(alpha);
    absgam = abs(*gamma);
    absest = abs(*sest);

    if (*job == 1) {

/*        Estimating largest singular value   

          special cases */

	if (*sest == 0.) {
	    s1 = ML_max(absgam,absalp);
	    if (s1 == 0.) {
		*s = 0.;
		*c = 1.;
		*sestpr = 0.;
	    } else {
		*s = alpha / s1;
		*c = *gamma / s1;
		tmp = sqrt(*s * *s + *c * *c);
		*s /= tmp;
		*c /= tmp;
		*sestpr = s1 * tmp;
	    }
	    return 0;
	} else if (absgam <= eps * absest) {
	    *s = 1.;
	    *c = 0.;
	    tmp = ML_max(absest,absalp);
	    s1 = absest / tmp;
	    s2 = absalp / tmp;
	    *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
	    return 0;
	} else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
		*s = 1.;
		*c = 0.;
		*sestpr = s2;
	    } else {
		*s = 0.;
		*c = 1.;
		*sestpr = s1;
	    }
	    return 0;
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
		tmp = s1 / s2;
		*s = sqrt(tmp * tmp + 1.);
		*sestpr = s2 * *s;
		*c = *gamma / s2 / *s;
		*s = ml_d_sign(&c_b5, &alpha) / *s;
	    } else {
		tmp = s2 / s1;
		*c = sqrt(tmp * tmp + 1.);
		*sestpr = s1 * *c;
		*s = alpha / s1 / *c;
		*c = ml_d_sign(&c_b5, gamma) / *c;
	    }
	    return 0;
	} else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *gamma / absest;

	    b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
	    *c = zeta1 * zeta1;
	    if (b > 0.) {
		t = *c / (b + sqrt(b * b + *c));
	    } else {
		t = sqrt(b * b + *c) - b;
	    }

	    sine = -zeta1 / t;
	    cosine = -zeta2 / (t + 1.);
	    tmp = sqrt(sine * sine + cosine * cosine);
	    *s = sine / tmp;
	    *c = cosine / tmp;
	    *sestpr = sqrt(t + 1.) * absest;
	    return 0;
	}

    } else if (*job == 2) {

/*        Estimating smallest singular value   

          special cases */

	if (*sest == 0.) {
	    *sestpr = 0.;
	    if (ML_max(absgam,absalp) == 0.) {
		sine = 1.;
		cosine = 0.;
	    } else {
		sine = -(*gamma);
		cosine = alpha;
	    }
/* Computing MAX */
	    d__1 = abs(sine), d__2 = abs(cosine);
	    s1 = ML_max(d__1,d__2);
	    *s = sine / s1;
	    *c = cosine / s1;
	    tmp = sqrt(*s * *s + *c * *c);
	    *s /= tmp;
	    *c /= tmp;
	    return 0;
	} else if (absgam <= eps * absest) {
	    *s = 0.;
	    *c = 1.;
	    *sestpr = absgam;
	    return 0;
	} else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
		*s = 0.;
		*c = 1.;
		*sestpr = s1;
	    } else {
		*s = 1.;
		*c = 0.;
		*sestpr = s2;
	    }
	    return 0;
	} else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
		tmp = s1 / s2;
		*c = sqrt(tmp * tmp + 1.);
		*sestpr = absest * (tmp / *c);
		*s = -(*gamma / s2) / *c;
		*c = ml_d_sign(&c_b5, &alpha) / *c;
	    } else {
		tmp = s2 / s1;
		*s = sqrt(tmp * tmp + 1.);
		*sestpr = absest / *s;
		*c = alpha / s1 / *s;
		*s = -ml_d_sign(&c_b5, gamma) / *s;
	    }
	    return 0;
	} else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *gamma / absest;

/* Computing MAX */
	    d__3 = zeta1 * zeta1 + 1. + (d__1 = zeta1 * zeta2, abs(d__1)), 
		    d__4 = (d__2 = zeta1 * zeta2, abs(d__2)) + zeta2 * zeta2;
	    norma = ML_max(d__3,d__4);

/*           See if root is closer to zero or to ONE */

	    test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
	    if (test >= 0.) {

/*              root is close to zero, compute directly */

		b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
		*c = zeta2 * zeta2;
		t = *c / (b + sqrt((d__1 = b * b - *c, abs(d__1))));
		sine = zeta1 / (1. - t);
		cosine = -zeta2 / t;
		*sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
	    } else {

/*              root is closer to ONE, shift by that amount */

		b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
		*c = zeta1 * zeta1;
		if (b >= 0.) {
		    t = -(*c) / (b + sqrt(b * b + *c));
		} else {
		    t = b - sqrt(b * b + *c);
		}
		sine = -zeta1 / t;
		cosine = -zeta2 / (t + 1.);
		*sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
	    }
	    tmp = sqrt(sine * sine + cosine * cosine);
	    *s = sine / tmp;
	    *c = cosine / tmp;
	    return 0;

	}
    }
    return 0;

/*     End of DLAIC1 */

} /* dlaic1_ */
#endif

#ifndef ML_DGETRI_FUNC
/* Subroutine */ int MLFORTRAN(dgetri)(integer *n, doublereal *a, integer *lda, integer 
	*ipiv, doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGETRI computes the inverse of a matrix using the LU factorization   
    computed by DGETRF.   

    This method inverts U and then computes inv(A) by solving the system 
  
    inv(A)*L = inv(U) for inv(A).   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the factors L and U from the factorization   
            A = P*L*U as computed by DGETRF.   
            On exit, if INFO = 0, the inverse of the original matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,N).   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices from DGETRF; for 1<=i<=N, row i of the   
            matrix was interchanged with row IPIV(i).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO=0, then WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= ML_max(1,N).   
            For optimal performance LWORK >= N*NB, where NB is   
            the optimal blocksize returned by ILAENV.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is   
                  singular and its inverse could not be computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__2 = 2;
    static doublereal c_b20 = -1.;
    static doublereal c_b22 = 1.;
    
    /* System generated locals */
    integer i__1, i__2, i__3;
    /* Local variables */
    static integer i, j;
    static integer nbmin;
    static integer jb, nb, jj, jp, nn;
    static integer ldwork;
    static integer iws;



#define IPIV(I) ipiv[(I)-1]
#undef WORK
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    WORK(1) = (doublereal) ML_max(*n,1);
    if (*n < 0) {
	*info = -1;
    } else if (*lda < ML_max(1,*n)) {
	*info = -3;
    } else if (*lwork < ML_max(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DGETRI", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,   
       and the inverse is not computed. */

    MLFORTRAN(dtrtri)("Upper", "Non-unit", n, &A(1,1), lda, info,
            strlen("Upper"), strlen("Non-unit"));
    if (*info > 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = MLFORTRAN(ilaenv)(&c__1, "DGETRI", " ", n, &c_n1, &c_n1, &c_n1, 6L, 1L);
    nbmin = 2;
    ldwork = *n;
    if (nb > 1 && nb < *n) {
/* Computing MAX */
	i__1 = ldwork * nb;
	iws = ML_max(i__1,1);
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
/* Computing MAX */
	    i__1 = 2, i__2 = MLFORTRAN(ilaenv)(&c__2, "DGETRI", " ", n, &c_n1, &c_n1, &
		    c_n1, 6L, 1L);
	    nbmin = ML_max(i__1,i__2);
	}
    } else {
	iws = *n;
    }

/*     Solve the equation inv(A)*L = inv(U) for inv(A). */

    if (nb < nbmin || nb >= *n) {

/*        Use unblocked code. */

	for (j = *n; j >= 1; --j) {

/*           Copy current column of L to WORK and replace with zer
os. */

	    i__1 = *n;
	    for (i = j + 1; i <= *n; ++i) {
		WORK(i) = A(i,j);
		A(i,j) = 0.;
/* L10: */
	    }

/*           Compute current column of inv(A). */

	    if (j < *n) {
		i__1 = *n - j;
		MLFORTRAN(dgemv)("No transpose", n, &i__1, &c_b20, &A(1,j+1), lda, &WORK(j + 1), &c__1, &c_b22, &A(1,j), &c__1, strlen("No transpose"));
	    }
/* L20: */
	}
    } else {

/*        Use blocked code. */

	nn = (*n - 1) / nb * nb + 1;
	i__1 = -nb;
	for (j = nn; -nb < 0 ? j >= 1 : j <= 1; j += -nb) {
/* Computing MIN */
	    i__2 = nb, i__3 = *n - j + 1;
	    jb = ML_min(i__2,i__3);

/*           Copy current block column of L to WORK and replace wi
th   
             zeros. */

	    i__2 = j + jb - 1;
	    for (jj = j; jj <= j+jb-1; ++jj) {
		i__3 = *n;
		for (i = jj + 1; i <= *n; ++i) {
		    WORK(i + (jj - j) * ldwork) = A(i,jj);
		    A(i,jj) = 0.;
/* L30: */
		}
/* L40: */
	    }

/*           Compute current block column of inv(A). */

	    if (j + jb <= *n) {
		i__2 = *n - j - jb + 1;
		MLFORTRAN(dgemm)("No transpose", "No transpose", n, &jb, &i__2, &c_b20, 
			&A(1,j+jb), lda, &WORK(j + jb), &
			ldwork, &c_b22, &A(1,j), lda,
		       strlen("No transpose"), strlen("No transpose"));
	    }
	    MLFORTRAN(dtrsm)("Right", "Lower", "No transpose", "Unit", n, &jb, &c_b22, &
		    WORK(j), &ldwork, &A(1,j), lda
         ,strlen("Right"),strlen("Lower"),strlen("No transpose"),strlen("Unit")
                  );
/* L50: */
	}
    }

/*     Apply column interchanges. */

    for (j = *n - 1; j >= 1; --j) {
	jp = IPIV(j);
	if (jp != j) {
	    MLFORTRAN(dswap)(n, &A(1,j), &c__1, &A(1,jp), &c__1);
	}
/* L60: */
    }

    WORK(1) = (doublereal) iws;
    return 0;

/*     End of DGETRI */

} /* dgetri_ */
#endif

#ifndef ML_DPOTRF_FUNC

/* Subroutine */ int MLFORTRAN(dpotrf)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info, int dummy1)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPOTRF computes the Cholesky factorization of a real symmetric   
    positive definite matrix A.   

    The factorization has the form   
       A = U**T * U,  if UPLO = 'U', or   
       A = L  * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the block version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U**T*U or A = L*L**T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static doublereal c_b13 = -1.;
    static doublereal c_b14 = 1.;
    
    /* System generated locals */
    integer i__1, i__3, i__4;
    /* Local variables */
    static integer j;
    static logical upper;
    static integer jb, nb;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = MLFORTRAN(lsame)(uplo, "U");
    if (! upper && ! MLFORTRAN(lsame)(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ML_max(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DPOTRF", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Determine the block size for this environment. */

    nb = MLFORTRAN(ilaenv)(&c__1, "DPOTRF", uplo, n, &c_n1, &c_n1, &c_n1, 6L, 1L);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code. */

	MLFORTRAN(dpotf2)(uplo, n, &A(1,1), lda, info, strlen(uplo));
    } else {

/*        Use blocked code. */

	if (upper) {

/*           Compute the Cholesky factorization A = U'*U. */

	    i__1 = *n;
	    for (j = 1; nb < 0 ? j >= *n : j <= *n; j += nb) {

/*              Update and factorize the current diagonal bloc
k and test   
                for non-positive-definiteness.   

   Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = ML_min(i__3,i__4);
		i__3 = j - 1;
		MLFORTRAN(dsyrk)("Upper", "Transpose", &jb, &i__3, &c_b13, &A(1,j), lda, &c_b14, &A(j,j), lda, strlen("Upper"), strlen("Transpose"));
		MLFORTRAN(dpotf2)("Upper", &jb, &A(j,j), lda, info, strlen("Upper"));
		if (*info != 0) {
		    goto L30;
		}
		if (j + jb <= *n) {

/*                 Compute the current block row. */

		    i__3 = *n - j - jb + 1;
		    i__4 = j - 1;
		    MLFORTRAN(dgemm)("Transpose", "No transpose", &jb, &i__3, &i__4, &
			    c_b13, &A(1,j), lda, &A(1,j+jb), lda, &c_b14, &A(j,j+jb), lda
                          ,strlen("Transpose"),strlen("No transpose")
);
		    i__3 = *n - j - jb + 1;
		    MLFORTRAN(dtrsm)("Left", "Upper", "Transpose", "Non-unit", &jb, &
			    i__3, &c_b14, &A(j,j), lda, &A(j,j+jb), lda
 ,strlen("Left"),strlen("Upper"),strlen("Transpose"),strlen("Non-unit")
);
		}
/* L10: */
	    }

	} else {

/*           Compute the Cholesky factorization A = L*L'. */

	    i__1 = nb;
	    for (j = 1; nb < 0 ? j >= *n : j <= *n; j += nb) {

/*              Update and factorize the current diagonal bloc
k and test   
                for non-positive-definiteness.   

   Computing MIN */
		i__3 = nb, i__4 = *n - j + 1;
		jb = ML_min(i__3,i__4);
		i__3 = j - 1;
		MLFORTRAN(dsyrk)("Lower", "No transpose", &jb, &i__3, &c_b13, &A(j,1), lda, &c_b14, &A(j,j), lda, strlen("Lower"), strlen("No transpose"));
		MLFORTRAN(dpotf2)("Lower", &jb, &A(j,j), lda, info, strlen("Lower"));
		if (*info != 0) {
		    goto L30;
		}
		if (j + jb <= *n) {

/*                 Compute the current block column. */

		    i__3 = *n - j - jb + 1;
		    i__4 = j - 1;
		    MLFORTRAN(dgemm)("No transpose", "Transpose", &i__3, &jb, &i__4, &
			    c_b13, &A(j+jb,1), lda, &A(j,1), 
			    lda, &c_b14, &A(j+jb,j), lda
                          ,strlen("No transpose"), strlen("Transpose")   );
		    i__3 = *n - j - jb + 1;
		    MLFORTRAN(dtrsm)("Right", "Lower", "Transpose", "Non-unit", &i__3, &
			    jb, &c_b14, &A(j,j), lda, &A(j+jb,j), lda
 ,strlen("Right"),strlen("Lower"),strlen("Transpose"),strlen("Non-unit"));
		}
/* L20: */
	    }
	}
    }
    goto L40;

L30:
    *info = *info + j - 1;

L40:
    return 0;

/*     End of DPOTRF */

} /* dpotrf_ */

#endif

#ifndef ML_DTRTRI_FUNC

/* Subroutine */ int MLFORTRAN(dtrtri)(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info, int dummy1, int dummy2)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DTRTRI computes the inverse of a real upper or lower triangular   
    matrix A.   

    This is the Level 3 BLAS version of the algorithm.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    DIAG    (input) CHARACTER*1   
            = 'N':  A is non-unit triangular;   
            = 'U':  A is unit triangular.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the triangular matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of the array A contains 
  
            the upper triangular matrix, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of the array A contains 
  
            the lower triangular matrix, and the strictly upper   
            triangular part of A is not referenced.  If DIAG = 'U', the   
            diagonal elements of A are also not referenced and are   
            assumed to be 1.   
            On exit, the (triangular) inverse of the original matrix, in 
  
            the same storage format.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, A(i,i) is exactly zero.  The triangular   
                 matrix is singular and its inverse can not be computed. 
  

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__2 = 2;
    static doublereal c_b18 = 1.;
    static doublereal c_b22 = -1.;
    
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2[2], i__4, i__5;
    char ch__1[2];

    /* Local variables */
    static integer j;
    static logical upper;
    static integer jb, nb, nn;
    static logical nounit;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = MLFORTRAN(lsame)(uplo, "U");
    nounit = MLFORTRAN(lsame)(diag, "N");
    if (! upper && ! MLFORTRAN(lsame)(uplo, "L")) {
	*info = -1;
    } else if (! nounit && ! MLFORTRAN(lsame)(diag, "U")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < ML_max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DTRTRI", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Check for singularity if non-unit. */

    if (nounit) {
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (A(*info,*info) == 0.) {
		return 0;
	    }
/* L10: */
	}
	*info = 0;
    }

/*     Determine the block size for this environment.   

   Writing concatenation */
    i__2[0] = 1, a__1[0] = uplo;
    i__2[1] = 1, a__1[1] = diag;
    ml_s_cat(ch__1, a__1, i__2, &c__2, 2L);
    nb = MLFORTRAN(ilaenv)(&c__1, "DTRTRI", ch__1, n, &c_n1, &c_n1, &c_n1, 6L, 2L);
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */

	MLFORTRAN(dtrti2)(uplo, diag, n, &A(1,1), lda, info,strlen(uplo),strlen(diag));
    } else {

/*        Use blocked code */

	if (upper) {

/*           Compute inverse of upper triangular matrix */

	    i__1 = *n;
	    for (j = 1; nb < 0 ? j >= *n : j <= *n; j += nb) {
/* Computing MIN */
		i__4 = nb, i__5 = *n - j + 1;
		jb = ML_min(i__4,i__5);

/*              Compute rows 1:j-1 of current block column */

		i__4 = j - 1;
		MLFORTRAN(dtrmm)("Left", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b18, &A(1,1), lda, &A(1,j), lda
 ,strlen("Left"),strlen("Upper"),strlen("No transpose"),strlen(diag));

		i__4 = j - 1;
		MLFORTRAN(dtrsm)("Right", "Upper", "No transpose", diag, &i__4, &jb, &
			c_b22, &A(j,j), lda, &A(1,j), 
			lda
 ,strlen("Right"),strlen("Upper"),strlen("No transpose"), strlen(diag));


/*              Compute inverse of current diagonal block */

		MLFORTRAN(dtrti2)("Upper", diag, &jb, &A(j,j), lda, info,strlen("Upper"),strlen(diag));
/* L20: */
	    }
	} else {

/*           Compute inverse of lower triangular matrix */

	    nn = (*n - 1) / nb * nb + 1;
	    for (j = nn; -nb < 0 ? j >= 1 : j <= 1; j += -nb) {
/* Computing MIN */
		i__1 = nb, i__4 = *n - j + 1;
		jb = ML_min(i__1,i__4);
		if (j + jb <= *n) {

/*                 Compute rows j+jb:n of current block co
lumn */

		    i__1 = *n - j - jb + 1;
		    MLFORTRAN(dtrmm)("Left", "Lower", "No transpose", diag, &i__1, &jb, 
			    &c_b18, &A(j+jb,j+jb), lda, &A(j+jb,j), lda
 ,strlen("Left"),strlen("Lower"),strlen("No transpose"),strlen(diag));

		    i__1 = *n - j - jb + 1;
		    MLFORTRAN(dtrsm)("Right", "Lower", "No transpose", diag, &i__1, &jb,
			     &c_b22, &A(j,j), lda, &A(j+jb,j), lda
 ,strlen("Right"),strlen("Lower"),strlen("No transpose"),strlen(diag));
		}

/*              Compute inverse of current diagonal block */

		MLFORTRAN(dtrti2)("Lower", diag, &jb, &A(j,j), lda, info,strlen("Lower"),strlen(diag));
/* L30: */
	    }
	}
    }

    return 0;

/*     End of DTRTRI */

} /* dtrtri_ */
#endif

#ifndef ML_DPOTF2_FUNC

/* Subroutine */ int MLFORTRAN(dpotf2)(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info, int dummy1)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DPOTF2 computes the Cholesky factorization of a real symmetric   
    positive definite matrix A.   

    The factorization has the form   
       A = U' * U ,  if UPLO = 'U', or   
       A = L  * L',  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the unblocked version of the algorithm, calling Level 2 BLAS. 
  

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored.   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            n by n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U'*U  or A = L*L'.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, the leading minor of order k is not   
                 positive definite, and the factorization could not be   
                 completed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b10 = -1.;
    static doublereal c_b12 = 1.;
    
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static logical upper;
    static doublereal ajj;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = MLFORTRAN(lsame)(uplo, "U");
    if (! upper && ! MLFORTRAN(lsame)(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < ML_max(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DPOTF2", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (upper) {

/*        Compute the Cholesky factorization A = U'*U. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness
. */

	    i__2 = j - 1;
	    ajj = A(j,j) - MLFORTRAN(ddot)(&i__2, &A(1,j), &c__1, 
		    &A(1,j), &c__1);
	    if (ajj <= 0.) {
		A(j,j) = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    A(j,j) = ajj;

/*           Compute elements J+1:N of row J. */

	    if (j < *n) {
		i__2 = j - 1;
		i__3 = *n - j;
		MLFORTRAN(dgemv)("Transpose", &i__2, &i__3, &c_b10, &A(1,j+1), lda, &A(1,j), &c__1, &c_b12, &A(j,j+1), lda
 ,strlen("Transpose"));
		i__2 = *n - j;
		d__1 = 1. / ajj;
		MLFORTRAN(dscal)(&i__2, &d__1, &A(j,j+1), lda);
	    }
/* L10: */
	}
    } else {

/*        Compute the Cholesky factorization A = L*L'. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness
. */

	    i__2 = j - 1;
	    ajj = A(j,j) - MLFORTRAN(ddot)(&i__2, &A(j,1), lda, &A(j,1), lda);
	    if (ajj <= 0.) {
		A(j,j) = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    A(j,j) = ajj;

/*           Compute elements J+1:N of column J. */

	    if (j < *n) {
		i__2 = *n - j;
		i__3 = j - 1;
		MLFORTRAN(dgemv)("No transpose", &i__2, &i__3, &c_b10, &A(j+1,1), lda, &A(j,1), lda, &c_b12, &A(j+1,j), &c__1, strlen("No transpose"));
		i__2 = *n - j;
		d__1 = 1. / ajj;
		MLFORTRAN(dscal)(&i__2, &d__1, &A(j+1,j), &c__1);
	    }
/* L20: */
	}
    }
    goto L40;

L30:
    *info = j;

L40:
    return 0;

/*     End of DPOTF2 */

} /* dpotf2_ */
#endif

#ifndef ML_DTRTI2_FUNC

/* Subroutine */ int MLFORTRAN(dtrti2)(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info, int dummy1, int dummy2)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DTRTI2 computes the inverse of a real upper or lower triangular   
    matrix.   

    This is the Level 2 BLAS version of the algorithm.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower triangular. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A is unit triangular.   
            = 'N':  Non-unit triangular   
            = 'U':  Unit triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the triangular matrix A.  If UPLO = 'U', the   
            leading n by n upper triangular part of the array A contains 
  
            the upper triangular matrix, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of the array A contains 
  
            the lower triangular matrix, and the strictly upper   
            triangular part of A is not referenced.  If DIAG = 'U', the   
            diagonal elements of A are also not referenced and are   
            assumed to be 1.   

            On exit, the (triangular) inverse of the original matrix, in 
  
            the same storage format.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= ML_max(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    static integer j;
    static logical upper;
    static logical nounit;
    static doublereal ajj;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = MLFORTRAN(lsame)(uplo, "U");
    nounit = MLFORTRAN(lsame)(diag, "N");
    if (! upper && ! MLFORTRAN(lsame)(uplo, "L")) {
	*info = -1;
    } else if (! nounit && ! MLFORTRAN(lsame)(diag, "U")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < ML_max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	MLFORTRAN(xerbla)("DTRTI2", &i__1);
	return 0;
    }

    if (upper) {

/*        Compute inverse of upper triangular matrix. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (nounit) {
		A(j,j) = 1. / A(j,j);
		ajj = -A(j,j);
	    } else {
		ajj = -1.;
	    }

/*           Compute elements 1:j-1 of j-th column. */

	    i__2 = j - 1;
	    MLFORTRAN(dtrmv)("Upper", "No transpose", diag, &i__2, &A(1,1), lda, &
		    A(1,j), &c__1);
	    i__2 = j - 1;
	    MLFORTRAN(dscal)(&i__2, &ajj, &A(1,j), &c__1);
/* L10: */
	}
    } else {

/*        Compute inverse of lower triangular matrix. */

	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		A(j,j) = 1. / A(j,j);
		ajj = -A(j,j);
	    } else {
		ajj = -1.;
	    }
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

		i__1 = *n - j;
		MLFORTRAN(dtrmv)("Lower", "No transpose", diag, &i__1, &A(j+1,j+1), lda, &A(j+1,j), &c__1);
		i__1 = *n - j;
		MLFORTRAN(dscal)(&i__1, &ajj, &A(j+1,j), &c__1);
	    }
/* L20: */
	}
    }

    return 0;

/*     End of DTRTI2 */

} /* dtrti2_ */

#endif
