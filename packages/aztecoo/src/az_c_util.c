/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

/* az_f_util.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

   http://www.netlib.org/f2c/libf2c.zip
*/

#include "AztecOO_config.h"

#define SSWAP_F77 F77_BLAS_MANGLE(sswap,SSWAP)
#define SDOT_F77 F77_BLAS_MANGLE(sdot,SDOT)
#define SLAMCH_F77 F77_BLAS_MANGLE(slamch,SLAMCH)

#define DDOT_F77 F77_BLAS_MANGLE(ddot,DDOT)
#define DLAMCH_F77 F77_BLAS_MANGLE(dlamch,DLAMCH)

/* From f2c.h  --  Standard Fortran to C header file */

typedef int integer;
typedef unsigned long int uinteger;
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
#ifdef INTEGER_STAR_8	/* Adjust for integer*8. */
typedef long long longint;		/* system-dependent */
typedef unsigned long long ulongint;	/* system-dependent */
#define qbit_clear(a,b)	((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)	((a) |  ((ulongint)1 << (b)))
#endif

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
typedef short ftnlen;
typedef short ftnint;
#else
typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;
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

/*typedef long int Long;*/	/* No longer used; formerly in Namelist */

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

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
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

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b5 = 1.;
static real c_b19 = 1.f;

double az_d_sign(doublereal * a, doublereal * b) {
  double x;
  x = (*a >= 0 ? *a : -*a);
  return (*b >= 0 ? x : -x);
}

/* Subroutine */ int az_dlaic1_c(integer *job, integer *j, doublereal *x, 
	doublereal *sest, doublereal *w, doublereal *gamma, doublereal *
	sestpr, doublereal *s, doublereal *c__)
{
  /* System generated locals */
  doublereal d__1, d__2, d__3, d__4;

  /* Builtin functions */
  double sqrt(doublereal);
  /* Local variables */
  static doublereal b, t, s1, s2, eps, tmp;
  extern doublereal DDOT_F77(integer *, doublereal *, integer *, doublereal *, 
    integer *);
  static doublereal sine, test, zeta1, zeta2, alpha, norma;
  extern doublereal DLAMCH_F77(char *, ftnlen);
  static doublereal absgam, absalp, cosine, absest;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  AZ_DLAIC1 applies one step of incremental condition estimation in */
/*  its simplest version: */

/*  Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/*  lower triangular matrix L, such that */
/*           twonorm(L*x) = sest */
/*  Then AZ_DLAIC1 computes sestpr, s, c such that */
/*  the vector */
/*                  [ s*x ] */
/*           xhat = [  c  ] */
/*  is an approximate singular vector of */
/*                  [ L     0  ] */
/*           Lhat = [ w' gamma ] */
/*  in the sense that */
/*           twonorm(Lhat*xhat) = sestpr. */

/*  Depending on JOB, an estimate for the largest or smallest singular */
/*  value is computed. */

/*  Note that [s c]' and sestpr**2 is an eigenpair of the system */

/*      diag(sest*sest, 0) + [alpha  gamma] * [ alpha ] */
/*                                            [ gamma ] */

/*  where  alpha =  x'*w. */

/*  Arguments */
/*  ========= */

/*  JOB     (input) INTEGER */
/*          = 1: an estimate for the largest singular value is computed. */
/*          = 2: an estimate for the smallest singular value is computed. */

/*  J       (input) INTEGER */
/*          Length of X and W */

/*  X       (input) DOUBLE PRECISION array, dimension (J) */
/*          The j-vector x. */

/*  SEST    (input) DOUBLE PRECISION */
/*          Estimated singular value of j by j matrix L */

/*  W       (input) DOUBLE PRECISION array, dimension (J) */
/*          The j-vector w. */

/*  GAMMA   (input) DOUBLE PRECISION */
/*          The diagonal element gamma. */

/*  SEDTPR  (output) DOUBLE PRECISION */
/*          Estimated singular value of (j+1) by (j+1) matrix Lhat. */

/*  S       (output) DOUBLE PRECISION */
/*          Sine needed in forming xhat. */

/*  C       (output) DOUBLE PRECISION */
/*          Cosine needed in forming xhat. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

  /* Parameter adjustments */
  --w;
  --x;

  /* Function Body */
  eps = DLAMCH_F77("Epsilon", (ftnlen)7);
  alpha = DDOT_F77(j, &x[1], &c__1, &w[1], &c__1);

  absalp = abs(alpha);
  absgam = abs(*gamma);
  absest = abs(*sest);

  if (*job == 1) {

/*        Estimating largest singular value */

/*        special cases */

    if (*sest == 0.) {
	    s1 = max(absgam,absalp);
	    if (s1 == 0.) {
        *s = 0.;
        *c__ = 1.;
        *sestpr = 0.;
	    } else {
        *s = alpha / s1;
        *c__ = *gamma / s1;
        tmp = sqrt(*s * *s + *c__ * *c__);
        *s /= tmp;
        *c__ /= tmp;
        *sestpr = s1 * tmp;
	    }
	    return 0;
    } else if (absgam <= eps * absest) {
	    *s = 1.;
	    *c__ = 0.;
	    tmp = max(absest,absalp);
	    s1 = absest / tmp;
	    s2 = absalp / tmp;
	    *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
	    return 0;
    } else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
        *s = 1.;
        *c__ = 0.;
        *sestpr = s2;
	    } else {
        *s = 0.;
        *c__ = 1.;
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
        *c__ = *gamma / s2 / *s;
        *s = az_d_sign(&c_b5, &alpha) / *s;
	    } else {
        tmp = s2 / s1;
        *c__ = sqrt(tmp * tmp + 1.);
        *sestpr = s1 * *c__;
        *s = alpha / s1 / *c__;
        *c__ = az_d_sign(&c_b5, gamma) / *c__;
	    }
	    return 0;
    } else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *gamma / absest;

	    b = (1. - zeta1 * zeta1 - zeta2 * zeta2) * .5;
	    *c__ = zeta1 * zeta1;
	    if (b > 0.) {
        t = *c__ / (b + sqrt(b * b + *c__));
	    } else {
        t = sqrt(b * b + *c__) - b;
	    }

	    sine = -zeta1 / t;
	    cosine = -zeta2 / (t + 1.);
	    tmp = sqrt(sine * sine + cosine * cosine);
	    *s = sine / tmp;
	    *c__ = cosine / tmp;
	    *sestpr = sqrt(t + 1.) * absest;
	    return 0;
    }

  } else if (*job == 2) {

/*        Estimating smallest singular value */

/*        special cases */

    if (*sest == 0.) {
	    *sestpr = 0.;
	    if (max(absgam,absalp) == 0.) {
        sine = 1.;
        cosine = 0.;
	    } else {
        sine = -(*gamma);
        cosine = alpha;
	    }
/* Computing MAX */
	    d__1 = abs(sine), d__2 = abs(cosine);
	    s1 = max(d__1,d__2);
	    *s = sine / s1;
	    *c__ = cosine / s1;
	    tmp = sqrt(*s * *s + *c__ * *c__);
	    *s /= tmp;
	    *c__ /= tmp;
	    return 0;
    } else if (absgam <= eps * absest) {
	    *s = 0.;
	    *c__ = 1.;
	    *sestpr = absgam;
	    return 0;
    } else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
        *s = 0.;
        *c__ = 1.;
        *sestpr = s1;
	    } else {
        *s = 1.;
        *c__ = 0.;
        *sestpr = s2;
	    }
	    return 0;
    } else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
        tmp = s1 / s2;
        *c__ = sqrt(tmp * tmp + 1.);
        *sestpr = absest * (tmp / *c__);
        *s = -(*gamma / s2) / *c__;
        *c__ = az_d_sign(&c_b5, &alpha) / *c__;
	    } else {
        tmp = s2 / s1;
        *s = sqrt(tmp * tmp + 1.);
        *sestpr = absest / *s;
        *c__ = alpha / s1 / *s;
        *s = -az_d_sign(&c_b5, gamma) / *s;
	    }
	    return 0;
    } else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *gamma / absest;

/* Computing MAX */
	    d__3 = zeta1 * zeta1 + 1. + (d__1 = zeta1 * zeta2, abs(d__1)), 
		    d__4 = (d__2 = zeta1 * zeta2, abs(d__2)) + zeta2 * zeta2;
	    norma = max(d__3,d__4);

/*           See if root is closer to zero or to ONE */

	    test = (zeta1 - zeta2) * 2. * (zeta1 + zeta2) + 1.;
	    if (test >= 0.) {

/*              root is close to zero, compute directly */

        b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.) * .5;
        *c__ = zeta2 * zeta2;
        t = *c__ / (b + sqrt((d__1 = b * b - *c__, abs(d__1))));
        sine = zeta1 / (1. - t);
        cosine = -zeta2 / t;
        *sestpr = sqrt(t + eps * 4. * eps * norma) * absest;
	    } else {

/*              root is closer to ONE, shift by that amount */

        b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.) * .5;
        *c__ = zeta1 * zeta1;
        if (b >= 0.) {
          t = -(*c__) / (b + sqrt(b * b + *c__));
        } else {
          t = b - sqrt(b * b + *c__);
        }
        sine = -zeta1 / t;
        cosine = -zeta2 / (t + 1.);
        *sestpr = sqrt(t + 1. + eps * 4. * eps * norma) * absest;
	    }
	    tmp = sqrt(sine * sine + cosine * cosine);
	    *s = sine / tmp;
	    *c__ = cosine / tmp;
	    return 0;

    }
  }
  return 0;

/*     End of AZ_DLAIC1 */

} /* az_dlaic1 */

/* Subroutine */ int az_dlaswp_c(integer *n, doublereal *a, integer *lda, 
	integer *k1, integer *k2, integer *ipiv, integer *incx)
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

  /* Local variables */
  static integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
  static doublereal temp;


/*  -- LAPACK auxiliary routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     June 30, 1999 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  AZ_DLASWP performs a series of row interchanges on the matrix A. */
/*  One row interchange is initiated for each of rows K1 through K2 of A. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A. */

/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*          On entry, the matrix of column dimension N to which the row */
/*          interchanges will be applied. */
/*          On exit, the permuted matrix. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. */

/*  K1      (input) INTEGER */
/*          The first element of IPIV for which a row interchange will */
/*          be done. */

/*  K2      (input) INTEGER */
/*          The last element of IPIV for which a row interchange will */
/*          be done. */

/*  IPIV    (input) INTEGER array, dimension (M*abs(INCX)) */
/*          The vector of pivot indices.  Only the elements in positions */
/*          K1 through K2 of IPIV are accessed. */
/*          IPIV(K) = L implies rows K and L are to be interchanged. */

/*  INCX    (input) INTEGER */
/*          The increment between successive values of IPIV.  If IPIV */
/*          is negative, the pivots are applied in reverse order. */

/*  Further Details */
/*  =============== */

/*  Modified by */
/*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Interchange row I with row IPIV(I) for each of rows K1 through K2. */

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --ipiv;

  /* Function Body */
  if (*incx > 0) {
    ix0 = *k1;
    i1 = *k1;
    i2 = *k2;
    inc = 1;
  } else if (*incx < 0) {
    ix0 = (1 - *k2) * *incx + 1;
    i1 = *k2;
    i2 = *k1;
    inc = -1;
  } else {
    return 0;
  }

  n32 = *n / 32 << 5;
  if (n32 != 0) {
    i__1 = n32;
    for (j = 1; j <= i__1; j += 32) {
	    ix = ix0;
	    i__2 = i2;
	    i__3 = inc;
	    for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
      {
        ip = ipiv[ix];
        if (ip != i__) {
          i__4 = j + 31;
          for (k = j; k <= i__4; ++k) {
            temp = a[i__ + k * a_dim1];
            a[i__ + k * a_dim1] = a[ip + k * a_dim1];
            a[ip + k * a_dim1] = temp;
/* L10: */
          }
        }
        ix += *incx;
/* L20: */
	    }
/* L30: */
    }
  }
  if (n32 != *n) {
    ++n32;
    ix = ix0;
    i__1 = i2;
    i__3 = inc;
    for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
	    ip = ipiv[ix];
	    if (ip != i__) {
        i__2 = *n;
        for (k = n32; k <= i__2; ++k) {
          temp = a[i__ + k * a_dim1];
          a[i__ + k * a_dim1] = a[ip + k * a_dim1];
          a[ip + k * a_dim1] = temp;
/* L40: */
        }
	    }
	    ix += *incx;
/* L50: */
    }
  }

  return 0;

/*     End of AZ_DLASWP */

} /* az_dlaswp */


double az_r_sign(real * a, real * b) {
  double x;
  x = (*a >= 0 ? *a : -*a);
  return (*b >= 0 ? x : -x);
}


/* Subroutine */ int az_slaic1_c(integer *job, integer *j, real *x, real *
	sest, real *w, real *gamma, real *sestpr, real *s, real *c__)
{
  /* System generated locals */
  real r__1, r__2, r__3, r__4;

  /* Builtin functions */
  double sqrt(doublereal);

  /* Local variables */
  static real b, t, s1, s2, eps, tmp, sine;
  extern doublereal SDOT_F77(integer *, real *, integer *, real *, integer *);
  static real test, zeta1, zeta2, alpha, norma, absgam, absalp;
  extern doublereal SLAMCH_F77(char *, ftnlen);
  static real cosine, absest;


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  AZ_SLAIC1 applies one step of incremental condition estimation in */
/*  its simplest version: */

/*  Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j */
/*  lower triangular matrix L, such that */
/*           twonorm(L*x) = sest */
/*  Then AZ_SLAIC1 computes sestpr, s, c such that */
/*  the vector */
/*                  [ s*x ] */
/*           xhat = [  c  ] */
/*  is an approximate singular vector of */
/*                  [ L     0  ] */
/*           Lhat = [ w' gamma ] */
/*  in the sense that */
/*           twonorm(Lhat*xhat) = sestpr. */

/*  Depending on JOB, an estimate for the largest or smallest singular */
/*  value is computed. */

/*  Note that [s c]' and sestpr**2 is an eigenpair of the system */

/*      diag(sest*sest, 0) + [alpha  gamma] * [ alpha ] */
/*                                            [ gamma ] */

/*  where  alpha =  x'*w. */

/*  Arguments */
/*  ========= */

/*  JOB     (input) INTEGER */
/*          = 1: an estimate for the largest singular value is computed. */
/*          = 2: an estimate for the smallest singular value is computed. */

/*  J       (input) INTEGER */
/*          Length of X and W */

/*  X       (input) REAL array, dimension (J) */
/*          The j-vector x. */

/*  SEST    (input) REAL */
/*          Estimated singular value of j by j matrix L */

/*  W       (input) REAL array, dimension (J) */
/*          The j-vector w. */

/*  GAMMA   (input) REAL */
/*          The diagonal element gamma. */

/*  SESTPR  (output) REAL */
/*          Estimated singular value of (j+1) by (j+1) matrix Lhat. */

/*  S       (output) REAL */
/*          Sine needed in forming xhat. */

/*  C       (output) REAL */
/*          Cosine needed in forming xhat. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

  /* Parameter adjustments */
  --w;
  --x;

  /* Function Body */
  eps = SLAMCH_F77("Epsilon", (ftnlen)7);
  alpha = SDOT_F77(j, &x[1], &c__1, &w[1], &c__1);

  absalp = dabs(alpha);
  absgam = dabs(*gamma);
  absest = dabs(*sest);

  if (*job == 1) {

/*        Estimating largest singular value */

/*        special cases */

    if (*sest == 0.f) {
	    s1 = dmax(absgam,absalp);
	    if (s1 == 0.f) {
        *s = 0.f;
        *c__ = 1.f;
        *sestpr = 0.f;
	    } else {
        *s = alpha / s1;
        *c__ = *gamma / s1;
        tmp = sqrt(*s * *s + *c__ * *c__);
        *s /= tmp;
        *c__ /= tmp;
        *sestpr = s1 * tmp;
	    }
	    return 0;
    } else if (absgam <= eps * absest) {
	    *s = 1.f;
	    *c__ = 0.f;
	    tmp = dmax(absest,absalp);
	    s1 = absest / tmp;
	    s2 = absalp / tmp;
	    *sestpr = tmp * sqrt(s1 * s1 + s2 * s2);
	    return 0;
    } else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
        *s = 1.f;
        *c__ = 0.f;
        *sestpr = s2;
	    } else {
        *s = 0.f;
        *c__ = 1.f;
        *sestpr = s1;
	    }
	    return 0;
    } else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
        tmp = s1 / s2;
        *s = sqrt(tmp * tmp + 1.f);
        *sestpr = s2 * *s;
        *c__ = *gamma / s2 / *s;
        *s = az_r_sign(&c_b19, &alpha) / *s;
	    } else {
        tmp = s2 / s1;
        *c__ = sqrt(tmp * tmp + 1.f);
        *sestpr = s1 * *c__;
        *s = alpha / s1 / *c__;
        *c__ = az_r_sign(&c_b19, gamma) / *c__;
	    }
	    return 0;
    } else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *gamma / absest;

	    b = (1.f - zeta1 * zeta1 - zeta2 * zeta2) * .5f;
	    *c__ = zeta1 * zeta1;
	    if (b > 0.f) {
        t = *c__ / (b + sqrt(b * b + *c__));
	    } else {
        t = sqrt(b * b + *c__) - b;
	    }

	    sine = -zeta1 / t;
	    cosine = -zeta2 / (t + 1.f);
	    tmp = sqrt(sine * sine + cosine * cosine);
	    *s = sine / tmp;
	    *c__ = cosine / tmp;
	    *sestpr = sqrt(t + 1.f) * absest;
	    return 0;
    }

  } else if (*job == 2) {

/*        Estimating smallest singular value */

/*        special cases */

    if (*sest == 0.f) {
	    *sestpr = 0.f;
	    if (dmax(absgam,absalp) == 0.f) {
        sine = 1.f;
        cosine = 0.f;
	    } else {
        sine = -(*gamma);
        cosine = alpha;
	    }
/* Computing MAX */
	    r__1 = dabs(sine), r__2 = dabs(cosine);
	    s1 = dmax(r__1,r__2);
	    *s = sine / s1;
	    *c__ = cosine / s1;
	    tmp = sqrt(*s * *s + *c__ * *c__);
	    *s /= tmp;
	    *c__ /= tmp;
	    return 0;
    } else if (absgam <= eps * absest) {
	    *s = 0.f;
	    *c__ = 1.f;
	    *sestpr = absgam;
	    return 0;
    } else if (absalp <= eps * absest) {
	    s1 = absgam;
	    s2 = absest;
	    if (s1 <= s2) {
        *s = 0.f;
        *c__ = 1.f;
        *sestpr = s1;
	    } else {
        *s = 1.f;
        *c__ = 0.f;
        *sestpr = s2;
	    }
	    return 0;
    } else if (absest <= eps * absalp || absest <= eps * absgam) {
	    s1 = absgam;
	    s2 = absalp;
	    if (s1 <= s2) {
        tmp = s1 / s2;
        *c__ = sqrt(tmp * tmp + 1.f);
        *sestpr = absest * (tmp / *c__);
        *s = -(*gamma / s2) / *c__;
        *c__ = az_r_sign(&c_b19, &alpha) / *c__;
	    } else {
        tmp = s2 / s1;
        *s = sqrt(tmp * tmp + 1.f);
        *sestpr = absest / *s;
        *c__ = alpha / s1 / *s;
        *s = -az_r_sign(&c_b19, gamma) / *s;
	    }
	    return 0;
    } else {

/*           normal case */

	    zeta1 = alpha / absest;
	    zeta2 = *gamma / absest;

/* Computing MAX */
	    r__3 = zeta1 * zeta1 + 1.f + (r__1 = zeta1 * zeta2, dabs(r__1)), 
		    r__4 = (r__2 = zeta1 * zeta2, dabs(r__2)) + zeta2 * zeta2;
	    norma = dmax(r__3,r__4);

/*           See if root is closer to zero or to ONE */

	    test = (zeta1 - zeta2) * 2.f * (zeta1 + zeta2) + 1.f;
	    if (test >= 0.f) {

/*              root is close to zero, compute directly */

        b = (zeta1 * zeta1 + zeta2 * zeta2 + 1.f) * .5f;
        *c__ = zeta2 * zeta2;
        t = *c__ / (b + sqrt((r__1 = b * b - *c__, dabs(r__1))));
        sine = zeta1 / (1.f - t);
        cosine = -zeta2 / t;
        *sestpr = sqrt(t + eps * 4.f * eps * norma) * absest;
	    } else {

/*              root is closer to ONE, shift by that amount */

        b = (zeta2 * zeta2 + zeta1 * zeta1 - 1.f) * .5f;
        *c__ = zeta1 * zeta1;
        if (b >= 0.f) {
          t = -(*c__) / (b + sqrt(b * b + *c__));
        } else {
          t = b - sqrt(b * b + *c__);
        }
        sine = -zeta1 / t;
        cosine = -zeta2 / (t + 1.f);
        *sestpr = sqrt(t + 1.f + eps * 4.f * eps * norma) * absest;
	    }
	    tmp = sqrt(sine * sine + cosine * cosine);
	    *s = sine / tmp;
	    *c__ = cosine / tmp;
	    return 0;

    }
  }
  return 0;

/*     End of AZ_SLAIC1 */

} /* az_slaic1 */

/* Subroutine */ int az_slaswp_c(integer *n, real *a, integer *lda, integer *
	k1, integer *k2, integer *ipiv, integer *incx)
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1;

  /* Local variables */
  static integer i__, ip, ix;
  extern /* Subroutine */ int sswap_c(integer *, real *, integer *, real *, 
    integer *);


/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     October 31, 1992 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  AZ_SLASWP performs a series of row interchanges on the matrix A. */
/*  One row interchange is initiated for each of rows K1 through K2 of A. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A. */

/*  A       (input/output) REAL array, dimension (LDA,N) */
/*          On entry, the matrix of column dimension N to which the row */
/*          interchanges will be applied. */
/*          On exit, the permuted matrix. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. */

/*  K1      (input) INTEGER */
/*          The first element of IPIV for which a row interchange will */
/*          be done. */

/*  K2      (input) INTEGER */
/*          The last element of IPIV for which a row interchange will */
/*          be done. */

/*  IPIV    (input) INTEGER array, dimension (M*abs(INCX)) */
/*          The vector of pivot indices.  Only the elements in positions */
/*          K1 through K2 of IPIV are accessed. */
/*          IPIV(K) = L implies rows K and L are to be interchanged. */

/*  INCX    (input) INTEGER */
/*          The increment between successive values of IPIV.  If IPIV */
/*          is negative, the pivots are applied in reverse order. */

/* ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Interchange row I with row IPIV(I) for each of rows K1 through K2. */

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --ipiv;

  /* Function Body */
  if (*incx == 0) {
    return 0;
  }
  if (*incx > 0) {
    ix = *k1;
  } else {
    ix = (1 - *k2) * *incx + 1;
  }
  if (*incx == 1) {
    i__1 = *k2;
    for (i__ = *k1; i__ <= i__1; ++i__) {
	    ip = ipiv[i__];
	    if (ip != i__) {
        SSWAP_F77(n, &a[i__ + a_dim1], lda, &a[ip + a_dim1], lda);
	    }
/* L10: */
    }
  } else if (*incx > 1) {
    i__1 = *k2;
    for (i__ = *k1; i__ <= i__1; ++i__) {
	    ip = ipiv[ix];
	    if (ip != i__) {
        SSWAP_F77(n, &a[i__ + a_dim1], lda, &a[ip + a_dim1], lda);
	    }
	    ix += *incx;
/* L20: */
    }
  } else if (*incx < 0) {
    i__1 = *k1;
    for (i__ = *k2; i__ >= i__1; --i__) {
	    ip = ipiv[ix];
	    if (ip != i__) {
        SSWAP_F77(n, &a[i__ + a_dim1], lda, &a[ip + a_dim1], lda);
	    }
	    ix += *incx;
/* L30: */
    }
  }

  return 0;

/*     End of AZ_SLASWP */

} /* az_slaswp */

