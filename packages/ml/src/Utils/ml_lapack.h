
/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/********************************************************************* */
/*          Utilities for Aztec/ML users                               */
/********************************************************************* */
#include "ml_defs.h"
#include "ml_utils.h"
#include <stdio.h>
#include <math.h>

#ifndef __MLLAPACK__
#define __MLLAPACK__


#ifdef f2c_i2
/* for -i2 */
typedef short ftnlen;
#else
typedef long ftnlen;
#endif


#ifndef FSUB_TYPE
#  define  FSUB_TYPE int
#  if defined(ncube)
#     define  FSUB_TYPE void
#  endif
#  if defined(paragon)
#     define  FSUB_TYPE void
#  endif
#  if defined(hp)
#     define  FSUB_TYPE void
#  endif
#endif


#ifdef __cplusplus
extern "C" {
#endif
extern FSUB_TYPE MLFORTRAN(daxpy)(int *n, double *da, double *dx,
				  int *incx, double *dy, int *incy);
extern FSUB_TYPE MLFORTRAN(dgetrs)(char *, int *, int *, double *, int *, int *,
                          double *, int *, int *, int);
extern double MLFORTRAN(ddot)(int *n1, double *v1, int *dum11, double *v2, int *dum21);
extern FSUB_TYPE MLFORTRAN(xerbla)(char *, int *);
extern FSUB_TYPE MLFORTRAN(dlaswp)(int *, double *, int *, 
				 int *, int *, int *, int *);
extern FSUB_TYPE MLFORTRAN(dgemm)(char *, char *, int *, int *, 
				int *, double *, double *, int *, double *, 
				int *, double *, double *, int *,int,int);
extern FSUB_TYPE MLFORTRAN(dtrsm)(char *, char *, char *, char *, 
				int *, int *, double *, double *, int *, 
				double *, int *,int,int,int,int);
extern FSUB_TYPE MLFORTRAN(dgetf2)(int *, int *, double *, int *, int *, int *);
extern FSUB_TYPE MLFORTRAN(dswap)(int *, double *, int *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dger)(int *, int *, double *, double *, int *, 
			       double *, int *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dscal)(int *, double *, double *, int *);
extern FSUB_TYPE MLFORTRAN(idamax)(int *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dgeqr2)(int *, int *, double *, 
				 int *, double *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dlarfb)(char *, char *, char *, char *, 
				 int *, int *, int *, double *, int *, 
				 double *, int *, double *, int *, double *, 
				 int *);
extern FSUB_TYPE MLFORTRAN(dlarft)(char *, char *, int *, int *, 
				 double *, int *, double *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dlarf)(char *, int *, int *, double *, int *, 
				double *, double *, int *, double *);
extern FSUB_TYPE MLFORTRAN(dlarfg)(int *, double *, double *, int *, double *);
extern FSUB_TYPE MLFORTRAN(dgemv)(char *, int *, int *, double *, double *, int *, 
				double *, int *, double *, double *, int *,int);
extern FSUB_TYPE MLFORTRAN(dtrmv)(char *, char *, char *, int *, double *, int *, 
				double *, int *);
extern FSUB_TYPE MLFORTRAN(dcopy)(int *, double *, int *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dtrmm)(char *, char *, char *, char *, int *, int *, double *, 
				double *, int *, double *, int *,int,int,int,int);
extern FSUB_TYPE MLFORTRAN(dlamc2)(int *, int *, long int *, 
				 double *, int *, double *, int *, double *);
extern FSUB_TYPE MLFORTRAN(dlamc1)(int *, int *, long int *, long int *);
extern FSUB_TYPE MLFORTRAN(dlamc4)(int *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dlamc5)(int *, int *, int *, long int *, int *, 
				 double *);
extern FSUB_TYPE MLFORTRAN(dorg2r)(int *, int *, int *, double *, int *, 
				 double *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dgelq2)(int *, int *, double *, 
				  int *, double *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dlabad)(double *, double *);
extern FSUB_TYPE MLFORTRAN(dgelqf)(int *, int *, double *, int *, double *, 
				 double *, int *, int *);
extern FSUB_TYPE MLFORTRAN(dlascl)(char *, int *, int *, double *, double *, 
				 int *, int *, double *, int *, int *);
extern FSUB_TYPE MLFORTRAN(dgeqrf)(int *, int *, double *, int *, 
				 double *, double *, int *, int *);
extern FSUB_TYPE MLFORTRAN(dlaset)(char *, int *, int *, double *, double *, 
				 double *, int *);
extern FSUB_TYPE MLFORTRAN(dormlq)(char *, char *, int *, int *, int *, double *, 
				 int *, double *, double *, int *, double *, 
				 int *, int *);
extern FSUB_TYPE MLFORTRAN(dormqr)(char *, char *, int *, int *, int *, 
				 double *, int *, double *, double *, int *, 
				 double *, int *, int *);
extern FSUB_TYPE MLFORTRAN(dlassq)(int *, double *, int *, double *, double *);
extern FSUB_TYPE MLFORTRAN(dorm2r)(char *, char *, int *, int *, int *, double *, 
				 int *, double *, double *, int *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dorml2)(char *, char *, int *, int *, int *, double *, 
				 int *, double *, double *, int *, double *, int *);
extern FSUB_TYPE MLFORTRAN(dtrsv)(char *uplo, char *trans, char *diag, int *n, 
				double *a, int *lda, double *x, int *incx);
extern FSUB_TYPE MLFORTRAN(dgetrf)(int *m, int *n, double *a, int *
				 lda, int *ipiv, int *info);
extern FSUB_TYPE MLFORTRAN(dpotrs)(char *uplo, int *n, int *nrhs, double *a, 
				 int *lda, double *b, int *ldb, int * info);
extern FSUB_TYPE MLFORTRAN(dgels)(char *trans, int *m, int *n, int *
				nrhs, double *a, int *lda, double *b, int *ldb, 
				double *work, int *lwork, int *info);
extern FSUB_TYPE MLFORTRAN(dorgqr)(int *m, int *n, int *k, double *
				 a, int *lda, double *tau, double *work, int *lwork, 
				 int *info);
extern FSUB_TYPE MLFORTRAN(dpotf2)(char *uplo, int *n, double *a, int *
		   lda, int *info, int);
extern FSUB_TYPE MLFORTRAN(dsyrk)(char *uplo, char *trans, int *n, int *k, 
	double *alpha, double *a, int *lda, double *beta, 
		  double *c, int *ldc, int, int);
extern FSUB_TYPE MLFORTRAN(dlaic1)(int *job, int *j, double *x, 
	double *sest, double *w, double *gamma, double *
		   sestpr, double *s, double *c);
extern FSUB_TYPE MLFORTRAN(dgetri)(int *n, double *a, int *lda, int 
		   *ipiv, double *work, int *lwork, int *info);
extern FSUB_TYPE MLFORTRAN(dpotrf)(char *uplo, int *n, double *a, int *
		   lda, int *info, int);
extern FSUB_TYPE MLFORTRAN(dtrtri)(char *uplo, char *diag, int *n, double *
		   a, int *lda, int *info, int, int);
extern FSUB_TYPE MLFORTRAN(dtrti2)(char *uplo, char *diag, int *n, double *
		   a, int *lda, int *info, int dummy1, int dummy2);


extern long int MLFORTRAN(lsame)(char *, char *);
extern FSUB_TYPE MLFORTRAN(ilaenv)(int *, char *, char *, int *, int *, 
				 int *, int *, ftnlen, ftnlen);
extern double MLFORTRAN(dlamc3)(double *, double *);
extern double MLFORTRAN(dlange)(char *, int *, int *, double *, int *, double *);
extern double MLFORTRAN(dnrm2)(int *, double *, int *);
extern double MLFORTRAN(dlapy2)(double *, double *);
extern double MLFORTRAN(dlamch)(char *);
extern double MLFORTRAN(dasum)(int*n, double *dx, int *incx);

#ifdef __cplusplus
}
#endif

#endif
