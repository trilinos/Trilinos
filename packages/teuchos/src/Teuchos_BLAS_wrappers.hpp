// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_BLAS_WRAPPERS_HPP_
#define _TEUCHOS_BLAS_WRAPPERS_HPP_

#include "Teuchos_ConfigDefs.hpp"

/*! \file Teuchos_BLAS_wrappers.h  
    \brief The Templated Petra BLAS Class.

    Lots of cryptic Fortranish #defines here.  More info on how and why 
    to be added at a later date when Paul gets around to it.
*/

/* Define fcd (Fortran Teuchos_fcd descriptor) for non-standard situations */

#if defined(CRAY_T3X) || defined(INTEL_CXML) || defined(INTEL_MKL)


#if defined(CRAY_T3X)

#include <fortran.h>
#define PREFIX
#define Teuchos_fcd fcd 

#define DASUM_F77   F77_FUNC(sasum,DASUM) 
#define DAXPY_F77   F77_FUNC(saxpy,DAXPY)
#define DCOPY_F77   F77_FUNC(scopy,DCOPY)
#define DDOT_F77    F77_FUNC(sdot,DDOT)
#define DNRM2_F77   F77_FUNC(snrm2,DNRM2)
#define DSCAL_F77   F77_FUNC(sscal,DSCAL)
#define IDAMAX_F77  F77_FUNC(isamax,IDAMAX)
#define DGEMV_F77   F77_FUNC(sgemv,DGEMV)
#define DGER_F77    F77_FUNC(sger,DGER)
#define DTRMV_F77   F77_FUNC(strmv,DTRMV)
#define DGEMM_F77   F77_FUNC(sgemm,DGEMM)
#define DSYMM_F77   F77_FUNC(ssymm,DSYMM)
#define DTRMM_F77   F77_FUNC(strmm,DTRMM)
#define DTRSM_F77   F77_FUNC(strsm,DTRSM)

#elif defined(INTEL_CXML)

#define PREFIX __stdcall 
#define Teuchos_fcd char *, unsigned int 

#define DASUM_F77   F77_FUNC(dasum,DASUM)
#define DAXPY_F77   F77_FUNC(daxpy,DAXPY)
#define DCOPY_F77   F77_FUNC(dcopy,DCOPY)
#define DDOT_F77    F77_FUNC(ddot,DDOT)  
#define DNRM2_F77   F77_FUNC(dnrm2,DNRM2)
#define DSCAL_F77   F77_FUNC(dscal,DSCAL)
#define IDAMAX_F77  F77_FUNC(idamax,IDAMAX)
#define DGEMV_F77   F77_FUNC(dgemv,DGEMV)
#define DGER_F77    F77_FUNC(dger,DGER)
#define DTRMV_F77   F77_FUNC(dtrmv,DTRMV)
#define DGEMM_F77   F77_FUNC(dgemm,DGEMM)
#define DSYMM_F77   F77_FUNC(dsymm,DSYMM)
#define DTRMM_F77   F77_FUNC(dtrmm,DTRMM)
#define DTRSM_F77   F77_FUNC(dtrsm,DTRSM)

#elif defined(INTEL_MKL)

#define PREFIX
#define Teuchos_fcd char *

#define DASUM_F77   F77_FUNC(dasum,DASUM)
#define DAXPY_F77   F77_FUNC(daxpy,DAXPY)
#define DCOPY_F77   F77_FUNC(dcopy,DCOPY)
#define DDOT_F77    F77_FUNC(ddot,DDOT)  
#define DNRM2_F77   F77_FUNC(dnrm2,DNRM2)
#define DSCAL_F77   F77_FUNC(dscal,DSCAL)
#define IDAMAX_F77  F77_FUNC(idamax,IDAMAX)
#define DGEMV_F77   F77_FUNC(dgemv,DGEMV)
#define DGER_F77    F77_FUNC(dger,DGER)
#define DTRMV_F77   F77_FUNC(dtrmv,DTRMV)
#define DGEMM_F77   F77_FUNC(dgemm,DGEMM)
#define DSYMM_F77   F77_FUNC(dsymm,DSYMM)
#define DTRMM_F77   F77_FUNC(dtrmm,DTRMM)
#define DTRSM_F77   F77_FUNC(dtrsm,DTRSM)

#endif 

/* All three of these machines use a simple uppercase mangling of Fortran names */

/* if F77_FUNC is defined undefine it because we want to redefine */

#ifdef F77_FUNC
#undef F77_FUNC
#endif


#define F77_FUNC(lcase,UCASE) PREFIX UCASE

#else /* Define Teuchos_fcd for all other machines */

#define PREFIX
#define Teuchos_fcd char * 

/* In the future use autoconf's definition of F77_FUNC */ 

#ifndef HAVE_CONFIG_H

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#ifdef TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE
#define F77_FUNC(lcase,UCASE) lcase
#else /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE not defined*/
#define F77_FUNC(lcase,UCASE) lcase ## _
#endif /* TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE */

#endif /* HAVE_CONFIG_H */

#define DASUM_F77   F77_FUNC(dasum,DASUM)
#define DAXPY_F77   F77_FUNC(daxpy,DAXPY)
#define DCOPY_F77   F77_FUNC(dcopy,DCOPY)
#define DDOT_F77    F77_FUNC(ddot,DDOT)
#define DNRM2_F77   F77_FUNC(dnrm2,DNRM2)
#define DSCAL_F77   F77_FUNC(dscal,DSCAL)
#define IDAMAX_F77  F77_FUNC(idamax,IDAMAX)
#define DGEMV_F77   F77_FUNC(dgemv,DGEMV)
#define DGER_F77    F77_FUNC(dger,DGER)
#define DTRMV_F77   F77_FUNC(dtrmv,DTRMV)
#define DGEMM_F77   F77_FUNC(dgemm,DGEMM) 
#define DSYMM_F77   F77_FUNC(dsymm,DSYMM)
#define DTRMM_F77   F77_FUNC(dtrmm,DTRMM)
#define DTRSM_F77   F77_FUNC(dtrsm,DTRSM)

#endif


/* Explicitly define each F77 name for all BLAS kernels */

#define SSCAL_F77   F77_FUNC(sscal,SSCAL) 
#define SCOPY_F77   F77_FUNC(scopy,SCOPY)
#define SAXPY_F77   F77_FUNC(saxpy,SAXPY)
#define SDOT_F77    F77_FUNC(sdot,SDOT)
#define SNRM2_F77   F77_FUNC(snrm2,SNRM2)
#define SASUM_F77   F77_FUNC(sasum,SASUM)
#define ISAMAX_F77 F77_FUNC(isamax,ISAMAX)

#define SGEMV_F77   F77_FUNC(sgemv,SGEMV)
#define SGER_F77    F77_FUNC(sger,SGER)
#define STRMV_F77   F77_FUNC(strmv,STRMV)
#define SGEMM_F77   F77_FUNC(sgemm,SGEMM)
#define SSYMM_F77   F77_FUNC(ssymm,SSYMM)
#define STRMM_F77   F77_FUNC(strmm,STRMM)
#define STRSM_F77   F77_FUNC(strsm,STRSM)

#ifdef __cplusplus
extern "C" {
#endif


/* Double precision BLAS 1 */
double PREFIX DASUM_F77(int* n, double x[], int* incx);
void PREFIX DAXPY_F77(int* n, double* alpha, double x[], int* incx, double y[], int* incy);
void PREFIX DCOPY_F77(int* n, double *x, int* incx, double *y, int* incy);
double PREFIX DDOT_F77(int* n, double x[], int* incx, double y[], int* incy);
double PREFIX DNRM2_F77(int* n, double x[], int* incx); 
void PREFIX DSCAL_F77(int* n, double* alpha, double *x, int* incx);
int PREFIX IDAMAX_F77(int* n, double *x, int* incx);

/* Single precision BLAS 1 */ 
float PREFIX SASUM_F77(int* n, float x[], int* incx);
void PREFIX SAXPY_F77(int* n, float* alpha, float x[], int* incx, float y[], int* incy);
void PREFIX SCOPY_F77(int* n, float *x, int* incx, float *y, int* incy);
float PREFIX SDOT_F77(int* n, float x[], int* incx, float y[], int* incy);
float PREFIX SNRM2_F77(int* n, float x[], int* incx); 
void PREFIX SSCAL_F77(int* n, float* alpha, float *x, int* incx);
int PREFIX ISAMAX_F77(int* n, float *x, int* incx);

/* Double precision BLAS 2 */
void PREFIX DGEMV_F77(Teuchos_fcd, int* m, int* n, double* alpha, double A[], int* lda,
                 double x[], int* incx, double* beta, double y[], int* incy);
void PREFIX DTRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, int *n, 
                double *a, int *lda, double *x, int *incx); 
void PREFIX DGER_F77(int *m, int *n, double *alpha, double *x, int *incx, double *y,
               int *incy, double *a, int *lda);


/* Single precision BLAS 2 */
void PREFIX SGEMV_F77(Teuchos_fcd, int* m, int* n, float* alpha, float A[], int* lda,
                 float x[], int* incx, float* beta, float y[], int* incy);
void PREFIX STRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, int *n,
                float *a, int *lda, float *x, int *incx); 
void PREFIX SGER_F77(int *m, int *n, float *alpha, float *x, int *incx, float *y,
               int *incy, float *a, int *lda);

/* Double precision BLAS 3 */
void PREFIX DGEMM_F77(Teuchos_fcd, Teuchos_fcd, int *m, int * 
                n, int *k, double *alpha, double *a, int *lda, 
                double *b, int *ldb, double *beta, double *c, int *ldc);
void PREFIX DSYMM_F77(Teuchos_fcd, Teuchos_fcd, int *m, int * n,
                double *alpha, double *a, int *lda, 
                double *b, int *ldb, double *beta, double *c, int *ldc);
void PREFIX DTRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,  
                int *m, int *n, double *alpha, double *a, int * lda, double *b, int *ldb);
void PREFIX DTRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, 
                int *m, int *n, double *alpha, double *a, int *
                lda, double *b, int *ldb);

/* Single precision BLAS 3 */
void PREFIX SGEMM_F77(Teuchos_fcd, Teuchos_fcd, int *m, int *
                n, int *k, float *alpha, float *a, int *lda, 
                float *b, int *ldb, float *beta, float *c, int *ldc);
void PREFIX SSYMM_F77(Teuchos_fcd, Teuchos_fcd, int *m, int * n,
                float *alpha, float *a, int *lda, 
                float *b, int *ldb, float *beta, float *c, int *ldc);
void PREFIX STRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, 
                int *m, int *n, float *alpha, float *a, int * lda, float *b, int *ldb);
void PREFIX STRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd,
                int *m, int *n, float *alpha, float *a, int *
                lda, float *b, int *ldb);


#ifdef __cplusplus
}
#endif

#endif // end of TEUCHOS_BLAS_WRAPPERS_HPP_
