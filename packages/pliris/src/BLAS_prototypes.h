/*
//@HEADER
// ************************************************************************
//
//               Pliris: Parallel Dense Solver Package
//                 Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER
*/

/* Fortran BLAS and C BLAS prototypes reside in this file */


/* Fortran BLAS prototypes: */

/* Single Precision Real Fortran BLAS: */
void scopy_(int *n, float *sx, int *incx, float *sy, int *incy);
void sscal_(int *n, float *sa, float *sx, int *incx);
void saxpy_(int *n, float *sa, float *sx, int *incx, float *sy, int *incy);
int isamax_(int *n, float *sx, int *incx);
void sasumbb_(float *sasumb, int *n, float *sx, int *incx);
void sdotbb_(float *sdotb, int *n, float *sx, int *incx, float *sy, int *incy);
void sgemmu_(char *transa, char *transb, int *m, int *n, int *k, float *alpha,
   float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);

/* Double Precision Real Fortran BLAS: */
void dcopy_(int *n, double *sx, int *incx, double *sy, int *incy);
void dscal_(int *n, double *sa, double *sx, int *incx);
void daxpy_(int *n, double *sa, double *sx, int *incx, double *sy, int *incy);
int idamax_(int *n, double *sx, int *incx);
double dasum_( int *n, double *sx, int *incx);
double  ddot_( int *n, double *sx, int *incx, double *sy, int *incy);
void dgemmu_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
   double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
   double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

/* Single Precision Complex Fortran BLAS: */
#ifdef SCPLX
	void ccopy_(int *n, scomplex *cx, int *incx, scomplex *cy, int *incy);
	void cscal_(int *n, scomplex *ca, scomplex *cx, int *incx);
	void caxpy_(int *n, scomplex *ca, scomplex *cx, int *incx, scomplex *cy,
   	   int *incy);
	int icamax_(int *n, scomplex *cx, int *incx);
	void scasumbb_(float *scasumb, int *n, scomplex *cx, int *incx);
	void cdotubb_(scomplex *cdotub, int *n, scomplex *cx, int *incx,
   	   scomplex *cy, int *incy);
	void cgemmu_(char *transa, char *transb, int *m, int *n, int *k,
   	   scomplex *alpha, scomplex *a, int *lda, scomplex *b, int *ldb,
   	   scomplex *beta, scomplex *c, int *ldc);
#endif

/* Double Precision Complex Fortran BLAS: */
#ifdef ZCPLX
	void zcopy_(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);
	void zscal_(int *n, dcomplex *ca, dcomplex *cx, int *incx);
	void zaxpy_(int *n, dcomplex *ca, dcomplex *cx, int *incx, dcomplex *cy,
   	   int *incy);
	int izamax_(int *n, dcomplex *cx, int *incx);
	void szasumbb_(double *scasumb, int *n, dcomplex *cx, int *incx);
	void zdotubb_(dcomplex *cdotub, int *n, dcomplex *cx, int *incx,
   	   dcomplex *cy, int *incy);
	void zgemm_(char *transa, char *transb, int *m, int *n, int *k,
   	   dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b, int *ldb,
   	   dcomplex *beta, dcomplex *c, int *ldc);
#endif





/* CBLAS prototypes: */

/* Single Precision Real CBLAS: */
void sscal(int n, float sa, float *sx, int incx);
void scopy(int n, float *sx, int incx, float *sy, int incy);
void saxpy(int n, float sa, float *sx, int incx, float *sy, int incy);
int isamax(int n, float *sx, int incx);
float sasum(int n, float *sx, int incx);
float sdot(int n, float *sx, int incx, float *sy, int incy);

/* Double Precision Real CBLAS: */
void dscal(int n, double da, double *dx, int incx);
void dcopy(int n, double *dx, int incx, double *dy, int incy);
void daxpy(int n, double da, double *dx, int incx, double *dy, int incy);
int idamax(int n, double *dx, int incx);
double dasum(int n, double *dx, int incx);
double ddot(int n, double *dx, int incx, double *dy, int incy);

/* Single Precision Complex CBLAS: */
#ifdef SCPLX
    	void cscal(int n, scomplex sa, scomplex *sx, int incx);
	void ccopy(int n, scomplex *sx, int incx, scomplex *sy, int incy);
	void caxpy(int n, scomplex sa, scomplex *sx, int incx, scomplex *sy,
		int incy);
	int icamax(int n, scomplex *sx, int incx);
	float scasum(int n, scomplex *sx, int incx);
	scomplex cdotu(int n, scomplex *sx, int incx, scomplex *sy, int incy);
#endif

/* Double Precision Complex CBLAS: */
#ifdef ZCPLX
	void zscal(int n, dcomplex sa, dcomplex *sx, int incx);
	void zcopy(int n, dcomplex *sx, int incx, dcomplex *sy, int incy);
	void zaxpy(int n, dcomplex sa, dcomplex *sx, int incx, dcomplex *sy,
		int incy);
	int izamax(int n, dcomplex *sx, int incx);
	double dzasum(int n, dcomplex *sx, int incx);
	dcomplex zdotu(int n, dcomplex *sx, int incx, dcomplex *sy, int incy);
#endif


#if defined(Pliris_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Pliris package is deprecated"
#endif
#endif
