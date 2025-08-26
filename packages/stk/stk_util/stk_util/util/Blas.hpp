// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STK_UTIL_UTIL_Blas_hpp
#define STK_UTIL_UTIL_Blas_hpp

#include <stk_util/stk_config.h>
#include "stk_util/util/Fortran.hpp"

extern "C"
{
void SIERRA_FORTRAN(daxpy)(const int *n, const double *dscale, const double x[], const int *incx, double y[],const int *incy); // y=y+dscale*x
void SIERRA_FORTRAN(dscal)(const int *n, const double *dscale, double *vect, const int *inc); //vect = dscale * vect
double SIERRA_FORTRAN(ddot)(const int * n, const double* x, const int * incx, const double* y, const int * incy); // < x , y >
void SIERRA_FORTRAN(sscal)(const int *n, const float *sscale, float *vect, const int *inc); //vect = sscale * vect
void SIERRA_FORTRAN(dswap)(const int* n, double* d, const int* inc, double* d1, const int* inc1); // switch d1 , d // D.N.E.
double SIERRA_FORTRAN(dasum)(const int * n,const double * x,const int * incx);
int SIERRA_FORTRAN(idamax)(const int *n, const double *vect, const int *inc);
void SIERRA_FORTRAN(saxpy)(const int *n, const float *xscale, const float x[], const int *incx, float y[],const int *incy); // y=y+sscale*x
float SIERRA_FORTRAN(sdot)(const int * n, const float* x, const int * incx, const float* y, const int * incy); // < x , y >
float SIERRA_FORTRAN(sasum)(const int * n,const float * x,const int * incx);
void SIERRA_FORTRAN(sswap)(const int* n, float* s, const int* inc, float* s1, const int* inc1); // switch s1 , s // D.N.E.
int SIERRA_FORTRAN(isamax)(const int *n, const float *vect, const int *inc);

void SIERRA_FORTRAN(dgemm)(const char* transa, const char* transb,
                           const int* m, const int* n, const int* k,
                           const double* alpha, const double* a,
                           const int* lda, const double* b, const int* ldb,
                           const double* beta, double* c, const int* ldc);

void SIERRA_FORTRAN(dgemv)(const char* trans, const int* m, const int* n,
                           const double* alpha, const double* a, const int* lda,
                           const double* x, const int* incx, const double* beta,
                           double* y, const int* incy);

void SIERRA_FORTRAN(dtrsm)(const char *side, const char *uplo, const char *transa, const char *diag,
                           const int *m, const int *n, const double *alpha, const double *a, const int *lda,
                           double *b, const int *ldb);
}

#endif // STK_UTIL_UTIL_Blas_hpp
