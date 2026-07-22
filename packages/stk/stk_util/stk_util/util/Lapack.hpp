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

#ifndef STK_UTIL_UTIL_Lapack_hpp
#define STK_UTIL_UTIL_Lapack_hpp

#include <stk_util/stk_config.h>
#include "stk_util/util/Fortran.hpp"

extern "C"
{
#if !defined(_MKL_LAPACK_H_)

void SIERRA_FORTRAN(dgels)(const char* trans, const int* m, const int* n,
                           const int* nrhs, double* a, const int* lda, double* b,
                           const int* ldb, double* work, const int* lwork,
                           int* info);

void SIERRA_FORTRAN(dgelsy)(const int* m, const int* n,
                            const int* nrhs, double* a, const int* lda, double* b,
                            const int* ldb, int* jpvt, const double* rcond, int* rank,
                            double* work, const int* lwork, int* info);

void SIERRA_FORTRAN(dgeqrf)(const int* m, const int* n, double* a,
                            const int* lda, double* tau, double* work,
                            const int* lwork, int* info );

void SIERRA_FORTRAN(dgetrf)(const int* m, const int* n, double* a,
                            const int* lda, int* ipiv, int* info );

void SIERRA_FORTRAN(dgetrs)(const char* trans, const int* n, const int* nrhs,
                            const double* a, const int* lda, const int* ipiv,
                            double* b, const int* ldb, int* info );

void SIERRA_FORTRAN(dormqr)(const char* side, const char* tran, const int* m,
                            const int* n, const int* k, double* a,
                            const int* lda, double* tau, double* c,
                            const int* ldc, double* work, const int* lwork,
                            int* info);

void SIERRA_FORTRAN(dgecon)(const char* NORM,const int* N, const double* A, const int* LDA,
                            const double* ANORM, double* RCOND, double* WORK, int* IWORK, int* INFO );

void SIERRA_FORTRAN(dgesvd)(const char* jobu, const char* jobvt, const int* m,
                            const int* n, double* a, const int* lda, double* s,
                            double* u, const int* ldu, double* vt, const int* ldvt,
                            double* work, const int* lwork, int* info );

void SIERRA_FORTRAN(dgeqp3)(int* m, int* n, double* A, int* lda, int* jpvt, double* tau, double* work, int* lwork, int* info);

void SIERRA_FORTRAN(dgeqp3)(int* m, int* n, double* A, int* lda, int* jpvt, double* tau, double* work, int* lwork, int* info);

#endif
}

#endif // STK_UTIL_UTIL_Lapack_hpp
