/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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

#ifndef TPETRA_DETAILS_CUBLASGEMM_HPP
#define TPETRA_DETAILS_CUBLASGEMM_HPP

/// \file Tpetra_Details_cublasGemm.hpp
/// \brief Implementation detail of Tpetra::MultiVector
///
/// \warning This file, and its contents, are an implementation detail
///   of Tpetra::MultiVector.  Either may disappear or change at any
///   time.
///
/// The point of this file is to wrap cuBLAS calls, so that
/// application code is not exposed to cublas.h or cublas_v2.h.  This
/// fixes the following issue relating to conflicts at build time
/// between the old and new cuBLAS API:
///
/// https://github.com/trilinos/Trilinos/issues/1194
///
/// It also generally improves encapsulation.

#include "TpetraCore_config.h"
#include "Kokkos_Complex.hpp"

namespace Tpetra {
namespace Details {
namespace Cublas {

/// \brief Wrapped version of cublasCgemm (v1 API).
///
/// See the cuBLAS documentation for details.
void
cgemm (const char char_transA,
       const char char_transB,
       const int m,
       const int n,
       const int k,
       const ::Kokkos::complex<float>& alpha,
       const ::Kokkos::complex<float> A[],
       const int lda,
       const ::Kokkos::complex<float> B[],
       const int ldb,
       const ::Kokkos::complex<float>& beta,
       ::Kokkos::complex<float> C[],
       const int ldc);

/// \brief Wrapped version of cublasDgemm (v1 API).
///
/// See the cuBLAS documentation for details.
void
dgemm (const char char_transA,
       const char char_transB,
       const int m,
       const int n,
       const int k,
       const double alpha,
       const double A[],
       const int lda,
       const double B[],
       const int ldb,
       const double beta,
       double C[],
       const int ldc);

/// \brief Wrapped version of cublasSgemm (v1 API).
///
/// See the cuBLAS documentation for details.
void
sgemm (const char char_transA,
       const char char_transB,
       const int m,
       const int n,
       const int k,
       const float alpha,
       const float A[],
       const int lda,
       const float B[],
       const int ldb,
       const float beta,
       float C[],
       const int ldc);

/// \brief Wrapped version of cublasZgemm (v1 API).
///
/// See the cuBLAS documentation for details.
void
zgemm (const char char_transA,
       const char char_transB,
       const int m,
       const int n,
       const int k,
       const ::Kokkos::complex<double>& alpha,
       const ::Kokkos::complex<double> A[],
       const int lda,
       const ::Kokkos::complex<double> B[],
       const int ldb,
       const ::Kokkos::complex<double>& beta,
       ::Kokkos::complex<double> C[],
       const int ldc);

} // namespace Cublas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CUBLASGEMM_HPP
