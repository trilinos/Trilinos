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

#include "Tpetra_Details_cublasGemm.hpp"
#include "Kokkos_Macros.hpp"
#ifdef KOKKOS_ENABLE_CUDA
#  include <cublas.h>
#endif // KOKKOS_ENABLE_CUDA
#include <sstream>
#include <stdexcept>

namespace Tpetra {
namespace Details {
namespace Cublas {

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
       const int ldc)
{
#ifdef KOKKOS_ENABLE_CUDA
  ::cublasDgemm (char_transA, char_transB, m, n, k,
                 alpha, A, lda, B, ldb, beta, C, ldc);
  cublasStatus info = cublasGetError ();
  if (info != CUBLAS_STATUS_SUCCESS) {
    std::ostringstream err;
    err << "cublasDgemm failed with status " << info << ".";
    throw std::runtime_error (err.str ());
  }
#else // NOT KOKKOS_ENABLE_CUDA
  throw std::runtime_error ("You must enable CUDA in your Trilinos build in "
                            "order to invoke cuBLAS functions in Tpetra.");
#endif // KOKKOS_ENABLE_CUDA
}

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
       const int ldc)
{
#ifdef KOKKOS_ENABLE_CUDA
  ::cublasSgemm (char_transA, char_transB, m, n, k,
                 alpha, A, lda, B, ldb, beta, C, ldc);
  cublasStatus info = cublasGetError ();
  if (info != CUBLAS_STATUS_SUCCESS) {
    std::ostringstream err;
    err << "cublasSgemm failed with status " << info << ".";
    throw std::runtime_error (err.str ());
  }
#else // NOT KOKKOS_ENABLE_CUDA
  throw std::runtime_error ("You must enable CUDA in your Trilinos build in "
                            "order to invoke cuBLAS functions in Tpetra.");
#endif // KOKKOS_ENABLE_CUDA
}

} // namespace Cublas
} // namespace Details
} // namespace Tpetra
