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

#ifndef TPETRA_DETAILS_GEMM_HPP
#define TPETRA_DETAILS_GEMM_HPP

/// \file Tpetra_Details_gemm.hpp
/// \brief Declaration and definition of Tpetra::Details::Blas::gemm,
///   an implementation detail of Tpetra::MultiVector.
///
/// \warning This file, and its contents, are an implementation detail
///   of Tpetra::MultiVector.  Either may disappear or change at any
///   time.
///
/// Search for "SKIP TO HERE FOR THE ACTUAL INTERFACE" (sans quotes)
/// to find the actual interface that Tpetra developers are supposed
/// to use.

#include "Tpetra_Details_Blas.hpp"
#include "Tpetra_Details_defaultGemm.hpp"
#include "Tpetra_Details_mklGemm.hpp"
#include "Tpetra_Details_cublasGemm.hpp"
#include "Tpetra_Details_libGemm.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {
namespace Blas {
namespace Impl {

/// \brief Implementation of ::Tpetra::Details::Blas::gemm.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType,
         const bool canUseBlasLibrary =
           Lib::GemmCanUseLib<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType>::value,
         const bool canUseCublas =
           Cublas::GemmCanUseCublas<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType>::value,
         const bool canUseMkl =
           Mkl::GemmCanUseMkl<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType>::value>
struct Gemm {
  static void
  gemm (const char transA,
        const char transB,
        const CoefficientType& alpha,
        const ViewType1& A,
        const ViewType2& B,
        const CoefficientType& beta,
        const ViewType3& C)
  {
    ::Tpetra::Details::Blas::Default::gemm (transA, transB, alpha, A, B, beta, C);
  }
};

//! Specialization that calls cuBLAS.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType>
struct Gemm<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType,
            true, true, false> {
  static void
  gemm (const char transA,
        const char transB,
        const CoefficientType& alpha,
        const ViewType1& A,
        const ViewType2& B,
        const CoefficientType& beta,
        const ViewType3& C)
  {
    ::Tpetra::Details::Blas::Cublas::gemm (transA, transB, alpha, A, B, beta, C);
  }
};

//! Specialization that calls the MKL.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType>
struct Gemm<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType,
            true, false, true> {
  static void
  gemm (const char transA,
        const char transB,
        const CoefficientType& alpha,
        const ViewType1& A,
        const ViewType2& B,
        const CoefficientType& beta,
        const ViewType3& C)
  {
    ::Tpetra::Details::Blas::Mkl::gemm (transA, transB, alpha, A, B, beta, C);
  }
};

//! Specialization that calls the BLAS library.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType>
struct Gemm<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType,
            true, false, false> {
  static void
  gemm (const char transA,
        const char transB,
        const CoefficientType& alpha,
        const ViewType1& A,
        const ViewType2& B,
        const CoefficientType& beta,
        const ViewType3& C)
  {
    ::Tpetra::Details::Blas::Lib::gemm (transA, transB, alpha, A, B, beta, C);
  }
};

} // namespace Impl

//
// SKIP TO HERE FOR THE ACTUAL INTERFACE
//

/// \brief Tpetra's wrapper for the BLAS' _GEMM (dense matrix-matrix
///   multiply, local to an MPI process).
///
/// \tparam ViewType1 Kokkos::View specialization; type of the A matrix.
/// \tparam ViewType2 Kokkos::View specialization; type of the B matrix.
/// \tparam ViewType3 Kokkos::View specialization; type of the C matrix.
/// \tparam CoefficientType Type of each of the coefficients alpha and beta.
/// \tparam IndexType Type of the index to use in loops.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType = int>
void
gemm (const char transA,
      const char transB,
      const CoefficientType& alpha,
      const ViewType1& A,
      const ViewType2& B,
      const CoefficientType& beta,
      const ViewType3& C)
{
  typedef Impl::Gemm<ViewType1, ViewType2, ViewType3,
    CoefficientType, IndexType> impl_type;
  impl_type::gemm (transA, transB, alpha, A, B, beta, C);
}

} // namespace Blas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GEMM_HPP
