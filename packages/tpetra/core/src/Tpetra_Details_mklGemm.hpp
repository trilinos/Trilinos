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

#ifndef TPETRA_DETAILS_MKLGEMM_HPP
#define TPETRA_DETAILS_MKLGEMM_HPP

/// \file Tpetra_Details_mklGemm.hpp
/// \brief Implementation detail of Tpetra::MultiVector
///
/// \warning This file, and its contents, are an implementation detail
///   of Tpetra::MultiVector.  Either may disappear or change at any
///   time.
///
/// The point of this file is to wrap MKL _GEMM calls, so that
/// application code is not exposed to MKL include files.

#include "TpetraCore_config.h"
#include "Tpetra_Details_Blas.hpp"
#include "Tpetra_Details_libGemm.hpp"

namespace Tpetra {
namespace Details {
namespace Blas {
namespace Mkl {

/// \brief For this set of template parameters, can and should we
///   implement Gemm (see below) using the MKL?
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType = typename ViewType1::non_const_value_type,
         class IndexType = int>
struct GemmCanUseMkl {
#ifdef HAVE_KOKKOSKERNELS_MKL
#  ifdef KOKKOS_ENABLE_CUDA
  static constexpr bool value =
    ::Tpetra::Details::Blas::Lib::GemmCanUseLib<ViewType1, ViewType2, ViewType3,
                          CoefficientType, IndexType>::value &&
    ! std::is_same<typename ViewType1::execution_space, ::Kokkos::Cuda>::value &&
    ! std::is_same<typename ViewType2::execution_space, ::Kokkos::Cuda>::value &&
    ! std::is_same<typename ViewType3::execution_space, ::Kokkos::Cuda>::value;
#  else // NOT KOKKOS_ENABLE_CUDA
  static constexpr bool value =
    ::Tpetra::Details::Blas::Lib::GemmCanUseLib<ViewType1, ViewType2, ViewType3,
                                                CoefficientType, IndexType>::value;
#  endif // KOKKOS_ENABLE_CUDA
#else // NOT HAVE_KOKKOSKERNELS_MKL
  static constexpr bool value = false;
#endif // NOT HAVE_KOKKOSKERNELS_MKL
};

namespace Impl {

/// \brief Wrapped version of MKL's cblas_cgemm.
///
/// See the MKL documentation for details.
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

/// \brief Wrapped version of MKL's cblas_dgemm.
///
/// See the MKL documentation for details.
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

/// \brief Wrapped version of MKL's cblas_sgemm.
///
/// See the MKL documentation for details.
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

/// \brief Wrapped version of MKL's cblas_zgemm.
///
/// See the MKL documentation for details.
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

/// \brief Wrapper for the above wrappers, templated on scalar type
///   (the type of each entry in the matrices).
template<class ScalarType> struct Gemm {};

template<>
struct Gemm< ::Kokkos::complex<float> > {
  typedef ::Kokkos::complex<float> scalar_type;

  static void
  gemm (const char transA,
        const char transB,
        const int m,
        const int n,
        const int k,
        const scalar_type& alpha,
        const scalar_type A[],
        const int lda,
        const scalar_type B[],
        const int ldb,
        const scalar_type& beta,
        scalar_type C[],
        const int ldc)
  {
    return cgemm (transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  }
};

template<>
struct Gemm<double> {
  typedef double scalar_type;

  static void
  gemm (const char transA,
        const char transB,
        const int m,
        const int n,
        const int k,
        const scalar_type& alpha,
        const scalar_type A[],
        const int lda,
        const scalar_type B[],
        const int ldb,
        const scalar_type& beta,
        scalar_type C[],
        const int ldc)
  {
    return dgemm (transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  }
};

template<>
struct Gemm<float> {
  typedef float scalar_type;

  static void
  gemm (const char transA,
        const char transB,
        const int m,
        const int n,
        const int k,
        const scalar_type& alpha,
        const scalar_type A[],
        const int lda,
        const scalar_type B[],
        const int ldb,
        const scalar_type& beta,
        scalar_type C[],
        const int ldc)
  {
    return sgemm (transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  }
};

template<>
struct Gemm< ::Kokkos::complex<double> > {
  typedef ::Kokkos::complex<double> scalar_type;

  static void
  gemm (const char transA,
        const char transB,
        const int m,
        const int n,
        const int k,
        const scalar_type& alpha,
        const scalar_type A[],
        const int lda,
        const scalar_type B[],
        const int ldb,
        const scalar_type& beta,
        scalar_type C[],
        const int ldc)
  {
    return zgemm (transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
  }
};

} // namespace Impl

template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType>
static void
gemm (const char transA,
      const char transB,
      const CoefficientType& alpha,
      const ViewType1& A,
      const ViewType2& B,
      const CoefficientType& beta,
      const ViewType3& C)
{
  typedef CoefficientType scalar_type;
  typedef Impl::Gemm<scalar_type> impl_type;
  typedef int index_type;

  const index_type lda = getStride2DView<ViewType1, index_type> (A);
  const index_type ldb = getStride2DView<ViewType2, index_type> (B);
  const index_type ldc = getStride2DView<ViewType3, index_type> (C);

  const index_type m = static_cast<index_type> (C.extent (0));
  const index_type n = static_cast<index_type> (C.extent (1));
  const bool noTransA = (transA == 'N' || transA == 'n');
  const index_type k = static_cast<index_type> (noTransA ?
                                                A.extent (1) :
                                                A.extent (0));
  impl_type::gemm (transA, transB, m, n, k,
                   alpha, A.data (), lda,
                   B.data (), ldb,
                   beta, C.data (), ldc);
}

} // namespace Mkl
} // namespace Blas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_MKLGEMM_HPP
