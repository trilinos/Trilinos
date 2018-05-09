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

#ifndef TPETRA_DETAILS_GEMV_HPP
#define TPETRA_DETAILS_GEMV_HPP

/// \file Tpetra_Details_gemv.hpp
/// \brief Declaration and definition of Tpetra::Details::Blas::gemv,
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
#include "Tpetra_Details_defaultGemv.hpp"
#include "Tpetra_Details_mklGemv.hpp"
#include "Tpetra_Details_cublasGemv.hpp"
#include "Tpetra_Details_libGemv.hpp"
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace Tpetra {
namespace Details {
namespace Blas {
namespace Impl {

/// \brief Implementation of ::Tpetra::Details::Blas::gemv.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType,
         const bool canUseBlasLibrary =
           Lib::GemvCanUseLib<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType>::value,
         const bool canUseCublas =
           Cublas::GemvCanUseCublas<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType>::value,
         const bool canUseMkl =
           Mkl::GemvCanUseMkl<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType>::value>
struct Gemv {
  static void
  gemv (const char trans,
        const CoefficientType& alpha,
        const ViewType1& A,
        const ViewType2& x,
        const CoefficientType& beta,
        const ViewType3& y)
  {
    ::Tpetra::Details::Blas::Default::gemv (trans, alpha, A, x, beta, y);
  }
};

//! Specialization that calls cuBLAS.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType>
struct Gemv<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType,
            true, true, false> {
  static void
  gemv (const char trans,
        const CoefficientType& alpha,
        const ViewType1& A,
        const ViewType2& x,
        const CoefficientType& beta,
        const ViewType3& y)
  {
    ::Tpetra::Details::Blas::Cublas::gemv (trans, alpha, A, x, beta, y);
  }
};

//! Specialization that calls the MKL.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType>
struct Gemv<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType,
            true, false, true> {
  static void
  gemv (const char trans,
        const CoefficientType& alpha,
        const ViewType1& A,
        const ViewType2& x,
        const CoefficientType& beta,
        const ViewType3& y)
  {
    ::Tpetra::Details::Blas::Mkl::gemv (trans, alpha, A, x, beta, y);
  }
};

//! Specialization that calls the BLAS library.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType>
struct Gemv<ViewType1, ViewType2, ViewType3, CoefficientType, IndexType,
            true, false, false> {
  static void
  gemv (const char trans,
        const CoefficientType& alpha,
        const ViewType1& A,
        const ViewType2& x,
        const CoefficientType& beta,
        const ViewType3& y)
  {
    ::Tpetra::Details::Blas::Lib::gemv (trans, alpha, A, x, beta, y);
  }
};

} // namespace Impl

//
// SKIP TO HERE FOR THE ACTUAL INTERFACE
//

/// \brief Tpetra's wrapper for the BLAS' _GEMV (dense matrix-vector
///   multiply, local to an MPI process).
///
/// \tparam ViewType1 Kokkos::View specialization; type of the A matrix.
/// \tparam ViewType2 Kokkos::View specialization; type of the x vector.
/// \tparam ViewType3 Kokkos::View specialization; type of the y vector.
/// \tparam CoefficientType Type of each of the coefficients alpha and beta.
/// \tparam IndexType Type of the index to use in loops.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType = int>
void
gemv (const char trans,
      const CoefficientType& alpha,
      const ViewType1& A,
      const ViewType2& x,
      const CoefficientType& beta,
      const ViewType3& y)
{
  // "Canonicalize" the input View types, in order to avoid excess
  // implementation type instantiations.
  typedef Kokkos::View<typename ViewType1::data_type,
    typename ViewType1::array_layout,
    typename ViewType1::device_type,
    Kokkos::MemoryUnmanaged> A_type;
  typedef Kokkos::View<typename ViewType2::data_type,
    typename ViewType2::array_layout,
    typename ViewType2::device_type,
    Kokkos::MemoryUnmanaged> x_type;
  typedef Kokkos::View<typename ViewType3::data_type,
    typename ViewType3::array_layout,
    typename ViewType3::device_type,
    Kokkos::MemoryUnmanaged> y_type;

  const char prefix[] = "Tpetra::Details::Blas::gemv: ";
  if (trans == 'n' || trans == 'N') {
    if (A.extent (1) != x.extent (0) ||
        A.extent (0) != y.extent (0)) {
      std::ostringstream os;
      os << prefix << "Dimensions don't match.  A is "
         << A.extent (0) << " x " << A.extent (1)
         << ", x is " << x.extent (0)
         << ", y is " << y.extent (0) << std::endl;
      throw std::invalid_argument (os.str ());
    }
  }
  else {
    if (A.extent (0) != x.extent (0) ||
        A.extent (1) != y.extent (0)) {
      std::ostringstream os;
      os << prefix << "Dimensions don't match.  A^T is "
         << A.extent (1) << " x " << A.extent (0)
         << ", x is " << x.extent (0)
         << ", y is " << y.extent (0) << std::endl;
      throw std::invalid_argument (os.str ());
    }
  }

  typedef Impl::Gemv<A_type, x_type, y_type,
    CoefficientType, IndexType> impl_type;
  impl_type::gemv (trans, alpha, A, x, beta, y);
}

} // namespace Blas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GEMV_HPP
