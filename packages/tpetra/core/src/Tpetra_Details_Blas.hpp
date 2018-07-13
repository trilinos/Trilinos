// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_DETAILS_BLAS_HPP
#define TPETRA_DETAILS_BLAS_HPP

/// \file Tpetra_Details_Blas.hpp
/// \brief Type traits for Tpetra's BLAS wrappers; an implementation
///   detail of Tpetra::MultiVector.
///
/// \warning This file, and its contents, are an implementation detail
///   of Tpetra::MultiVector.  Either may disappear or change at any
///   time.

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Complex.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {
namespace Blas {

/// \brief Do BLAS libraries (all that are compliant with the BLAS
///   Standard) support the given "scalar" (matrix entry) type?
///
/// Use the class like this:
/// <tt> BlasSupportsScalar<YourScalarType>::value; </tt>
template<class ScalarType>
struct BlasSupportsScalar {
  static constexpr bool value =
    std::is_same<ScalarType, float>::value ||
    std::is_same<ScalarType, double>::value ||
    std::is_same<ScalarType, ::Kokkos::complex<float> >::value ||
    std::is_same<ScalarType, ::Kokkos::complex<double> >::value;
};

/// \brief Do BLAS libraries (all that are compliant with the BLAS
///   Standard) support the given Kokkos array layout?
///
/// Use the class like this:
/// <tt> BlasSupportsLayout<typename SomeKokkosViewType::array_layout>::value; </tt>
template<class LayoutType>
struct BlasSupportsLayout {
  static constexpr bool value =
    std::is_same<LayoutType, ::Kokkos::LayoutLeft>::value;
};

//! Get the stride (leading dimension) of the 2-D Kokkos::View A.
template<class ViewType,
         class IndexType = int>
IndexType
getStride2DView (const ViewType& A)
{
  static_assert (ViewType::Rank == 2, "A must be a rank-2 Kokkos::View.");
  static_assert (std::is_same<typename ViewType::array_layout, Kokkos::LayoutLeft>::value ||
                 std::is_same<typename ViewType::array_layout, Kokkos::LayoutRight>::value ||
                 std::is_same<typename ViewType::array_layout, Kokkos::LayoutStride>::value,
                 "A's layout must be either LayoutLeft, LayoutRight, or LayoutStride.");
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be a built-in integer type.");
  IndexType stride[8];
  A.stride (stride);
  return (A.extent (1) > 1) ? stride[1] : A.extent (0);
}

namespace Impl {

template<class ViewType,
         class ArrayLayout,
         class IndexType>
struct GetStride1DView {
  typedef ArrayLayout array_layout;

  static IndexType getStride (const ViewType& x)
  {
    static_assert (ViewType::Rank == 1, "x must be a rank-1 Kokkos::View.");
    static_assert (std::is_same<typename ViewType::array_layout, Kokkos::LayoutLeft>::value ||
                   std::is_same<typename ViewType::array_layout, Kokkos::LayoutRight>::value ||
                   std::is_same<typename ViewType::array_layout, Kokkos::LayoutStride>::value,
                   "x's layout must be either LayoutLeft, LayoutRight, or LayoutStride.");
    static_assert (std::is_same<typename ViewType::array_layout, array_layout>::value,
                   "ViewType::array_layout must be the same as array_layout.");
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be a built-in integer type.");
    IndexType stride[8];
    x.stride (stride);
    return stride[0];
  }
};

template<class ViewType,
         class IndexType>
struct GetStride1DView<ViewType, Kokkos::LayoutLeft, IndexType> {
  typedef Kokkos::LayoutLeft array_layout;

  static IndexType getStride (const ViewType&)
  {
    static_assert (ViewType::Rank == 1, "x must be a rank-1 Kokkos::View.");
    static_assert (std::is_same<typename ViewType::array_layout, array_layout>::value,
                   "ViewType::array_layout must be the same as array_layout.");
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be a built-in integer type.");
    return static_cast<IndexType> (1);
  }
};

template<class ViewType,
         class IndexType>
struct GetStride1DView<ViewType, Kokkos::LayoutRight, IndexType> {
  typedef Kokkos::LayoutRight array_layout;

  static IndexType getStride (const ViewType&)
  {
    static_assert (ViewType::Rank == 1, "x must be a rank-1 Kokkos::View.");
    static_assert (std::is_same<typename ViewType::array_layout, array_layout>::value,
                   "ViewType::array_layout must be the same as array_layout.");
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be a built-in integer type.");
    return static_cast<IndexType> (1);
  }
};

} // namespace Impl

//! Get the stride ("INCX" in BLAS terms) of the 1-D Kokkos::View x.
template<class ViewType,
         class IndexType = int>
IndexType
getStride1DView (const ViewType& x)
{
  typedef Impl::GetStride1DView<ViewType, typename ViewType::array_layout, IndexType> impl_type;
  return impl_type::getStride (x);
}

} // namespace Blas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_BLAS_HPP
