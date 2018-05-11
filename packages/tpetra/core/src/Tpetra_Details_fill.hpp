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

#ifndef TPETRA_DETAILS_FILL_HPP
#define TPETRA_DETAILS_FILL_HPP

/// \file Tpetra_Details_fill.hpp
/// \brief Declaration and definition of Tpetra::Details::Blas::fill,
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
#include <type_traits>

namespace Tpetra {
namespace Details {
namespace Blas {
namespace Impl {

//! Wrap std::memset, to avoid exposing unnecessary includes.
void*
memsetWrapper (void* dest, int ch, std::size_t count);

/// \brief Implementation of ::Tpetra::Details::Blas::fill.
template<class ViewType,
         class ValueType,
         class ExecutionSpace,
         class IndexType,
         const int rank = ViewType::Rank>
class Fill {
public:
  static void
  fill  (const ExecutionSpace& execSpace,
         const ViewType& X,
         const ValueType& alpha,
         const IndexType /* numRows */,
         const IndexType /* numCols */ )
  {
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be a built-in integer type.");
    Kokkos::deep_copy (execSpace, X, alpha);
  }
};

//! Specialization for rank-1 Views.
template<class ViewType,
         class ValueType,
         class ExecutionSpace,
         class IndexType>
class Fill<ViewType, ValueType, ExecutionSpace, IndexType, 1> {
public:
  Fill (const ViewType& X, const ValueType& alpha) :
    X_ (X), alpha_ (alpha)
  {
    static_assert (ViewType::Rank == 1,
                   "ViewType must be a rank-1 Kokkos::View.");
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be a built-in integer type.");
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const IndexType& i) const
  {
    X_[i] = alpha_;
  }

  static void
  fill (const ExecutionSpace& execSpace,
        const ViewType& X,
        const ValueType& alpha,
        const IndexType numRows,
        const IndexType /* numCols */)
  {
    typedef Kokkos::RangePolicy<ExecutionSpace, IndexType> range_type;
    Kokkos::parallel_for ("fill", range_type (0, numRows),
                          Fill (X, alpha));
  }

private:
  ViewType X_;
  ValueType alpha_;
};

//! Specialization for rank-2 Views.
template<class ViewType,
         class ValueType,
         class ExecutionSpace,
         class IndexType>
class Fill<ViewType, ValueType, ExecutionSpace, IndexType, 2> {
public:
  Fill (const ViewType& X,
        const ValueType& alpha,
        const IndexType numCols) :
    X_ (X), alpha_ (alpha), numCols_ (numCols)
  {
    static_assert (ViewType::Rank == 2,
                   "ViewType must be a rank-2 Kokkos::View.");
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be a built-in integer type.");
  }

  KOKKOS_INLINE_FUNCTION void
  operator () (const IndexType& i) const
  {
    for (IndexType j = 0; j < numCols_; ++j) {
      X_(i,j) = alpha_;
    }
  }

  static void
  fill (const ExecutionSpace& execSpace,
        const ViewType& X,
        const ValueType& alpha,
        const IndexType numRows,
        const IndexType numCols)
  {
    typedef Kokkos::RangePolicy<ExecutionSpace, IndexType> range_type;
    Kokkos::parallel_for ("fill", range_type (0, numRows),
                          Fill (X, alpha, numCols));
  }

private:
  ViewType X_;
  ValueType alpha_;
  IndexType numCols_;
};

#if defined(KOKKOS_ENABLE_SERIAL)
/// \brief Specialization for ExecutionSpace = Kokkos::Serial
///   and rank = 1.
template<class ViewType,
         class ValueType,
         class IndexType>
struct Fill<ViewType,
            ValueType,
            Kokkos::Serial,
            IndexType,
            1>
{
  static void
  fill (const Kokkos::Serial& /* execSpace */,
        const ViewType& X,
        const ValueType& alpha,
        const IndexType numRows,
        const IndexType /* numCols */ )
  {
    static_assert (ViewType::Rank == 1,
                   "ViewType must be a rank-1 Kokkos::View.");
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be a built-in integer type.");
    using ::Tpetra::Details::Blas::BlasSupportsScalar;
    typedef typename ViewType::non_const_value_type view_value_type;

    // Do sizeof(view_value_type) and taking the address of a
    // value_type instance work correctly with memset?
    constexpr bool podType = BlasSupportsScalar<view_value_type>::value;

    if (podType && X.span_is_contiguous () && alpha == ValueType (0.0)) {
      memsetWrapper (X.data (), 0, X.span () * sizeof (view_value_type));
    }
    else {
      for (IndexType k = 0; k < numRows; ++k) {
        X[k] = alpha;
      }
    }
  }
};
#endif // defined(KOKKOS_ENABLE_SERIAL)

#if defined(KOKKOS_ENABLE_SERIAL)
/// \brief Specialization for ExecutionSpace = Kokkos::Serial
///   and rank = 2.
template<class ViewType,
         class ValueType,
         class IndexType>
struct Fill<ViewType,
            ValueType,
            Kokkos::Serial,
            IndexType,
            2>
{
  static void
  fill (const Kokkos::Serial& /* execSpace */,
        const ViewType& X,
        const ValueType& alpha,
        const IndexType numRows,
        const IndexType numCols)
  {
    static_assert (ViewType::Rank == 2,
                   "ViewType must be a rank-2 Kokkos::View.");
    static_assert (std::is_integral<IndexType>::value,
                   "IndexType must be a built-in integer type.");
    using ::Tpetra::Details::Blas::BlasSupportsScalar;
    typedef typename ViewType::non_const_value_type view_value_type;
    typedef typename ViewType::array_layout array_layout;

    // Do sizeof(view_value_type) and taking the address of a
    // value_type instance work correctly with memset?
    constexpr bool podType = BlasSupportsScalar<view_value_type>::value;

    if (podType && alpha == ValueType (0.0)) {
      if (X.span_is_contiguous ()) {
        memsetWrapper (X.data (), 0, X.span () * sizeof (view_value_type));
      }
      else if (std::is_same<array_layout, Kokkos::LayoutLeft>::value) {
        // Tpetra::MultiVector needs to optimize for LayoutLeft.
        for (IndexType j = 0; j < numCols; ++j) {
          auto X_j = Kokkos::subview (X, Kokkos::ALL (), j);
          memsetWrapper (X_j.data (), 0,
                         X_j.extent (0) * sizeof (view_value_type));
        }
      }
      else {
        Kokkos::deep_copy (X, view_value_type (0.0));
      }
    }
    else {
      Kokkos::deep_copy (X, alpha);
    }
  }
};
#endif // defined(KOKKOS_ENABLE_SERIAL)

} // namespace Impl

//
// SKIP TO HERE FOR THE ACTUAL INTERFACE
//

/// \brief Fill the entries of the given 1-D or 2-D Kokkos::View with
///   the given scalar value alpha.
///
/// \tparam ViewType Kokkos::View specialization.
/// \tparam ValueType Type of the scalar value alpha to assign to each
///   entry of X.
/// \tparam IndexType Type of the index to use in loops.
template<class ViewType,
         class ValueType,
         class IndexType,
         class ExecutionSpace>
void
fill (const ExecutionSpace& execSpace,
      const ViewType& X,
      const ValueType& alpha,
      const IndexType numRows,
      const IndexType numCols)
{
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be a built-in integer type.");
  typedef Impl::Fill<ViewType, ValueType, ExecutionSpace,
    IndexType, ViewType::Rank> impl_type;
  impl_type::fill (execSpace, X, alpha, numRows, numCols);
}

template<class ViewType,
         class ValueType,
         class IndexType,
         class ExecutionSpace>
void
fill (const ExecutionSpace& execSpace,
      const ViewType& X,
      const ValueType& alpha,
      const IndexType numRows,
      const IndexType numCols,
      const size_t whichVectors[])
{
  static_assert (ViewType::Rank == 2, "ViewType must be a rank-2 "
                 "Kokkos::View in order to call the \"whichVectors\" "
                 "specialization of fill.");
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be a built-in integer type.");
  for (IndexType k = 0; k < numCols; ++k) {
    const IndexType j = whichVectors[k];
    auto X_j = Kokkos::subview (X, Kokkos::ALL (), j);
    typedef decltype (X_j) one_d_view_type;
    typedef Impl::Fill<one_d_view_type, ValueType, ExecutionSpace,
      IndexType, 1> impl_type;
    impl_type::fill (execSpace, X_j, alpha, numRows, IndexType (1));
  }
}

} // namespace Blas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_FILL_HPP
