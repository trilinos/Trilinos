// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include "Tpetra_Details_Blas.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {
namespace Blas {

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
    auto X_j = Kokkos::subview (X, Kokkos::make_pair(IndexType(0), numRows), Kokkos::make_pair(IndexType(0), numCols));
    Kokkos::deep_copy(execSpace, X_j, alpha);
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
  static_assert (ViewType::rank == 2, "ViewType must be a rank-2 "
                 "Kokkos::View in order to call the \"whichVectors\" "
                 "specialization of fill.");
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be a built-in integer type.");
  for (IndexType k = 0; k < numCols; ++k) {
    const IndexType j = whichVectors[k];
    auto X_j = Kokkos::subview (X, Kokkos::make_pair(IndexType(0), numRows), j);
    Kokkos::deep_copy(execSpace, X_j, alpha);
  }
}

} // namespace Blas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_FILL_HPP
