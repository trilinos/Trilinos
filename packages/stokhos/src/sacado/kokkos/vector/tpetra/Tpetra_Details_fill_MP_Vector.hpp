// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_FILL_MP_VECTOR_HPP
#define TPETRA_DETAILS_FILL_MP_VECTOR_HPP

#include "Tpetra_Details_fill.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"

namespace Tpetra {
namespace Details {
namespace Blas {

template<class DT, class ... DP,
         class ValueType,
         class IndexType,
         class ExecutionSpace>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<DT,DP...> >::value >::type
fill (const ExecutionSpace& execSpace,
      const Kokkos::View<DT,DP...>& X,
      const ValueType& alpha,
      const IndexType numRows,
      const IndexType numCols)
{
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be a built-in integer type.");
  Kokkos::deep_copy(execSpace, X, alpha);
}

template<class DT, class ... DP,
         class ValueType,
         class IndexType,
         class ExecutionSpace>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View<DT,DP...> >::value >::type
fill (const ExecutionSpace& execSpace,
      const Kokkos::View<DT,DP...>& X,
      const ValueType& alpha,
      const IndexType numRows,
      const IndexType numCols,
      const size_t whichVectors[])
{
  typedef Kokkos::View<DT,DP...> ViewType;
  static_assert (ViewType::rank == 2, "ViewType must be a rank-2 "
                 "Kokkos::View in order to call the \"whichVectors\" "
                 "specialization of fill.");
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be a built-in integer type.");
  for (IndexType k = 0; k < numCols; ++k) {
    const IndexType j = whichVectors[k];
    auto X_j = Kokkos::subview (X, Kokkos::ALL (), j);
    Kokkos::deep_copy(execSpace, X_j, alpha);
  }
}

} // namespace Blas
} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_FILL_MP_VECTOR_HPP
