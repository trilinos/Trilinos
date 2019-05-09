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
  static_assert (ViewType::Rank == 2, "ViewType must be a rank-2 "
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
