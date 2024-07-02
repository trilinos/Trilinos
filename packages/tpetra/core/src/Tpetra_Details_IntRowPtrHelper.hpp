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
// ************************************************************************
// @HEADER

#ifndef TPETRA_INTROWPTR_HELPER_HPP
#define TPETRA_INTROWPTR_HELPER_HPP

#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Kokkos_Core.hpp"
#include "KokkosSparse_spmv_handle.hpp"

namespace Tpetra {
namespace Details {

/// Helper for CrsMatrix BlockCrsMatrix. Converts rowptrs
/// of the local matrix to int if it's not already int-valued, and overflow won't occur.
/// This enables TPLs to be used in KokkosSparse::spmv and KokkosSparse::spgemm
///
/// This functionality is temporary. When KokkosKernels 4.4 is released, the default offset (rowptrs element type)
/// in KK's ETI will switch to int, and Tpetra's matrix types will also switch to int. This will make the int rowptrs handling
/// here unnecessary, and CrsMatrix/BlockCrsMatrix can simply store the SPMVHandle themselves.
/// \tparam IntLocalMatrix: if LocalMatrix is a CrsMatrix, this parameter is automatically a CrsMatrix with integer-typed rows pointers. Otherwise (BlockCrsMatrix), a corresponding integer-typed BlockCrsMatrix will need to be provided by the 
template<typename LocalMatrix, 
         typename IntLocalMatrix = KokkosSparse::CrsMatrix<typename LocalMatrix::value_type,
                              typename LocalMatrix::ordinal_type,
                              typename LocalMatrix::device_type,
                              void, int>>
struct IntRowPtrHelper
{
  static_assert(is_crs_matrix_v<LocalMatrix> == is_crs_matrix_v<IntLocalMatrix>, "Expected both to be CrsMatrix or BlockCrsMatrix");

  using Rowptrs = typename LocalMatrix::row_map_type;
  using IntRowptrs = typename IntLocalMatrix::row_map_type;

  IntRowPtrHelper(size_t nnz, const typename LocalMatrix::row_map_type& rowptrs)
    : nnz_(nnz)
  {
    if(shouldUseIntRowptrs())
      fillRowptrsInt(rowptrs);
    // otherwise, leave rowptrs_int empty
  }

  bool shouldUseIntRowptrs() const
  {
    // Use int-valued rowptrs if:
    // - number of entries doesn't overflow int
    // - LocalMatrix's rowptrs is not already int-typed
    return nnz_ <= size_t(INT_MAX) && !std::is_same_v<int, typename Rowptrs::non_const_value_type>;
  }

  void fillRowptrsInt(const Rowptrs& rowptrs)
  {
    // Populate rowptrs_int as a deep copy of rowptrs, converted to int.
    // This needs to be in its own public function because KOKKOS_LAMBDA cannot be used inside constructor or protected/private
    typename IntRowptrs::non_const_type rowptrs_temp(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowptrs_int"), rowptrs.extent(0));
    Kokkos::parallel_for(Kokkos::RangePolicy<typename LocalMatrix::execution_space>(0, rowptrs.extent(0)),
      KOKKOS_LAMBDA(int i)
      {
        rowptrs_temp(i) = int(rowptrs(i));
      });
    rowptrs_int_ = rowptrs_temp;
  }

  IntLocalMatrix getIntRowptrMatrix(const LocalMatrix& A)
  {
    if constexpr(KokkosSparse::is_crs_matrix_v<LocalMatrix>)
    {
      return IntLocalMatrix("", A.numRows(), A.numCols(), A.nnz(),
          A.values, rowptrs_int_, A.graph.entries);
    }
    else
    {
      return IntLocalMatrix("", A.numRows(), A.numCols(), A.nnz(),
          A.values, rowptrs_int_, A.graph.entries, A.blockDim());
    }
  }

  size_t nnz_;
  IntRowptrs rowptrs_int_;
};
}
}

#endif // TPETRA_INTROWPTR_HELPER_HPP

