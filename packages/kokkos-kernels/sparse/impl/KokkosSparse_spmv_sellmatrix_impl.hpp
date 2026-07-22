// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSSPARSE_SPMV_SELLMATRIX_IMPL_HPP
#define KOKKOSSPARSE_SPMV_SELLMATRIX_IMPL_HPP

#include "KokkosSparse_SellMatrix.hpp"
#include "KokkosKernels_ArithTraits.hpp"

namespace KokkosSparse {
namespace Impl {

template <class ExecutionSpace, class AMatrix, class XVector, class YVector>
void spmv_sellmatrix_impl(const ExecutionSpace& space, [[maybe_unused]] const char mode[],
                          const typename YVector::non_const_value_type& alpha, const AMatrix& A, const XVector& x,
                          const typename YVector::non_const_value_type& beta, const YVector& y) {
  using ordinal_type = typename AMatrix::non_const_ordinal_type;
  using y_value_type = typename YVector::value_type;

  const ordinal_type row_length = A.sell_nnz / A.num_rows;
  Kokkos::parallel_for(
      "KokkosSparse::spmv SellMatrix", Kokkos::RangePolicy<ExecutionSpace, ordinal_type>(space, 0, A.num_rows),
      KOKKOS_LAMBDA(const ordinal_type rowIdx) {
        y_value_type sum = KokkosKernels::ArithTraits<y_value_type>::zero();
        ordinal_type entryIdx, colIdx = 0;
        for (ordinal_type idx = 0; idx < row_length; ++idx) {
          entryIdx = rowIdx + idx * A.num_rows_per_slice;
          if (A.entries(entryIdx) > -1) {
            colIdx = A.entries(entryIdx);
          }
          sum += A.values(entryIdx) * x(colIdx);
        }
        y(rowIdx) = beta * y(rowIdx) + alpha * sum;
      });
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_SPMV_SELLMATRIX_IMPL_HPP
