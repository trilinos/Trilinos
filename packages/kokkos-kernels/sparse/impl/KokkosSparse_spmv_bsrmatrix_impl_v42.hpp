//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSSPARSE_BSRMATRIX_SPMV_IMPL_V42_HPP
#define KOKKOSSPARSE_BSRMATRIX_SPMV_IMPL_V42_HPP

#include <Kokkos_Core.hpp>

#include <KokkosKernels_ViewUtils.hpp>

namespace KokkosSparse {
namespace Impl {

/* One thread for each entry in the product multivector

   Each thread accumulates the partial products for its entry, and writes it
   out.
*/
template <typename Alpha, typename AMatrix, typename XVector, typename Beta, typename YVector>
class BsrSpmvV42NonTrans {
  Alpha alpha_;
  AMatrix a_;
  XVector x_;
  Beta beta_;
  YVector y_;

 public:
  BsrSpmvV42NonTrans(const Alpha &alpha, const AMatrix &a, const XVector &x, const Beta &beta, const YVector &y)
      : alpha_(alpha), a_(a), x_(x), beta_(beta), y_(y) {}

  template <unsigned BLOCK_SIZE = 0>
  KOKKOS_INLINE_FUNCTION void impl(const size_t k) const {
    using a_ordinal_type   = typename AMatrix::non_const_ordinal_type;
    using a_size_type      = typename AMatrix::non_const_size_type;
    using y_value_type     = typename YVector::non_const_value_type;
    using const_block_type = typename AMatrix::const_block_type;

    const a_ordinal_type irhs = k / y_.extent(0);
    const a_ordinal_type row  = k % y_.extent(0);

    // scale by beta
    if (0 == beta_) {
      y_(row, irhs) = 0;  // convert NaN to 0
    } else if (1 != beta_) {
      y_(row, irhs) *= beta_;
    }

    // for non-zero template instantiations,
    // constant propagation should optimize divmod
    a_ordinal_type blocksz;
    if constexpr (0 == BLOCK_SIZE) {
      blocksz = a_.blockDim();
    } else {
      blocksz = BLOCK_SIZE;
    }

    if (0 != alpha_) {
      const a_ordinal_type blockRow = row / blocksz;
      const a_ordinal_type lclrow   = row % blocksz;
      y_value_type accum            = 0;
      const a_size_type j_begin     = a_.graph.row_map(blockRow);
      const a_size_type j_end       = a_.graph.row_map(blockRow + 1);
      for (a_size_type j = j_begin; j < j_end; ++j) {
        const_block_type b            = a_.unmanaged_block_const(j);
        const a_ordinal_type blockcol = a_.graph.entries(j);
        const a_ordinal_type x_start  = blockcol * blocksz;

        const auto x_lcl = Kokkos::subview(x_, Kokkos::make_pair(x_start, x_start + blocksz), irhs);
        for (a_ordinal_type i = 0; i < blocksz; ++i) {
          accum += b(lclrow, i) * x_lcl(i);
        }
      }
      y_(row, irhs) += alpha_ * accum;
    }
  }

  KOKKOS_INLINE_FUNCTION void operator()(const size_t k) const {
    if (false) {
    }
    // clang-format off
    else if ( 1 == a_.blockDim()) { impl< 1>(k); }
    else if ( 2 == a_.blockDim()) { impl< 2>(k); }
    else if ( 3 == a_.blockDim()) { impl< 3>(k); }
    else if ( 4 == a_.blockDim()) { impl< 4>(k); }
    else if ( 5 == a_.blockDim()) { impl< 5>(k); }
    else if ( 6 == a_.blockDim()) { impl< 6>(k); }
    else if ( 7 == a_.blockDim()) { impl< 7>(k); }
    else if ( 8 == a_.blockDim()) { impl< 8>(k); }
    else if ( 9 == a_.blockDim()) { impl< 9>(k); }
    else if (10 == a_.blockDim()) { impl<10>(k); }
    else if (11 == a_.blockDim()) { impl<11>(k); }
    // clang-format on
    else {
      impl<0>(k);
    }
  }
};

template <typename Alpha, typename AMatrix, typename XVector, typename Beta, typename YVector>
void apply_v42(const typename AMatrix::execution_space &exec, const Alpha &alpha, const AMatrix &a, const XVector &x,
               const Beta &beta, const YVector &y) {
  using execution_space = typename AMatrix::execution_space;

  Kokkos::RangePolicy<execution_space> policy(exec, 0, y.size());
  if constexpr (YVector::rank == 1) {
    // Implementation expects a 2D view, so create an unmanaged 2D view
    // with extent 1 in the second dimension
    using Y2D = KokkosKernels::Impl::with_unmanaged_t<
        Kokkos::View<typename YVector::value_type *[1], typename YVector::array_layout, typename YVector::device_type,
                     typename YVector::memory_traits>>;
    using X2D = KokkosKernels::Impl::with_unmanaged_t<
        Kokkos::View<typename XVector::value_type *[1], typename XVector::array_layout, typename XVector::device_type,
                     typename XVector::memory_traits>>;
    const Y2D yu(y.data(), y.extent(0), 1);
    const X2D xu(x.data(), x.extent(0), 1);
    BsrSpmvV42NonTrans op(alpha, a, xu, beta, yu);
    Kokkos::parallel_for(policy, op);
  } else {
    BsrSpmvV42NonTrans op(alpha, a, x, beta, y);
    Kokkos::parallel_for(policy, op);
  }
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_BSRMATRIX_SPMV_IMPL_V42_HPP
