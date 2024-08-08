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
#ifndef KOKKOSSPARSE_BSPGEMM_DEBUG_HPP_
#define KOKKOSSPARSE_BSPGEMM_DEBUG_HPP_
#include "KokkosKernels_helpers.hpp"
#include "KokkosBatched_Gemm_Serial_Internal.hpp"
#include <cstring>

namespace KokkosSparse {

namespace Impl {

template <typename data_view_t>
using kk_subview1d = decltype(Kokkos::subview(data_view_t(), Kokkos::make_pair(0, 0)));

// Returns subview
template <typename data_view_t, typename size_type, typename lno_t>
KOKKOS_INLINE_FUNCTION kk_subview1d<data_view_t> get_block(data_view_t data, size_type block_index, lno_t block_size) {
  const auto i = block_index * block_size;
  return Kokkos::subview(data, Kokkos::make_pair(i, i + block_size));
}

template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_, typename ascalar_nnz_view_t_,
          typename blno_row_view_t_, typename blno_nnz_view_t_, typename bscalar_nnz_view_t_, typename clno_row_view_t_,
          typename clno_nnz_view_t_, typename cscalar_nnz_view_t_>
void bspgemm_debug_numeric(KernelHandle* /* handle */, typename KernelHandle::nnz_lno_t m,
                           typename KernelHandle::nnz_lno_t /* n */, typename KernelHandle::nnz_lno_t k,
                           typename KernelHandle::nnz_lno_t block_dim, alno_row_view_t_ row_mapA,
                           alno_nnz_view_t_ entriesA, ascalar_nnz_view_t_ valuesA,

                           bool /* transposeA */, blno_row_view_t_ row_mapB, blno_nnz_view_t_ entriesB,
                           bscalar_nnz_view_t_ valuesB, bool /* transposeB */, clno_row_view_t_ row_mapC,
                           clno_nnz_view_t_ entriesC, cscalar_nnz_view_t_ valuesC) {
  typename alno_row_view_t_::HostMirror h_rma = Kokkos::create_mirror_view(row_mapA);
  Kokkos::deep_copy(h_rma, row_mapA);
  typename alno_nnz_view_t_::HostMirror h_enta = Kokkos::create_mirror_view(entriesA);
  Kokkos::deep_copy(h_enta, entriesA);
  typename ascalar_nnz_view_t_::HostMirror h_vala = Kokkos::create_mirror_view(valuesA);
  Kokkos::deep_copy(h_vala, valuesA);

  typename blno_row_view_t_::HostMirror h_rmb = Kokkos::create_mirror_view(row_mapB);
  Kokkos::deep_copy(h_rmb, row_mapB);
  typename blno_nnz_view_t_::HostMirror h_entb = Kokkos::create_mirror_view(entriesB);
  Kokkos::deep_copy(h_entb, entriesB);
  typename bscalar_nnz_view_t_::HostMirror h_valb = Kokkos::create_mirror_view(valuesB);
  Kokkos::deep_copy(h_valb, valuesB);
  typename clno_row_view_t_::HostMirror h_rmc = Kokkos::create_mirror_view(row_mapC);
  Kokkos::deep_copy(h_rmc, row_mapC);

  typename clno_nnz_view_t_::HostMirror h_entc    = Kokkos::create_mirror_view(entriesC);
  typename cscalar_nnz_view_t_::HostMirror h_valc = Kokkos::create_mirror_view(valuesC);
  Kokkos::fence();

  typedef typename KernelHandle::nnz_lno_t lno_t;
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::nnz_scalar_t scalar_t;
  typedef KokkosBatched::SerialGemmInternal<KokkosBatched::Algo::Gemm::Unblocked> GEMM;

  const auto block_size = block_dim * block_dim;
  const auto ZERO       = static_cast<scalar_t>(0);
  const auto ONE        = static_cast<scalar_t>(1);

  typename cscalar_nnz_view_t_::HostMirror accumulator("acc", k * block_size);
  Kokkos::deep_copy(accumulator, ZERO);
  Kokkos::fence();
  std::vector<bool> acc_flag(k, false);

  h_rmc(0) = 0;
  for (lno_t i = 0; i < m; ++i) {
    const size_type a_row_begin = h_rma(i);
    const size_type a_row_end   = h_rma(i + 1);
    lno_t a_row_size            = a_row_end - a_row_begin;

    size_type c_row_begin    = h_rmc(i);
    lno_t c_row_size         = h_rmc(i + 1) - c_row_begin;
    lno_t c_row_size_counter = 0;

    for (lno_t j = 0; j < a_row_size; ++j) {
      size_type a_ind             = a_row_begin + j;
      lno_t col                   = h_enta(a_ind);
      auto a_val                  = &h_vala(a_ind * block_size);
      const size_type b_row_begin = h_rmb(col);
      const size_type b_row_end   = h_rmb(col + 1);
      lno_t b_row_size            = b_row_end - b_row_begin;
      for (lno_t z = 0; z < b_row_size; ++z) {
        size_type b_ind = b_row_begin + z;
        lno_t b_col     = h_entb(b_ind);
        auto b_val      = &h_valb(b_ind * block_size);

        if (acc_flag[b_col] == false) {
          acc_flag[b_col]                            = true;
          h_entc(c_row_begin + c_row_size_counter++) = b_col;
        }
        // accumulator(b_col) += a_val * b_val
        auto acc = get_block(accumulator, b_col, block_size);
        GEMM::invoke(block_dim, block_dim, block_dim, ONE, a_val, block_dim, 1, b_val, block_dim, 1, ONE, acc.data(),
                     block_dim, 1);
      }
    }

    // if (i == 0) std::cout << "result_cols" << std::endl;

    for (lno_t j = 0; j < c_row_size; ++j) {
      size_type c_ind  = c_row_begin + j;
      lno_t result_col = h_entc(c_ind);
      auto acc         = get_block(accumulator, result_col, block_size);
      Kokkos::deep_copy(get_block(h_valc, c_ind, block_size), acc);
      Kokkos::deep_copy(acc, ZERO);
      Kokkos::fence();
      acc_flag[result_col] = false;
    }
  }

  Kokkos::deep_copy(entriesC, h_entc);
  Kokkos::deep_copy(valuesC, h_valc);
  Kokkos::fence();
}

}  // namespace Impl
}  // namespace KokkosSparse
#endif
