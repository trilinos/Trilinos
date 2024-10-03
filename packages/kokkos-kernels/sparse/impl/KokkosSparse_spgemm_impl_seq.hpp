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
#ifndef KOKKOSSPARSE_SPGEMM_DEBUG_HPP_
#define KOKKOSSPARSE_SPGEMM_DEBUG_HPP_
#include "KokkosKernels_helpers.hpp"
namespace KokkosSparse {

namespace Impl {

template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_, typename blno_row_view_t_,
          typename blno_nnz_view_t_, typename clno_row_view_t_>
void spgemm_debug_symbolic(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,
                           typename KernelHandle::nnz_lno_t /* n */, typename KernelHandle::nnz_lno_t k,
                           alno_row_view_t_ row_mapA, alno_nnz_view_t_ entriesA,

                           bool /* transposeA */, blno_row_view_t_ row_mapB, blno_nnz_view_t_ entriesB,
                           bool /* transposeB */, clno_row_view_t_ row_mapC) {
  Kokkos::Profiling::pushRegion("KokkosSparse::spgemm_symbolic[NATIVE/DEBUG]");

  typename alno_row_view_t_::HostMirror h_rma = Kokkos::create_mirror_view(row_mapA);
  Kokkos::deep_copy(h_rma, row_mapA);
  typename alno_nnz_view_t_::HostMirror h_enta = Kokkos::create_mirror_view(entriesA);
  Kokkos::deep_copy(h_enta, entriesA);

  typename blno_row_view_t_::HostMirror h_rmb = Kokkos::create_mirror_view(row_mapB);
  Kokkos::deep_copy(h_rmb, row_mapB);
  typename blno_nnz_view_t_::HostMirror h_entb = Kokkos::create_mirror_view(entriesB);
  Kokkos::deep_copy(h_entb, entriesB);
  typename clno_row_view_t_::HostMirror h_rmc = Kokkos::create_mirror_view(row_mapC);
  Kokkos::fence();

  typedef typename KernelHandle::nnz_lno_t lno_t;
  typedef typename KernelHandle::size_type size_type;
  // typedef typename KernelHandle::nnz_scalar_t scalar_t;

  std::vector<bool> acc_flag(k, false);

  std::vector<lno_t> result_c_col_indices(k);

  size_type result_index = 0;

  h_rmc(0) = 0;
  for (lno_t i = 0; i < m; ++i) {
    const size_type a_row_begin = h_rma(i);
    const size_type a_row_end   = h_rma(i + 1);
    lno_t a_row_size            = a_row_end - a_row_begin;
    lno_t row_size              = 0;

    for (lno_t j = 0; j < a_row_size; ++j) {
      size_type a_ind = a_row_begin + j;
      lno_t col       = h_enta(a_ind);
      // scalar_t val = h_vala(a_ind);

      const size_type b_row_begin = h_rmb(col);
      const size_type b_row_end   = h_rmb(col + 1);
      lno_t b_row_size            = b_row_end - b_row_begin;
      for (lno_t z = 0; z < b_row_size; ++z) {
        size_type b_ind = b_row_begin + z;
        lno_t b_col     = h_entb(b_ind);
        // scalar_t b_val = h_valb(b_ind);
        // if (i == 0) std::cout << "\tb col:" <<  b_col << std::endl;
        if (acc_flag[b_col] == false) {
          acc_flag[b_col]                  = true;
          result_c_col_indices[row_size++] = b_col;
        }
      }
    }
    result_index += row_size;
    h_rmc(i + 1) = result_index;
    // size_type c_row_begin = h_rmc(i);

    // if (i == 0) std::cout << "result_cols" << std::endl;

    for (lno_t j = 0; j < row_size; ++j) {
      lno_t result_col     = result_c_col_indices[j];
      acc_flag[result_col] = false;
    }
  }

  handle->get_spgemm_handle()->set_c_nnz(result_index);
  Kokkos::deep_copy(row_mapC, h_rmc);
  Kokkos::fence();
  Kokkos::Profiling::popRegion();
}

template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_, typename ascalar_nnz_view_t_,
          typename blno_row_view_t_, typename blno_nnz_view_t_, typename bscalar_nnz_view_t_, typename clno_row_view_t_,
          typename clno_nnz_view_t_, typename cscalar_nnz_view_t_>
void spgemm_debug_numeric(KernelHandle * /* handle */, typename KernelHandle::nnz_lno_t m,
                          typename KernelHandle::nnz_lno_t /* n */, typename KernelHandle::nnz_lno_t k,
                          alno_row_view_t_ row_mapA, alno_nnz_view_t_ entriesA, ascalar_nnz_view_t_ valuesA,

                          bool /* transposeA */, blno_row_view_t_ row_mapB, blno_nnz_view_t_ entriesB,
                          bscalar_nnz_view_t_ valuesB, bool /* transposeB */, clno_row_view_t_ row_mapC,
                          clno_nnz_view_t_ entriesC, cscalar_nnz_view_t_ valuesC) {
  Kokkos::Profiling::pushRegion("KokkosSparse::spgemm_numeric[NATIVE/DEBUG]");
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

  std::vector<scalar_t> accumulator(k, 0);
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
      scalar_t val                = h_vala(a_ind);
      const size_type b_row_begin = h_rmb(col);
      const size_type b_row_end   = h_rmb(col + 1);
      lno_t b_row_size            = b_row_end - b_row_begin;
      for (lno_t z = 0; z < b_row_size; ++z) {
        size_type b_ind = b_row_begin + z;
        lno_t b_col     = h_entb(b_ind);
        scalar_t b_val  = h_valb(b_ind);

        if (acc_flag[b_col] == false) {
          acc_flag[b_col]                            = true;
          h_entc(c_row_begin + c_row_size_counter++) = b_col;
        }
        accumulator[b_col] += b_val * val;
      }
    }

    // if (i == 0) std::cout << "result_cols" << std::endl;

    for (lno_t j = 0; j < c_row_size; ++j) {
      size_type c_ind         = c_row_begin + j;
      lno_t result_col        = h_entc(c_ind);
      h_valc(c_ind)           = accumulator[result_col];
      accumulator[result_col] = 0;
      acc_flag[result_col]    = false;
    }
  }

  Kokkos::deep_copy(entriesC, h_entc);
  Kokkos::deep_copy(valuesC, h_valc);
  Kokkos::fence();
  Kokkos::Profiling::popRegion();
}

}  // namespace Impl
}  // namespace KokkosSparse
#endif
