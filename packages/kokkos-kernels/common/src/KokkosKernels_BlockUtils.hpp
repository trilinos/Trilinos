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
#ifndef _KOKKOSKERNELS_BLOCKUTILS_HPP
#define _KOKKOSKERNELS_BLOCKUTILS_HPP

// #include <Kokkos_Atomic.hpp>
// #include <atomic>
#include "KokkosBatched_Gemm_Serial_Internal.hpp"

namespace KokkosSparse {
namespace Impl {

// Initializes block: A = [val, val, val, ....]
template <typename size_type, typename value_type>
KOKKOS_INLINE_FUNCTION void kk_block_init(const size_type block_dim, value_type *dst,
                                          const value_type val = static_cast<value_type>(
                                              0)) {  // Note: replaces __host__ std::fill() not to be called from GPU
  for (auto end = dst + (block_dim * block_dim); dst < end; ++dst) {
    *dst = val;
  }
}

// Initializes block: A = B
template <typename size_type, typename value_type>
KOKKOS_INLINE_FUNCTION void kk_block_set(const size_type block_dim, value_type *dst, const value_type *val) {
  memcpy((void *)dst, val, block_dim * block_dim * sizeof(value_type));
}

// Performs A += B on blocks
template <typename size_type, typename value_type>
KOKKOS_INLINE_FUNCTION void kk_block_add(const size_type block_dim, value_type *dst, const value_type *val) {
  const auto end = dst + block_dim * block_dim;
  while (dst < end) {
    *(dst++) += *(val++);
  }
}

// Performs C += A * B on blocks
// Note: block is assumed to be row-major, dense matrix (no extra padding)
// Note: set clear=true to set C = 0 before increment
template <typename size_type, typename value_type,
          typename DGEMM = KokkosBatched::SerialGemmInternal<KokkosBatched::Algo::Gemm::Unblocked>>
KOKKOS_INLINE_FUNCTION void kk_block_dgemm(const size_type block_dim, value_type *dst, const value_type *valA,
                                           const value_type *valB, const bool clear = false) {
  const auto ZERO = static_cast<value_type>(0);
  const auto ONE  = static_cast<value_type>(1);
  DGEMM::invoke(block_dim, block_dim, block_dim, ONE, valA, block_dim, 1, valB, block_dim, 1, clear ? ZERO : ONE, dst,
                block_dim, 1);
}

// dgemm: C = A * B
template <typename size_type, typename value_type>
KOKKOS_INLINE_FUNCTION void kk_block_set_mul(const size_type block_dim, value_type *c_val, const value_type *a_val,
                                             const value_type *b_val) {
  kk_block_dgemm(block_dim, c_val, a_val, b_val, true);
}

// dgemm: C += A * B
template <typename size_type, typename value_type>
KOKKOS_INLINE_FUNCTION void kk_block_add_mul(const size_type block_dim, value_type *c_val, const value_type *a_val,
                                             const value_type *b_val) {
  kk_block_dgemm(block_dim, c_val, a_val, b_val, false);
}

// Performs C += A * B (dense GEMM) on blocks
// Note: all pointers reference dense row-major blocks (no extra padding)
template <typename size_type, typename value_type>
KOKKOS_INLINE_FUNCTION void kk_vector_block_add_mul(const size_type block_dim, value_type *dst, const value_type *valA,
                                                    const value_type *valB) {
  // NOTE: this should be replaced by batched DGEMM
  //       once atomic increment is supported there
  for (size_type row = 0; row < block_dim; ++row) {
    auto const row_offset = row * block_dim;
    for (size_type col = 0; col < block_dim; ++col) {
      auto v  = &dst[row_offset + col];
      auto vb = valB + col;
      for (const value_type *va = valA + row_offset, *end = va + block_dim; va < end; ++va) {
        Kokkos::atomic_add(v, (*va) * (*vb));
        vb += block_dim;
      }
    }
  }
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif  //  _KOKKOSKERNELS_BLOCKUTILS_HPP
