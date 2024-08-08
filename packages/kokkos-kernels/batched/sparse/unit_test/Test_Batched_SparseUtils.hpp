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
#ifndef TEST_BATCHED_SPARSE_HELPER_HPP
#define TEST_BATCHED_SPARSE_HELPER_HPP

namespace KokkosBatched {
template <typename IntView, typename VectorViewType>
void create_tridiagonal_batched_matrices(const int nnz, const int BlkSize, const int N, const IntView &r,
                                         const IntView &c, const VectorViewType &D, const VectorViewType &X,
                                         const VectorViewType &B) {
  Kokkos::Random_XorShift64_Pool<typename VectorViewType::device_type::execution_space> random(13718);
  Kokkos::fill_random(X, random, Kokkos::reduction_identity<typename VectorViewType::value_type>::prod());
  Kokkos::fill_random(B, random, Kokkos::reduction_identity<typename VectorViewType::value_type>::prod());

  auto D_host = Kokkos::create_mirror_view(D);
  auto r_host = Kokkos::create_mirror_view(r);
  auto c_host = Kokkos::create_mirror_view(c);

  r_host(0) = 0;

  int current_col = 0;

  for (int i = 0; i < BlkSize; ++i) {
    r_host(i + 1) = r_host(i) + (i == 0 || i == (BlkSize - 1) ? 2 : 3);
  }
  for (int i = 0; i < nnz; ++i) {
    if (i % 3 == 0) {
      for (int l = 0; l < N; ++l) {
        D_host(l, i) = typename VectorViewType::value_type(2.0);
      }
      c_host(i) = current_col;
      ++current_col;
    } else {
      for (int l = 0; l < N; ++l) {
        D_host(l, i) = typename VectorViewType::value_type(-1.0);
      }
      c_host(i) = current_col;
      if (i % 3 == 1)
        --current_col;
      else
        ++current_col;
    }
  }

  Kokkos::fence();

  Kokkos::deep_copy(D, D_host);
  Kokkos::deep_copy(r, r_host);
  Kokkos::deep_copy(c, c_host);

  Kokkos::fence();
}
}  // namespace KokkosBatched

#endif  // TEST_BATCHED_SPARSE_HELPER_HPP
