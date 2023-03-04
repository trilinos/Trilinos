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
#ifndef TEST_BATCHED_DENSE_HELPER_HPP
#define TEST_BATCHED_DENSE_HELPER_HPP

namespace KokkosBatched {
template <typename MatrixViewType, typename VectorViewType>
void create_tridiagonal_batched_matrices(const MatrixViewType &A,
                                         const VectorViewType &B) {
  Kokkos::Random_XorShift64_Pool<
      typename VectorViewType::device_type::execution_space>
      random(13718);
  Kokkos::fill_random(
      B, random,
      Kokkos::reduction_identity<typename VectorViewType::value_type>::prod());

  auto A_host = Kokkos::create_mirror_view(A);

  const int N       = A.extent(0);
  const int BlkSize = A.extent(1);

  for (int l = 0; l < N; ++l) {
    for (int i = 0; i < BlkSize; ++i) {
      for (int j = i; j < BlkSize; ++j) {
        if (i == j)
          A_host(l, i, j) = typename VectorViewType::value_type(2.0);
        else if (i == j - 1) {
          A_host(l, i, j) = typename VectorViewType::value_type(-1.0);
          A_host(l, j, i) = typename VectorViewType::value_type(-1.0);
        } else {
          A_host(l, i, j) = typename VectorViewType::value_type(0.0);
          A_host(l, j, i) = typename VectorViewType::value_type(0.0);
        }
      }
    }
  }

  Kokkos::fence();

  Kokkos::deep_copy(A, A_host);

  Kokkos::fence();
}
}  // namespace KokkosBatched

#endif  // TEST_BATCHED_DENSE_HELPER_HPP
