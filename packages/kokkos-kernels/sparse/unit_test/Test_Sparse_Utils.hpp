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

#ifndef TEST_SPARSE_UTILS_HPP
#define TEST_SPARSE_UTILS_HPP

#include "KokkosSparse_spmv.hpp"

namespace Test {

template <typename crsMat_t, typename vector_t>
vector_t create_random_y_vector(crsMat_t crsMat, vector_t x_vector) {
  vector_t y_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Y VECTOR"),
                    crsMat.numRows());
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
}

template <typename crsMat_t, typename vector_t>
vector_t create_random_y_vector_mv(crsMat_t crsMat, vector_t x_vector) {
  vector_t y_vector(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Y VECTOR"),
                    crsMat.numRows(), x_vector.extent(1));
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 0, y_vector);
  return y_vector;
}

}  // namespace Test

#endif  // TEST_SPARSE_UTILS_HPP
