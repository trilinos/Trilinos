// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_SPMV_DATA_HPP_
#define KOKKOSKERNELS_SPMV_DATA_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
#include "armpl.h"
#endif

struct spmv_additional_data {
  int test;

#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
  armpl_spmat_t A;
#endif

  spmv_additional_data() = default;

  spmv_additional_data(int test_) : test(test_) {}

#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
  void set_armpl_spmat(int numRows, int numCols, const int* rowptrs, const int* entries, const float* values) {
    armpl_spmat_create_csr_s(&A, numRows, numCols, rowptrs, entries, values, 0);
  }

  void set_armpl_spmat(int numRows, int numCols, const int* rowptrs, const int* entries, const double* values) {
    armpl_spmat_create_csr_d(&A, numRows, numCols, rowptrs, entries, values, 0);
  }
#endif

  ~spmv_additional_data() {
#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
    if (test == 2) {
      armpl_spmat_destroy(A);
    }
#endif
  }
};

#endif /* KOKKOSKERNELS_SPMV_DATA_HPP_ */
