// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef ARMPL_SPMV_HPP_
#define ARMPL_SPMV_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
#include <armpl.h>

void armpl_matvec_wrapper(armpl_spmat_t A, float* x, float* y) {
  const float alpha = 1.0;
  const float beta  = 0.0;
  armpl_spmv_exec_s(ARMPL_SPARSE_OPERATION_NOTRANS, alpha, A, x, beta, y);
}

void armpl_matvec_wrapper(armpl_spmat_t A, double* x, double* y) {
  const double alpha = 1.0;
  const double beta  = 0.0;
  armpl_spmv_exec_d(ARMPL_SPARSE_OPERATION_NOTRANS, alpha, A, x, beta, y);
}

template <typename AType, typename XType, typename YType>
void armpl_matvec(AType /*A*/, XType x, YType y, spmv_additional_data* data) {
  // using Scalar = typename AType::non_const_value_type;
  // Run armpl spmv corresponding to scalar type
  armpl_matvec_wrapper(data->A, x.data(), y.data());
}

#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL
#endif  // ARMPL_SPMV_HPP_
