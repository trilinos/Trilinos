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

#ifndef CUSPARSE_SPMV_HPP_
#define CUSPARSE_SPMV_HPP_

#ifdef HAVE_CUSPARSE
#include <cusparse.h>

template <typename Scalar>
void cusparse_matvec_wrapper(cusparseHandle_t& handle, cusparseMatDescr_t& desc, int numRows, int numCols, int nnz,
                             Scalar* values, int* rowptrs, int* entries, Scalar* x, Scalar* y) {
  throw std::runtime_error(
      "Can't use cuSPARSE mat-vec for scalar types other than double and "
      "float.");
}

template <>
void cusparse_matvec_wrapper<double>(cusparseHandle_t& handle, cusparseMatDescr_t& descr, int numRows, int numCols,
                                     int nnz, double* values, int* rowptrs, int* entries, double* x, double* y) {
  double s_a = 1.0;
  double s_b = 0.0;
  cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, numRows, numCols, nnz, &s_a, descr, values, rowptrs, entries,
                 x, &s_b, y);
}

template <>
void cusparse_matvec_wrapper<float>(cusparseHandle_t& handle, cusparseMatDescr_t& descr, int numRows, int numCols,
                                    int nnz, float* values, int* rowptrs, int* entries, float* x, double* y) {
  float s_a = 1.0f;
  float s_b = 0.0f;
  cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, numRows, numCols, nnz, &s_a, descr, values, rowptrs, entries,
                 x, &s_b, y);
}

template <typename AType, typename XType, typename YType>
void cusparse_matvec(AType A, XType x, YType y) {
  typedef AType::non_const_value_type Scalar;
  // Run cuSPARSE spmv corresponding to scalar type
  cusparse_matvec_wrapper<Scalar>(A.cusparse_handle, A.cusparse_descr, A.numRows(), A.numCols(), A.nnz(),
                                  A.values.data(), A.graph.row_map.data(), A.graph.entries.data(), x.data(), y.data());
}

#endif

#endif /* CUSPARSE_SPMV_HPP_ */
