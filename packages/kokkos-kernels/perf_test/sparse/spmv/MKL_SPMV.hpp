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

#ifndef MKL_SPMV_HPP_
#define MKL_SPMV_HPP_

#ifdef HAVE_MKL
#include <mkl.h>

template <typename Scalar>
void mkl_matvec_wrapper<Scalar>(int numRows, int numCols, int nnz, Scalar* values, int* rowptrs, int* entries,
                                Scalar* x, Scalar* y) {
  throw std::runtime_error(
      "Can't use cuSPARSE mat-vec for scalar types other than double and "
      "float.");
}

template <>
void mkl_matvec_wrapper<double>(int numRows, int numCols, int nnz, double* values, int* rowptrs, int* entries,
                                double* x, double* y) {
  double s_a        = 1.0;
  double s_b        = 0.0;
  char matdescra[6] = "GLNC0";
  char transa       = 'N';
  mkl_dcsrmv(&transa, &numRows, &numCols, &s_a, matdescra, values, entries, rowptrs, rowptrs + 1, x, &s_b, y);
}

template <>
void mkl_matvec_wrapper<float>(int numRows, int numCols, int nnz, float* values, int* rowptrs, int* entries, float* x,
                               float* y) {
  float s_a         = 1.0;
  float s_b         = 0.0;
  char matdescra[6] = "GLNC0";
  char transa       = 'N';
  mkl_scsrmv(&transa, &numRows, &numCols, &s_a, matdescra, values, entries, rowptrs, rowptrs + 1, x, &s_b, y);
}

template <typename AType, typename XType, typename YType>
void mkl_matvec(AType A, XType x, YType y) {
  typedef AType::non_const_value_type Scalar;
  mkl_matvec_wrapper<Scalar>(A.cusparse_handle, A.cusparse_descr, A.numRows(), A.numCols(), A.nnz(), A.values.data(),
                             A.graph.row_map.data(), A.graph.entries.data(), x.data(), y.data());
}
#endif

#endif /* MKL_SPMV_HPP_ */
