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

#ifndef _KOKKOSKERNELS_SPARSEUTILS_MKL_HPP
#define _KOKKOSKERNELS_SPARSEUTILS_MKL_HPP

#include "KokkosKernels_config.h"

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL

#include <mkl.h>

namespace KokkosSparse {
namespace Impl {

inline void mkl_internal_safe_call(sparse_status_t mkl_status, const char *name, const char *file = nullptr,
                                   const int line = 0) {
  if (SPARSE_STATUS_SUCCESS != mkl_status) {
    std::ostringstream oss;
    oss << "MKL call \"" << name << "\" at " << file << ":" << line << " encountered error: ";
    switch (mkl_status) {
      case SPARSE_STATUS_NOT_INITIALIZED: oss << "SPARSE_STATUS_NOT_INITIALIZED (empty handle or matrix arrays)"; break;
      case SPARSE_STATUS_ALLOC_FAILED:
        oss << "SPARSE_STATUS_ALLOC_FAILED (internal error: memory allocation "
               "failed)";
        break;
      case SPARSE_STATUS_INVALID_VALUE: oss << "SPARSE_STATUS_INVALID_VALUE (invalid input value)"; break;
      case SPARSE_STATUS_EXECUTION_FAILED:
        oss << "SPARSE_STATUS_EXECUTION_FAILED (e.g. 0-diagonal element for "
               "triangular solver)";
        break;
      case SPARSE_STATUS_INTERNAL_ERROR: oss << "SPARSE_STATUS_INTERNAL_ERROR"; break;
      case SPARSE_STATUS_NOT_SUPPORTED:
        oss << "SPARSE_STATUS_NOT_SUPPORTED (e.g. operation for double "
               "precision doesn't support other types)";
        break;
      default: oss << "unknown (code " << (int)mkl_status << ")"; break;
    }
    oss << '\n';
    Kokkos::abort(oss.str().c_str());
  }
}

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSKERNELS_MKL_SAFE_CALL(call) KokkosSparse::Impl::mkl_internal_safe_call(call, #call, __FILE__, __LINE__)

namespace KokkosSparse {
namespace Impl {

inline sparse_operation_t mode_kk_to_mkl(char mode_kk) {
  switch (toupper(mode_kk)) {
    case 'N': return SPARSE_OPERATION_NON_TRANSPOSE;
    case 'T': return SPARSE_OPERATION_TRANSPOSE;
    case 'H': return SPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    default:;
  }
  throw std::invalid_argument("Invalid mode for MKL (should be one of N, T, H)");
}

template <typename value_type>
struct mkl_is_supported_value_type : std::false_type {};

template <>
struct mkl_is_supported_value_type<float> : std::true_type {};
template <>
struct mkl_is_supported_value_type<double> : std::true_type {};
template <>
struct mkl_is_supported_value_type<Kokkos::complex<float>> : std::true_type {};
template <>
struct mkl_is_supported_value_type<Kokkos::complex<double>> : std::true_type {};

// Helper to:
//  - define the MKL type equivalent to a given Kokkos scalar type
//  - provide an easy implicit conversion to that MKL type
template <typename Scalar>
struct KokkosToMKLScalar {
  static_assert(mkl_is_supported_value_type<Scalar>::value, "Scalar type not supported by MKL");
  using type = Scalar;
  KokkosToMKLScalar(Scalar val_) : val(val_) {}
  operator Scalar() const { return val; }
  Scalar val;
};

template <>
struct KokkosToMKLScalar<Kokkos::complex<float>> {
  using type = MKL_Complex8;
  KokkosToMKLScalar(Kokkos::complex<float> val_) : val(val_) {}
  operator MKL_Complex8() const { return {val.real(), val.imag()}; }
  Kokkos::complex<float> val;
};

template <>
struct KokkosToMKLScalar<Kokkos::complex<double>> {
  using type = MKL_Complex16;
  KokkosToMKLScalar(Kokkos::complex<double> val_) : val(val_) {}
  operator MKL_Complex16() const { return {val.real(), val.imag()}; }
  Kokkos::complex<double> val;
};

template <typename Scalar>
struct KokkosToOneMKLScalar {
  // Note: we happen to use the same set of types in classic MKL and OneMKL.
  // If that changes, update this logic.
  static_assert(mkl_is_supported_value_type<Scalar>::value, "Scalar type not supported by OneMKL");
  using type = Scalar;
};

template <typename Mag>
struct KokkosToOneMKLScalar<Kokkos::complex<Mag>> {
  static_assert(mkl_is_supported_value_type<Kokkos::complex<Mag>>::value, "Scalar type not supported by OneMKL");
  using type = std::complex<Mag>;
};

// MKLSparseMatrix provides thin wrapper around MKL matrix handle
// (sparse_matrix_t) and encapsulates MKL call dispatches related to details
// like value_type, allowing simple client code in kernels.
template <typename value_type>
class MKLSparseMatrix {
  static_assert(mkl_is_supported_value_type<value_type>::value, "Provided value_type type not supported by MKL");
  sparse_matrix_t mtx;

 public:
  inline MKLSparseMatrix(sparse_matrix_t mtx_) : mtx(mtx_) {}

  // Constructs MKL sparse matrix from KK sparse views (m rows x n cols)
  inline MKLSparseMatrix(const MKL_INT num_rows, const MKL_INT num_cols, MKL_INT *xadj, MKL_INT *adj,
                         value_type *values) {}

  // Allows using MKLSparseMatrix directly in MKL calls
  inline operator sparse_matrix_t() const { return mtx; }

  // Exports MKL sparse matrix contents into KK views
  inline void export_data(MKL_INT &num_rows, MKL_INT &num_cols, MKL_INT *&rows_start, MKL_INT *&columns,
                          value_type *&values) {}

  inline void destroy() { KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_destroy(mtx)); }
};

template <>
inline MKLSparseMatrix<float>::MKLSparseMatrix(const MKL_INT rows, const MKL_INT cols, MKL_INT *xadj, MKL_INT *adj,
                                               float *values) {
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_s_create_csr(&mtx, SPARSE_INDEX_BASE_ZERO, rows, cols, xadj, xadj + 1, adj, values));
}

template <>
inline MKLSparseMatrix<double>::MKLSparseMatrix(const MKL_INT rows, const MKL_INT cols, MKL_INT *xadj, MKL_INT *adj,
                                                double *values) {
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_d_create_csr(&mtx, SPARSE_INDEX_BASE_ZERO, rows, cols, xadj, xadj + 1, adj, values));
}

template <>
inline MKLSparseMatrix<Kokkos::complex<float>>::MKLSparseMatrix(const MKL_INT rows, const MKL_INT cols, MKL_INT *xadj,
                                                                MKL_INT *adj, Kokkos::complex<float> *values) {
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_create_csr(&mtx, SPARSE_INDEX_BASE_ZERO, rows, cols, xadj, xadj + 1, adj,
                                                      reinterpret_cast<MKL_Complex8 *>(values)));
}

template <>
inline MKLSparseMatrix<Kokkos::complex<double>>::MKLSparseMatrix(const MKL_INT rows, const MKL_INT cols, MKL_INT *xadj,
                                                                 MKL_INT *adj, Kokkos::complex<double> *values) {
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_create_csr(&mtx, SPARSE_INDEX_BASE_ZERO, rows, cols, xadj, xadj + 1, adj,
                                                      reinterpret_cast<MKL_Complex16 *>(values)));
}

template <>
inline void MKLSparseMatrix<float>::export_data(MKL_INT &num_rows, MKL_INT &num_cols, MKL_INT *&rows_start,
                                                MKL_INT *&columns, float *&values) {
  sparse_index_base_t indexing;
  MKL_INT *rows_end;
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_s_export_csr(mtx, &indexing, &num_rows, &num_cols, &rows_start, &rows_end, &columns, &values));
  if (SPARSE_INDEX_BASE_ZERO != indexing) {
    throw std::runtime_error("Expected zero based indexing in exported MKL sparse matrix\n");
    return;
  }
}

template <>
inline void MKLSparseMatrix<double>::export_data(MKL_INT &num_rows, MKL_INT &num_cols, MKL_INT *&rows_start,
                                                 MKL_INT *&columns, double *&values) {
  sparse_index_base_t indexing;
  MKL_INT *rows_end;
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_d_export_csr(mtx, &indexing, &num_rows, &num_cols, &rows_start, &rows_end, &columns, &values));
  if (SPARSE_INDEX_BASE_ZERO != indexing) {
    throw std::runtime_error("Expected zero based indexing in exported MKL sparse matrix\n");
    return;
  }
}

template <>
inline void MKLSparseMatrix<Kokkos::complex<float>>::export_data(MKL_INT &num_rows, MKL_INT &num_cols,
                                                                 MKL_INT *&rows_start, MKL_INT *&columns,
                                                                 Kokkos::complex<float> *&values) {
  sparse_index_base_t indexing;
  MKL_INT *rows_end;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_export_csr(mtx, &indexing, &num_rows, &num_cols, &rows_start, &rows_end,
                                                      &columns, reinterpret_cast<MKL_Complex8 **>(&values)));
  if (SPARSE_INDEX_BASE_ZERO != indexing) {
    throw std::runtime_error("Expected zero based indexing in exported MKL sparse matrix\n");
    return;
  }
}

template <>
inline void MKLSparseMatrix<Kokkos::complex<double>>::export_data(MKL_INT &num_rows, MKL_INT &num_cols,
                                                                  MKL_INT *&rows_start, MKL_INT *&columns,
                                                                  Kokkos::complex<double> *&values) {
  sparse_index_base_t indexing;
  MKL_INT *rows_end;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_export_csr(mtx, &indexing, &num_rows, &num_cols, &rows_start, &rows_end,
                                                      &columns, reinterpret_cast<MKL_Complex16 **>(&values)));
  if (SPARSE_INDEX_BASE_ZERO != indexing) {
    throw std::runtime_error("Expected zero based indexing in exported MKL sparse matrix\n");
    return;
  }
}

}  // namespace Impl
}  // namespace KokkosSparse

// Utilities for oneMKL SYCL code
#ifdef KOKKOS_ENABLE_SYCL
#include "oneapi/mkl/spblas.hpp"

namespace KokkosSparse {
namespace Impl {}
}  // namespace KokkosSparse
#endif  // KOKKOS_ENABLE_SYCL

#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

#endif  // _KOKKOSKERNELS_SPARSEUTILS_MKL_HPP
