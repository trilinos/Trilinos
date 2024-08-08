/*
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
*/

#ifndef KOKKOSPARSE_SPGEMM_NOREUSE_TPL_SPEC_DECL_HPP_
#define KOKKOSPARSE_SPGEMM_NOREUSE_TPL_SPEC_DECL_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#include "KokkosSparse_Utils_cusparse.hpp"
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include "KokkosSparse_Utils_mkl.hpp"
#include "mkl_spblas.h"
#endif

namespace KokkosSparse {
namespace Impl {

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && (CUDA_VERSION >= 11000)

template <typename Matrix, typename MatrixConst>
Matrix spgemm_noreuse_cusparse(const MatrixConst &A, const MatrixConst &B) {
  using Scalar                = typename Matrix::value_type;
  cudaDataType cudaScalarType = Impl::cuda_data_type_from<Scalar>();
  KokkosKernels::Experimental::Controls kkControls;
  cusparseHandle_t cusparseHandle = kkControls.getCusparseHandle();
  cusparseOperation_t op          = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cusparseSpMatDescr_t descr_A, descr_B, descr_C;
  cusparseSpGEMMDescr_t spgemmDescr;
  cusparseSpGEMMAlg_t alg = CUSPARSE_SPGEMM_DEFAULT;
  size_t bufferSize1 = 0, bufferSize2 = 0;
  void *buffer1 = nullptr, *buffer2 = nullptr;
  // A is m*n, B is n*k, C is m*k
  int m            = A.numRows();
  int n            = B.numRows();
  int k            = B.numCols();
  const auto alpha = Kokkos::ArithTraits<Scalar>::one();
  const auto beta  = Kokkos::ArithTraits<Scalar>::zero();
  typename Matrix::row_map_type::non_const_type row_mapC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C rowmap"),
                                                         m + 1);

  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_createDescr(&spgemmDescr));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(
      &descr_A, m, n, A.graph.entries.extent(0), (void *)A.graph.row_map.data(), (void *)A.graph.entries.data(),
      (void *)A.values.data(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, cudaScalarType));

  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(
      &descr_B, n, k, B.graph.entries.extent(0), (void *)B.graph.row_map.data(), (void *)B.graph.entries.data(),
      (void *)B.values.data(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, cudaScalarType));

  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&descr_C, m, k, 0, (void *)row_mapC.data(), nullptr, nullptr,
                                              CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                                              cudaScalarType));

  //----------------------------------------------------------------------
  // query workEstimation buffer size, allocate, then call again with buffer.
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_workEstimation(cusparseHandle, op, op, &alpha, descr_A, descr_B, &beta,
                                                          descr_C, cudaScalarType, alg, spgemmDescr, &bufferSize1,
                                                          nullptr));
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&buffer1, bufferSize1));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_workEstimation(cusparseHandle, op, op, &alpha, descr_A, descr_B, &beta,
                                                          descr_C, cudaScalarType, alg, spgemmDescr, &bufferSize1,
                                                          buffer1));

  //----------------------------------------------------------------------
  // query compute buffer size, allocate, then call again with buffer.

  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_compute(cusparseHandle, op, op, &alpha, descr_A, descr_B, &beta, descr_C,
                                                   cudaScalarType, alg, spgemmDescr, &bufferSize2, nullptr));
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&buffer2, bufferSize2));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_compute(cusparseHandle, op, op, &alpha, descr_A, descr_B, &beta, descr_C,
                                                   cudaScalarType, alg, spgemmDescr, &bufferSize2, buffer2));
  int64_t unused1, unused2, c_nnz;
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMatGetSize(descr_C, &unused1, &unused2, &c_nnz));

  typename Matrix::index_type entriesC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C entries"), c_nnz);
  typename Matrix::values_type valuesC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C values"), c_nnz);

  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseCsrSetPointers(descr_C, (void *)row_mapC.data(), (void *)entriesC.data(), (void *)valuesC.data()));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_compute(cusparseHandle, op, op, &alpha, descr_A, descr_B, &beta, descr_C,
                                                   cudaScalarType, alg, spgemmDescr, &bufferSize2, buffer2));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_copy(cusparseHandle, op, op, &alpha, descr_A, descr_B, &beta, descr_C,
                                                cudaScalarType, alg, spgemmDescr));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(descr_A));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(descr_B));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(descr_C));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_destroyDescr(spgemmDescr));
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(buffer1));
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(buffer2));
  return Matrix("C", m, k, c_nnz, valuesC, row_mapC, entriesC);
}

#define SPGEMM_NOREUSE_DECL_CUSPARSE(SCALAR, MEMSPACE, TPL_AVAIL)                                                   \
  template <>                                                                                                       \
  struct SPGEMM_NOREUSE<KokkosSparse::CrsMatrix<SCALAR, int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, void, int>,    \
                        KokkosSparse::CrsMatrix<const SCALAR, const int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,    \
                                                Kokkos::MemoryTraits<Kokkos::Unmanaged>, const int>,                \
                        KokkosSparse::CrsMatrix<const SCALAR, const int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,    \
                                                Kokkos::MemoryTraits<Kokkos::Unmanaged>, const int>,                \
                        true, TPL_AVAIL> {                                                                          \
    using Matrix      = KokkosSparse::CrsMatrix<SCALAR, int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, void, int>;    \
    using ConstMatrix = KokkosSparse::CrsMatrix<const SCALAR, const int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,    \
                                                Kokkos::MemoryTraits<Kokkos::Unmanaged>, const int>;                \
    static KokkosSparse::CrsMatrix<SCALAR, int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, void, int> spgemm_noreuse(  \
        const ConstMatrix &A, bool, const ConstMatrix &B, bool) {                                                   \
      std::string label = "KokkosSparse::spgemm_noreuse[TPL_CUSPARSE," + Kokkos::ArithTraits<SCALAR>::name() + "]"; \
      Kokkos::Profiling::pushRegion(label);                                                                         \
      Matrix C = spgemm_noreuse_cusparse<Matrix>(A, B);                                                             \
      Kokkos::Profiling::popRegion();                                                                               \
      return C;                                                                                                     \
    }                                                                                                               \
  };

#define SPGEMM_NOREUSE_DECL_CUSPARSE_S(SCALAR, TPL_AVAIL)            \
  SPGEMM_NOREUSE_DECL_CUSPARSE(SCALAR, Kokkos::CudaSpace, TPL_AVAIL) \
  SPGEMM_NOREUSE_DECL_CUSPARSE(SCALAR, Kokkos::CudaUVMSpace, TPL_AVAIL)

SPGEMM_NOREUSE_DECL_CUSPARSE_S(float, true)
SPGEMM_NOREUSE_DECL_CUSPARSE_S(double, true)
SPGEMM_NOREUSE_DECL_CUSPARSE_S(Kokkos::complex<float>, true)
SPGEMM_NOREUSE_DECL_CUSPARSE_S(Kokkos::complex<double>, true)

SPGEMM_NOREUSE_DECL_CUSPARSE_S(float, false)
SPGEMM_NOREUSE_DECL_CUSPARSE_S(double, false)
SPGEMM_NOREUSE_DECL_CUSPARSE_S(Kokkos::complex<float>, false)
SPGEMM_NOREUSE_DECL_CUSPARSE_S(Kokkos::complex<double>, false)

#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
template <typename Matrix, typename MatrixConst>
Matrix spgemm_noreuse_mkl(const MatrixConst &A, const MatrixConst &B) {
  using size_type   = typename Matrix::non_const_size_type;
  using index_type  = typename Matrix::non_const_ordinal_type;
  using scalar_type = typename Matrix::non_const_value_type;
  using ExecSpace   = typename Matrix::execution_space;
  using MKLMatrix   = MKLSparseMatrix<scalar_type>;
  auto m            = A.numRows();
  auto n            = A.numCols();
  auto k            = B.numCols();
  MKLMatrix Amkl(m, n, const_cast<size_type *>(A.graph.row_map.data()),
                 const_cast<index_type *>(A.graph.entries.data()), const_cast<scalar_type *>(A.values.data()));
  MKLMatrix Bmkl(n, k, const_cast<size_type *>(B.graph.row_map.data()),
                 const_cast<index_type *>(B.graph.entries.data()), const_cast<scalar_type *>(B.values.data()));
  sparse_matrix_t C;
  matrix_descr generalDescr;
  generalDescr.type = SPARSE_MATRIX_TYPE_GENERAL;
  generalDescr.mode = SPARSE_FILL_MODE_FULL;
  generalDescr.diag = SPARSE_DIAG_NON_UNIT;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Amkl, Bmkl, &C));
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_order(C));
  MKLMatrix wrappedC(C);
  MKL_INT nrows = 0, ncols = 0;
  MKL_INT *rowmapRaw     = nullptr;
  MKL_INT *entriesRaw    = nullptr;
  scalar_type *valuesRaw = nullptr;
  wrappedC.export_data(nrows, ncols, rowmapRaw, entriesRaw, valuesRaw);
  if (nrows != m || ncols != k)
    throw std::runtime_error(
        "KokkosSparse::spgemm: matrix returned by MKL has incorrect "
        "dimensions");
  MKL_INT c_nnz = rowmapRaw[m];
  Kokkos::View<size_type *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> rowmapRawView(rowmapRaw, m + 1);
  Kokkos::View<index_type *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> entriesRawView(entriesRaw,
                                                                                                        c_nnz);
  Kokkos::View<scalar_type *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> valuesRawView(valuesRaw,
                                                                                                        c_nnz);

  typename Matrix::row_map_type::non_const_type row_mapC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C rowmap"),
                                                         m + 1);
  typename Matrix::index_type entriesC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C entries"), c_nnz);
  typename Matrix::values_type valuesC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C values"), c_nnz);

  Kokkos::deep_copy(ExecSpace(), row_mapC, rowmapRawView);
  Kokkos::deep_copy(ExecSpace(), entriesC, entriesRawView);
  Kokkos::deep_copy(ExecSpace(), valuesC, valuesRawView);
  // Now, done with the copy of C owned by MKL
  wrappedC.destroy();
  return Matrix("C", m, k, c_nnz, valuesC, row_mapC, entriesC);
}

#define SPGEMM_NOREUSE_DECL_MKL(SCALAR, EXEC, TPL_AVAIL)                                                              \
  template <>                                                                                                         \
  struct SPGEMM_NOREUSE<                                                                                              \
      KokkosSparse::CrsMatrix<SCALAR, MKL_INT, Kokkos::Device<EXEC, Kokkos::HostSpace>, void, MKL_INT>,               \
      KokkosSparse::CrsMatrix<const SCALAR, const MKL_INT, Kokkos::Device<EXEC, Kokkos::HostSpace>,                   \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const MKL_INT>,                                \
      KokkosSparse::CrsMatrix<const SCALAR, const MKL_INT, Kokkos::Device<EXEC, Kokkos::HostSpace>,                   \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const MKL_INT>,                                \
      true, TPL_AVAIL> {                                                                                              \
    using Matrix = KokkosSparse::CrsMatrix<SCALAR, MKL_INT, Kokkos::Device<EXEC, Kokkos::HostSpace>, void, MKL_INT>;  \
    using ConstMatrix = KokkosSparse::CrsMatrix<const SCALAR, const MKL_INT, Kokkos::Device<EXEC, Kokkos::HostSpace>, \
                                                Kokkos::MemoryTraits<Kokkos::Unmanaged>, const MKL_INT>;              \
    static KokkosSparse::CrsMatrix<SCALAR, MKL_INT, Kokkos::Device<EXEC, Kokkos::HostSpace>, void, MKL_INT>           \
    spgemm_noreuse(const ConstMatrix &A, bool, const ConstMatrix &B, bool) {                                          \
      std::string label = "KokkosSparse::spgemm_noreuse[TPL_MKL," + Kokkos::ArithTraits<SCALAR>::name() + "]";        \
      Kokkos::Profiling::pushRegion(label);                                                                           \
      Matrix C = spgemm_noreuse_mkl<Matrix>(A, B);                                                                    \
      Kokkos::Profiling::popRegion();                                                                                 \
      return C;                                                                                                       \
    }                                                                                                                 \
  };

#define SPGEMM_NOREUSE_DECL_MKL_SE(SCALAR, EXEC) \
  SPGEMM_NOREUSE_DECL_MKL(SCALAR, EXEC, true)    \
  SPGEMM_NOREUSE_DECL_MKL(SCALAR, EXEC, false)

#define SPGEMM_NOREUSE_DECL_MKL_E(EXEC)                    \
  SPGEMM_NOREUSE_DECL_MKL_SE(float, EXEC)                  \
  SPGEMM_NOREUSE_DECL_MKL_SE(double, EXEC)                 \
  SPGEMM_NOREUSE_DECL_MKL_SE(Kokkos::complex<float>, EXEC) \
  SPGEMM_NOREUSE_DECL_MKL_SE(Kokkos::complex<double>, EXEC)

#ifdef KOKKOS_ENABLE_SERIAL
SPGEMM_NOREUSE_DECL_MKL_E(Kokkos::Serial)
#endif
#ifdef KOKKOS_ENABLE_OPENMP
SPGEMM_NOREUSE_DECL_MKL_E(Kokkos::OpenMP)
#endif
#endif

}  // namespace Impl
}  // namespace KokkosSparse

#endif
