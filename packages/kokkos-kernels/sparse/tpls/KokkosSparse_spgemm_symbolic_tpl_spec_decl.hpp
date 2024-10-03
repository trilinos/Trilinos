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

#ifndef KOKKOSPARSE_SPGEMM_SYMBOLIC_TPL_SPEC_DECL_HPP_
#define KOKKOSPARSE_SPGEMM_SYMBOLIC_TPL_SPEC_DECL_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#include "KokkosSparse_Utils_cusparse.hpp"
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
#include "rocsparse/rocsparse.h"
#include "KokkosSparse_Utils_rocsparse.hpp"
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include "KokkosSparse_Utils_mkl.hpp"
#include "mkl_spblas.h"
#endif

namespace KokkosSparse {
namespace Impl {

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
// NOTE: all versions of cuSPARSE 10.x and 11.x support exactly the same matrix
// types, so there is no ifdef'ing on versions needed in avail. Offset and
// Ordinal must both be 32-bit. Even though the "generic" API lets you specify
// offsets and ordinals independently as either 16, 32 or 64-bit integers,
// cusparse will just fail at runtime if you don't use 32 for both.

#if (CUDA_VERSION >= 11040)
// 11.4+ supports generic API with reuse (full symbolic/numeric separation)
// However, its "symbolic" (cusparseSpGEMMreuse_nnz) does not populate C's
// rowptrs.
template <typename KernelHandle, typename lno_t, typename ConstRowMapType, typename ConstEntriesType,
          typename RowMapType>
void spgemm_symbolic_cusparse(KernelHandle *handle, lno_t m, lno_t n, lno_t k, const ConstRowMapType &row_mapA,
                              const ConstEntriesType &entriesA, const ConstRowMapType &row_mapB,
                              const ConstEntriesType &entriesB, const RowMapType &row_mapC, bool computeRowptrs) {
  // Split symbolic into two sub-phases: handle/buffer setup and nnz(C), and
  // then rowptrs (if requested). That way, calling symbolic once with
  // computeRowptrs=false, and then again with computeRowptrs=true will not
  // duplicate any work.
  if (!handle->is_symbolic_called()) {
    handle->create_cusparse_spgemm_handle(false, false);
    auto h = handle->get_cusparse_spgemm_handle();

    // Follow
    // https://github.com/NVIDIA/CUDALibrarySamples/tree/master/cuSPARSE/spgemm_reuse
    void *buffer1      = nullptr;
    void *buffer2      = nullptr;
    size_t bufferSize1 = 0;
    size_t bufferSize2 = 0;

    // When nnz is not zero, cusparseCreateCsr insists non-null a value pointer,
    // which however is not available in this function. So we fake it with the
    // entries instead. Fortunately, it seems cupsarse does not access that in
    // the symbolic phase.
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&h->descr_A, m, n, entriesA.extent(0), (void *)row_mapA.data(),
                                                (void *)entriesA.data(), (void *)entriesA.data() /*fake*/,
                                                CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                                                h->scalarType));

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&h->descr_B, n, k, entriesB.extent(0), (void *)row_mapB.data(),
                                                (void *)entriesB.data(), (void *)entriesB.data() /*fake*/,
                                                CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                                                h->scalarType));

#if CUDA_VERSION >= 12020
    // at some point cusparseCreateCsr started to need a non-null row-pointer
    // array, even if the operation that consumed the handle doesn't need to
    // read it. This was observed on a system with CUDA 12.2, but it may have
    // started earlier.
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&h->descr_C, m, k, 0, (void *)row_mapC.data(), nullptr, nullptr,
                                                CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                                                h->scalarType));
#else
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&h->descr_C, m, k, 0, nullptr, nullptr, nullptr, CUSPARSE_INDEX_32I,
                                                CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, h->scalarType));
#endif

    //----------------------------------------------------------------------
    // ask bufferSize1 bytes for external memory
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMMreuse_workEstimation(h->cusparseHandle, h->opA, h->opB, h->descr_A,
                                                                 h->descr_B, h->descr_C, h->alg, h->spgemmDescr,
                                                                 &bufferSize1, nullptr));

    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&buffer1, bufferSize1));
    // inspect matrices A and B to understand the memory requirement for the
    // next step
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMMreuse_workEstimation(h->cusparseHandle, h->opA, h->opB, h->descr_A,
                                                                 h->descr_B, h->descr_C, h->alg, h->spgemmDescr,
                                                                 &bufferSize1, buffer1));

    //----------------------------------------------------------------------
    // Compute nnz of C
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMMreuse_nnz(h->cusparseHandle, h->opA, h->opB, h->descr_A, h->descr_B,
                                                      h->descr_C, h->alg, h->spgemmDescr, &bufferSize2, nullptr,
                                                      &h->bufferSize3, nullptr, &h->bufferSize4, nullptr));

    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&buffer2, bufferSize2));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&h->buffer3, h->bufferSize3));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&h->buffer4, h->bufferSize4));

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMMreuse_nnz(h->cusparseHandle, h->opA, h->opB, h->descr_A, h->descr_B,
                                                      h->descr_C, h->alg, h->spgemmDescr, &bufferSize2, buffer2,
                                                      &h->bufferSize3, h->buffer3, &h->bufferSize4, h->buffer4));

    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(buffer2));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(buffer1));

    int64_t C_nrow, C_ncol, C_nnz;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMatGetSize(h->descr_C, &C_nrow, &C_ncol, &C_nnz));
    if (C_nnz > std::numeric_limits<int>::max()) {
      throw std::runtime_error("nnz of C overflowed over 32-bit int\n");
    }
    handle->set_c_nnz(C_nnz);
    handle->set_call_symbolic();
  }
  if (computeRowptrs && !handle->are_rowptrs_computed()) {
    using Scalar  = typename KernelHandle::nnz_scalar_t;
    using Ordinal = typename KernelHandle::nnz_lno_t;
    using Offset  = typename KernelHandle::size_type;
    Ordinal *dummyEntries;
    Scalar *dummyValues;
    auto C_nnz = handle->get_c_nnz();
    auto h     = handle->get_cusparse_spgemm_handle();
    // We just want rowptrs, but since C's entries/values are not yet allocated,
    // we must use dummy versions and then discard them.
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&dummyEntries, C_nnz * sizeof(Ordinal)));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&dummyValues, C_nnz * sizeof(Scalar)));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCsrSetPointers(h->descr_C, row_mapC.data(), dummyEntries, dummyValues));
    //--------------------------------------------------------------------------

    cusparseSpGEMMreuse_copy(h->cusparseHandle, h->opA, h->opB, h->descr_A, h->descr_B, h->descr_C, h->alg,
                             h->spgemmDescr, &h->bufferSize5, nullptr);
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&h->buffer5, h->bufferSize5));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMMreuse_copy(h->cusparseHandle, h->opA, h->opB, h->descr_A, h->descr_B,
                                                       h->descr_C, h->alg, h->spgemmDescr, &h->bufferSize5,
                                                       h->buffer5));
    if (!handle->get_c_nnz()) {
      // cuSPARSE does not populate C rowptrs if C has no entries
      Kokkos::deep_copy(typename KernelHandle::HandleExecSpace(), row_mapC, Offset(0));
    }
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(dummyValues));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(dummyEntries));
    handle->set_computed_rowptrs();
  }
}

#elif (CUDA_VERSION >= 11000)
// 11.0-11.3 supports only the generic API, but not reuse.
template <typename KernelHandle, typename lno_t, typename ConstRowMapType, typename ConstEntriesType,
          typename RowMapType>
void spgemm_symbolic_cusparse(KernelHandle *handle, lno_t m, lno_t n, lno_t k, const ConstRowMapType &row_mapA,
                              const ConstEntriesType &entriesA, const ConstRowMapType &row_mapB,
                              const ConstEntriesType &entriesB, const RowMapType &row_mapC, bool computeRowptrs) {
  using scalar_type      = typename KernelHandle::nnz_scalar_t;
  using ordinal_type     = typename KernelHandle::nnz_lno_t;
  const auto alpha       = Kokkos::ArithTraits<scalar_type>::one();
  const auto beta        = Kokkos::ArithTraits<scalar_type>::zero();
  void *dummyValues_AB   = nullptr;
  bool firstSymbolicCall = false;
  if (!handle->is_symbolic_called()) {
    handle->create_cusparse_spgemm_handle(false, false);
    auto h = handle->get_cusparse_spgemm_handle();

    // Follow
    // https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuSPARSE/spgemm

    // In non-reuse interface, forced to give A,B dummy values to
    // cusparseSpGEMM_compute. And it actually reads them, so they must be
    // allocated and of the correct type. This compute will be called again in
    // numeric with the real values.
    //
    // The dummy values can be uninitialized. cusparseSpGEMM_compute does
    // not remove numerical zeros from the sparsity pattern.
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        cudaMalloc(&dummyValues_AB, sizeof(scalar_type) * std::max(entriesA.extent(0), entriesB.extent(0))));

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&h->descr_A, m, n, entriesA.extent(0), (void *)row_mapA.data(),
                                                (void *)entriesA.data(), dummyValues_AB, CUSPARSE_INDEX_32I,
                                                CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, h->scalarType));

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&h->descr_B, n, k, entriesB.extent(0), (void *)row_mapB.data(),
                                                (void *)entriesB.data(), dummyValues_AB, CUSPARSE_INDEX_32I,
                                                CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, h->scalarType));

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&h->descr_C, m, k, 0, row_mapC.data(), nullptr, nullptr,
                                                CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                                                h->scalarType));

    //----------------------------------------------------------------------
    // query workEstimation buffer size, allocate, then call again with buffer.
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_workEstimation(h->cusparseHandle, h->opA, h->opB, &alpha, h->descr_A,
                                                            h->descr_B, &beta, h->descr_C, h->scalarType, h->alg,
                                                            h->spgemmDescr, &h->bufferSize3, nullptr));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&h->buffer3, h->bufferSize3));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_workEstimation(h->cusparseHandle, h->opA, h->opB, &alpha, h->descr_A,
                                                            h->descr_B, &beta, h->descr_C, h->scalarType, h->alg,
                                                            h->spgemmDescr, &h->bufferSize3, h->buffer3));

    //----------------------------------------------------------------------
    // query compute buffer size, allocate, then call again with buffer.

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_compute(h->cusparseHandle, h->opA, h->opB, &alpha, h->descr_A, h->descr_B,
                                                     &beta, h->descr_C, h->scalarType, CUSPARSE_SPGEMM_DEFAULT,
                                                     h->spgemmDescr, &h->bufferSize4, nullptr));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&h->buffer4, h->bufferSize4));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_compute(h->cusparseHandle, h->opA, h->opB, &alpha, h->descr_A, h->descr_B,
                                                     &beta, h->descr_C, h->scalarType, CUSPARSE_SPGEMM_DEFAULT,
                                                     h->spgemmDescr, &h->bufferSize4, h->buffer4));
    int64_t C_nrow, C_ncol, C_nnz;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMatGetSize(h->descr_C, &C_nrow, &C_ncol, &C_nnz));
    if (C_nnz > std::numeric_limits<int>::max()) {
      throw std::runtime_error("nnz of C overflowed over 32-bit int\n");
    }
    handle->set_c_nnz(C_nnz);
    handle->set_call_symbolic();
    firstSymbolicCall = true;
  }

  if (computeRowptrs && !handle->are_rowptrs_computed()) {
    auto h     = handle->get_cusparse_spgemm_handle();
    auto C_nnz = handle->get_c_nnz();
    if (!firstSymbolicCall) {
      // This is not the first call to symbolic, so dummyValues_AB was not
      // allocated above. But, descr_A and descr_B will have been saved in the
      // handle, so we can reuse those.
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          cudaMalloc(&dummyValues_AB, sizeof(scalar_type) * std::max(entriesA.extent(0), entriesB.extent(0))));
      KOKKOS_CUSPARSE_SAFE_CALL(
          cusparseCsrSetPointers(h->descr_A, (void *)row_mapA.data(), (void *)entriesA.data(), dummyValues_AB));
      KOKKOS_CUSPARSE_SAFE_CALL(
          cusparseCsrSetPointers(h->descr_B, (void *)row_mapB.data(), (void *)entriesB.data(), dummyValues_AB));
    }
    void *dummyEntries_C, *dummyValues_C;
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc(&dummyEntries_C, sizeof(ordinal_type) * C_nnz));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc(&dummyValues_C, sizeof(scalar_type) * C_nnz));
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCsrSetPointers(h->descr_C, (void *)row_mapC.data(), dummyEntries_C, dummyValues_C));

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_copy(h->cusparseHandle, h->opA, h->opB, &alpha, h->descr_A, h->descr_B,
                                                  &beta, h->descr_C, h->scalarType, CUSPARSE_SPGEMM_DEFAULT,
                                                  h->spgemmDescr));

    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(dummyValues_C));
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(dummyEntries_C));
    handle->set_computed_rowptrs();
  }
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(dummyValues_AB));
}

#else
// 10.x supports the pre-generic interface (cusparseXcsrgemmNnz). It always
// populates C rowptrs.
template <typename KernelHandle, typename lno_t, typename ConstRowMapType, typename ConstEntriesType,
          typename RowMapType>
void spgemm_symbolic_cusparse(KernelHandle *handle, lno_t m, lno_t n, lno_t k, const ConstRowMapType &row_mapA,
                              const ConstEntriesType &entriesA, const ConstRowMapType &row_mapB,
                              const ConstEntriesType &entriesB, const RowMapType &row_mapC, bool /* computeRowptrs */) {
  // using scalar_type = typename KernelHandle::nnz_scalar_t;
  using size_type = typename KernelHandle::size_type;
  if (handle->are_rowptrs_computed()) return;
  handle->create_cusparse_spgemm_handle(false, false);
  auto h   = handle->get_cusparse_spgemm_handle();
  int nnzA = entriesA.extent(0);
  int nnzB = entriesB.extent(0);

  int baseC, nnzC;
  int *nnzTotalDevHostPtr = &nnzC;

  // In empty (zero entries) matrix case, cusparse does not populate rowptrs to
  // zeros
  if (m == 0 || n == 0 || k == 0 || entriesA.extent(0) == size_type(0) || entriesB.extent(0) == size_type(0)) {
    Kokkos::deep_copy(typename KernelHandle::HandleExecSpace(), row_mapC, size_type(0));
    nnzC = 0;
  } else {
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseXcsrgemmNnz(h->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, m, k,
                            n, h->generalDescr, nnzA, row_mapA.data(), entriesA.data(), h->generalDescr, nnzB,
                            row_mapB.data(), entriesB.data(), h->generalDescr, row_mapC.data(), nnzTotalDevHostPtr));
    if (nullptr != nnzTotalDevHostPtr) {
      nnzC = *nnzTotalDevHostPtr;
    } else {
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMemcpy(&nnzC, row_mapC.data() + m, sizeof(int), cudaMemcpyDeviceToHost));
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMemcpy(&baseC, row_mapC.data(), sizeof(int), cudaMemcpyDeviceToHost));
      nnzC -= baseC;
    }
  }
  handle->set_c_nnz(nnzC);
  handle->set_call_symbolic();
  handle->set_computed_rowptrs();
}

#endif

#define SPGEMM_SYMBOLIC_DECL_CUSPARSE(SCALAR, MEMSPACE, TPL_AVAIL)                                                     \
  template <>                                                                                                          \
  struct SPGEMM_SYMBOLIC<KokkosKernels::Experimental::KokkosKernelsHandle<const int, const int, const SCALAR,          \
                                                                          Kokkos::Cuda, MEMSPACE, MEMSPACE>,           \
                         Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,             \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,             \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,             \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,             \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                   \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         true, TPL_AVAIL> {                                                                            \
    using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<const int, const int, const SCALAR,          \
                                                                          Kokkos::Cuda, MEMSPACE, MEMSPACE>;           \
    using c_int_view_t = Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,             \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    using int_view_t   = Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                   \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    static void spgemm_symbolic(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,                              \
                                typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,                \
                                c_int_view_t row_mapA, c_int_view_t entriesA, bool, c_int_view_t row_mapB,             \
                                c_int_view_t entriesB, bool, int_view_t row_mapC, bool computeRowptrs) {               \
      std::string label = "KokkosSparse::spgemm_symbolic[TPL_CUSPARSE," + Kokkos::ArithTraits<SCALAR>::name() + "]";   \
      Kokkos::Profiling::pushRegion(label);                                                                            \
      spgemm_symbolic_cusparse(handle->get_spgemm_handle(), m, n, k, row_mapA, entriesA, row_mapB, entriesB, row_mapC, \
                               computeRowptrs);                                                                        \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define SPGEMM_SYMBOLIC_DECL_CUSPARSE_S(SCALAR, TPL_AVAIL)            \
  SPGEMM_SYMBOLIC_DECL_CUSPARSE(SCALAR, Kokkos::CudaSpace, TPL_AVAIL) \
  SPGEMM_SYMBOLIC_DECL_CUSPARSE(SCALAR, Kokkos::CudaUVMSpace, TPL_AVAIL)

SPGEMM_SYMBOLIC_DECL_CUSPARSE_S(float, true)
SPGEMM_SYMBOLIC_DECL_CUSPARSE_S(double, true)
SPGEMM_SYMBOLIC_DECL_CUSPARSE_S(Kokkos::complex<float>, true)
SPGEMM_SYMBOLIC_DECL_CUSPARSE_S(Kokkos::complex<double>, true)

SPGEMM_SYMBOLIC_DECL_CUSPARSE_S(float, false)
SPGEMM_SYMBOLIC_DECL_CUSPARSE_S(double, false)
SPGEMM_SYMBOLIC_DECL_CUSPARSE_S(Kokkos::complex<float>, false)
SPGEMM_SYMBOLIC_DECL_CUSPARSE_S(Kokkos::complex<double>, false)

#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
//=============================================================================
// Overload rocsparse_Xcsrgemm_buffer_size() over scalar types
#define ROCSPARSE_XCSRGEMM_BUFFER_SIZE_SPEC(scalar_type, TOKEN)                                                 \
  inline rocsparse_status rocsparse_Xcsrgemm_buffer_size(                                                       \
      rocsparse_handle handle, rocsparse_operation trans_A, rocsparse_operation trans_B, rocsparse_int m,       \
      rocsparse_int n, rocsparse_int k, const scalar_type *alpha, const rocsparse_mat_descr descr_A,            \
      rocsparse_int nnz_A, const rocsparse_int *csr_row_ptr_A, const rocsparse_int *csr_col_ind_A,              \
      const rocsparse_mat_descr descr_B, rocsparse_int nnz_B, const rocsparse_int *csr_row_ptr_B,               \
      const rocsparse_int *csr_col_ind_B, const scalar_type *beta, const rocsparse_mat_descr descr_D,           \
      rocsparse_int nnz_D, const rocsparse_int *csr_row_ptr_D, const rocsparse_int *csr_col_ind_D,              \
      rocsparse_mat_info info_C, size_t *buffer_size) {                                                         \
    return rocsparse_##TOKEN##csrgemm_buffer_size(                                                              \
        handle, trans_A, trans_B, m, n, k, alpha, descr_A, nnz_A, csr_row_ptr_A, csr_col_ind_A, descr_B, nnz_B, \
        csr_row_ptr_B, csr_col_ind_B, beta, descr_D, nnz_D, csr_row_ptr_D, csr_col_ind_D, info_C, buffer_size); \
  }

ROCSPARSE_XCSRGEMM_BUFFER_SIZE_SPEC(float, s)
ROCSPARSE_XCSRGEMM_BUFFER_SIZE_SPEC(double, d)
ROCSPARSE_XCSRGEMM_BUFFER_SIZE_SPEC(rocsparse_float_complex, c)
ROCSPARSE_XCSRGEMM_BUFFER_SIZE_SPEC(rocsparse_double_complex, z)

template <typename KernelHandle, typename ain_row_index_view_type, typename ain_nonzero_index_view_type,
          typename bin_row_index_view_type, typename bin_nonzero_index_view_type, typename cin_row_index_view_type>
void spgemm_symbolic_rocsparse(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,
                               typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,
                               ain_row_index_view_type rowptrA, ain_nonzero_index_view_type colidxA,
                               bin_row_index_view_type rowptrB, bin_nonzero_index_view_type colidxB,
                               cin_row_index_view_type rowptrC) {
  using index_type            = typename KernelHandle::nnz_lno_t;
  using scalar_type           = typename KernelHandle::nnz_scalar_t;
  using rocsparse_scalar_type = typename kokkos_to_rocsparse_type<scalar_type>::type;

  auto nnz_A = colidxA.extent(0);
  auto nnz_B = colidxB.extent(0);

  if (handle->is_symbolic_called()) {
    return;
  }
  handle->create_rocsparse_spgemm_handle(false, false);
  typename KernelHandle::rocSparseSpgemmHandleType *h = handle->get_rocsparse_spgemm_handle();

  // alpha, beta are on host, but since we use singleton on the rocsparse
  // handle, we save/restore the pointer mode to not interference with
  // others' use
  const auto alpha = Kokkos::ArithTraits<scalar_type>::one();
  const auto beta  = Kokkos::ArithTraits<scalar_type>::zero();
  rocsparse_pointer_mode oldPtrMode;

  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_get_pointer_mode(h->rocsparseHandle, &oldPtrMode));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_set_pointer_mode(h->rocsparseHandle, rocsparse_pointer_mode_host));

  // C = alpha * OpA(A) * OpB(B) + beta * D
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_Xcsrgemm_buffer_size(
      h->rocsparseHandle, h->opA, h->opB, m, k, n, reinterpret_cast<const rocsparse_scalar_type *>(&alpha), h->descr_A,
      nnz_A, rowptrA.data(), colidxA.data(), h->descr_B, nnz_B, rowptrB.data(), colidxB.data(),
      reinterpret_cast<const rocsparse_scalar_type *>(&beta), h->descr_D, 0, nullptr, nullptr, h->info_C,
      &h->bufferSize));

  KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(&h->buffer, h->bufferSize));

  rocsparse_int nnz_C = 0;
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_csrgemm_nnz(h->rocsparseHandle, h->opA, h->opB, m, k, n, h->descr_A, nnz_A,
                                                        rowptrA.data(), colidxA.data(), h->descr_B, nnz_B,
                                                        rowptrB.data(), colidxB.data(), h->descr_D, 0, nullptr, nullptr,
                                                        h->descr_C, rowptrC.data(), &nnz_C, h->info_C, h->buffer));
  // If C has zero rows, its rowptrs are not populated
  if (m == 0) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipMemset(rowptrC.data(), 0, rowptrC.extent(0) * sizeof(index_type)));
  }
  // Restore previous pointer mode
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_set_pointer_mode(h->rocsparseHandle, oldPtrMode));

  handle->set_c_nnz(nnz_C);
  handle->set_call_symbolic();
  handle->set_computed_rowptrs();
}

#define SPGEMM_SYMBOLIC_DECL_ROCSPARSE(SCALAR, TPL_AVAIL)                                                             \
  template <>                                                                                                         \
  struct SPGEMM_SYMBOLIC<KokkosKernels::Experimental::KokkosKernelsHandle<                                            \
                             const int, const int, const SCALAR, Kokkos::HIP, Kokkos::HIPSpace, Kokkos::HIPSpace>,    \
                         Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,     \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                         Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,     \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                         Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,     \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                         Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,     \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                         Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,           \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                         true, TPL_AVAIL> {                                                                           \
    using KernelHandle =                                                                                              \
        KokkosKernels::Experimental::KokkosKernelsHandle<const int, const int, const SCALAR, Kokkos::HIP,             \
                                                         Kokkos::HIPSpace, Kokkos::HIPSpace>;                         \
    using c_int_view_t = Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,     \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                       \
    using int_view_t   = Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,           \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                       \
    static void spgemm_symbolic(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,                             \
                                typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,               \
                                c_int_view_t row_mapA, c_int_view_t entriesA, bool, c_int_view_t row_mapB,            \
                                c_int_view_t entriesB, bool, int_view_t row_mapC, bool) {                             \
      std::string label = "KokkosSparse::spgemm_symbolic[TPL_ROCSPARSE," + Kokkos::ArithTraits<SCALAR>::name() + "]"; \
      Kokkos::Profiling::pushRegion(label);                                                                           \
      spgemm_symbolic_rocsparse(handle->get_spgemm_handle(), m, n, k, row_mapA, entriesA, row_mapB, entriesB,         \
                                row_mapC);                                                                            \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

SPGEMM_SYMBOLIC_DECL_ROCSPARSE(float, false)
SPGEMM_SYMBOLIC_DECL_ROCSPARSE(double, false)
SPGEMM_SYMBOLIC_DECL_ROCSPARSE(Kokkos::complex<float>, false)
SPGEMM_SYMBOLIC_DECL_ROCSPARSE(Kokkos::complex<double>, false)

SPGEMM_SYMBOLIC_DECL_ROCSPARSE(float, true)
SPGEMM_SYMBOLIC_DECL_ROCSPARSE(double, true)
SPGEMM_SYMBOLIC_DECL_ROCSPARSE(Kokkos::complex<float>, true)
SPGEMM_SYMBOLIC_DECL_ROCSPARSE(Kokkos::complex<double>, true)
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
template <typename KernelHandle, typename ain_row_index_view_type, typename ain_nonzero_index_view_type,
          typename bin_row_index_view_type, typename bin_nonzero_index_view_type, typename cin_row_index_view_type>
void spgemm_symbolic_mkl(KernelHandle *handle, typename KernelHandle::nnz_lno_t m, typename KernelHandle::nnz_lno_t n,
                         typename KernelHandle::nnz_lno_t k, ain_row_index_view_type rowptrA,
                         ain_nonzero_index_view_type colidxA, bin_row_index_view_type rowptrB,
                         bin_nonzero_index_view_type colidxB, cin_row_index_view_type rowptrC) {
  using ExecSpace   = typename KernelHandle::HandleExecSpace;
  using index_type  = typename KernelHandle::nnz_lno_t;
  using size_type   = typename KernelHandle::size_type;
  using scalar_type = typename KernelHandle::nnz_scalar_t;
  using MKLMatrix   = MKLSparseMatrix<scalar_type>;
  if (m == 0 || n == 0 || k == 0 || colidxA.extent(0) == size_type(0) || colidxB.extent(0) == size_type(0)) {
    Kokkos::deep_copy(ExecSpace(), rowptrC, size_type(0));
    handle->set_call_symbolic();
    handle->set_computed_rowptrs();
    handle->set_c_nnz(0);
    return;
  }
  MKLMatrix A(m, n, (MKL_INT *)rowptrA.data(), (MKL_INT *)colidxA.data(), nullptr);
  MKLMatrix B(n, k, (MKL_INT *)rowptrB.data(), (MKL_INT *)colidxB.data(), nullptr);
  sparse_matrix_t C;
  matrix_descr generalDescr;
  generalDescr.type = SPARSE_MATRIX_TYPE_GENERAL;
  generalDescr.mode = SPARSE_FILL_MODE_FULL;
  generalDescr.diag = SPARSE_DIAG_NON_UNIT;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, generalDescr, A,
                                              SPARSE_OPERATION_NON_TRANSPOSE, generalDescr, B, SPARSE_STAGE_NNZ_COUNT,
                                              &C));
  MKLMatrix wrappedC(C);
  MKL_INT nrows = 0, ncols = 0;
  MKL_INT *rowptrRaw     = nullptr;
  MKL_INT *colidxRaw     = nullptr;
  scalar_type *valuesRaw = nullptr;
  wrappedC.export_data(nrows, ncols, rowptrRaw, colidxRaw, valuesRaw);
  Kokkos::View<index_type *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> rowptrRawView(rowptrRaw,
                                                                                                       nrows + 1);
  Kokkos::deep_copy(ExecSpace(), rowptrC, rowptrRawView);
  handle->create_mkl_spgemm_handle(C);
  handle->set_call_symbolic();
  handle->set_computed_rowptrs();
  handle->set_c_nnz(rowptrC(m));
}

#define SPGEMM_SYMBOLIC_DECL_MKL(SCALAR, EXEC, TPL_AVAIL)                                                              \
  template <>                                                                                                          \
  struct SPGEMM_SYMBOLIC<KokkosKernels::Experimental::KokkosKernelsHandle<const MKL_INT, const MKL_INT, const SCALAR,  \
                                                                          EXEC, Kokkos::HostSpace, Kokkos::HostSpace>, \
                         Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,        \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,        \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,        \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,        \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         Kokkos::View<MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,              \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                         true, TPL_AVAIL> {                                                                            \
    using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<const MKL_INT, const MKL_INT, const SCALAR,  \
                                                                          EXEC, Kokkos::HostSpace, Kokkos::HostSpace>; \
    using c_int_view_t = Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,        \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    using int_view_t   = Kokkos::View<MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,              \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    static void spgemm_symbolic(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,                              \
                                typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,                \
                                c_int_view_t row_mapA, c_int_view_t entriesA, bool, c_int_view_t row_mapB,             \
                                c_int_view_t entriesB, bool, int_view_t row_mapC, bool) {                              \
      std::string label = "KokkosSparse::spgemm_symbolic[TPL_MKL," + Kokkos::ArithTraits<SCALAR>::name() + "]";        \
      Kokkos::Profiling::pushRegion(label);                                                                            \
      spgemm_symbolic_mkl(handle->get_spgemm_handle(), m, n, k, row_mapA, entriesA, row_mapB, entriesB, row_mapC);     \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define SPGEMM_SYMBOLIC_DECL_MKL_SE(SCALAR, EXEC) \
  SPGEMM_SYMBOLIC_DECL_MKL(SCALAR, EXEC, true)    \
  SPGEMM_SYMBOLIC_DECL_MKL(SCALAR, EXEC, false)

#define SPGEMM_SYMBOLIC_DECL_MKL_E(EXEC)                    \
  SPGEMM_SYMBOLIC_DECL_MKL_SE(float, EXEC)                  \
  SPGEMM_SYMBOLIC_DECL_MKL_SE(double, EXEC)                 \
  SPGEMM_SYMBOLIC_DECL_MKL_SE(Kokkos::complex<float>, EXEC) \
  SPGEMM_SYMBOLIC_DECL_MKL_SE(Kokkos::complex<double>, EXEC)

#ifdef KOKKOS_ENABLE_SERIAL
SPGEMM_SYMBOLIC_DECL_MKL_E(Kokkos::Serial)
#endif
#ifdef KOKKOS_ENABLE_OPENMP
SPGEMM_SYMBOLIC_DECL_MKL_E(Kokkos::OpenMP)
#endif
#endif

}  // namespace Impl
}  // namespace KokkosSparse

#endif
