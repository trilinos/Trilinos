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

#ifndef KOKKOSPARSE_SPGEMM_NUMERIC_TPL_SPEC_DECL_HPP_
#define KOKKOSPARSE_SPGEMM_NUMERIC_TPL_SPEC_DECL_HPP_

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
#if (CUDA_VERSION >= 11040)

// 11.4+ supports generic API with reuse (full symbolic/numeric separation)
template <typename KernelHandle, typename lno_t, typename ConstRowMapType, typename ConstEntriesType,
          typename ConstValuesType, typename EntriesType, typename ValuesType>
void spgemm_numeric_cusparse(KernelHandle *handle, lno_t /*m*/, lno_t /*n*/, lno_t /*k*/,
                             const ConstRowMapType &row_mapA, const ConstEntriesType &entriesA,
                             const ConstValuesType &valuesA, const ConstRowMapType &row_mapB,
                             const ConstEntriesType &entriesB, const ConstValuesType &valuesB,
                             const ConstRowMapType &row_mapC, const EntriesType &entriesC, const ValuesType &valuesC) {
  using scalar_type = typename KernelHandle::nnz_scalar_t;
  using size_type   = typename KernelHandle::size_type;
  auto h            = handle->get_cusparse_spgemm_handle();

  if (handle->get_c_nnz() == size_type(0)) {
    // Handle empty C case. entriesC and valuesC have extent 0 so nothing needs
    // to be done for them, but we must populate row_mapC to zeros if not
    // already done.
    if (!handle->are_rowptrs_computed()) {
      Kokkos::View<size_type *, typename ConstRowMapType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          row_mapC_nonconst(const_cast<size_type *>(row_mapC.data()), row_mapC.extent(0));
      Kokkos::deep_copy(typename KernelHandle::HandleExecSpace(), row_mapC_nonconst, size_type(0));
      handle->set_computed_rowptrs();
    }
    handle->set_computed_entries();
    handle->set_call_numeric();
    return;
  }

  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseCsrSetPointers(h->descr_A, (void *)row_mapA.data(), (void *)entriesA.data(), (void *)valuesA.data()));

  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseCsrSetPointers(h->descr_B, (void *)row_mapB.data(), (void *)entriesB.data(), (void *)valuesB.data()));

  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseCsrSetPointers(h->descr_C, (void *)row_mapC.data(), (void *)entriesC.data(), (void *)valuesC.data()));

  if (!handle->are_entries_computed()) {
    if (!h->buffer5) {
      // If symbolic was previously called with computeRowptrs=true, then
      // buffer5 will have already been allocated to the correct size. Otherwise
      // size and allocate it here.
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMMreuse_copy(h->cusparseHandle, h->opA, h->opB, h->descr_A, h->descr_B,
                                                         h->descr_C, h->alg, h->spgemmDescr, &h->bufferSize5, nullptr));
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&h->buffer5, h->bufferSize5));
    }
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMMreuse_copy(h->cusparseHandle, h->opA, h->opB, h->descr_A, h->descr_B,
                                                       h->descr_C, h->alg, h->spgemmDescr, &h->bufferSize5,
                                                       h->buffer5));
    handle->set_computed_rowptrs();
    handle->set_computed_entries();
  }

  // C' = alpha * opA(A) * opB(B) + beta * C
  const auto alpha = Kokkos::ArithTraits<scalar_type>::one();
  const auto beta  = Kokkos::ArithTraits<scalar_type>::zero();

  // alpha, beta are on host, but since we use singleton on the cusparse
  // handle, we save/restore the pointer mode to not interference with
  // others' use
  cusparsePointerMode_t oldPtrMode;
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseGetPointerMode(h->cusparseHandle, &oldPtrMode));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetPointerMode(h->cusparseHandle, CUSPARSE_POINTER_MODE_HOST));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMMreuse_compute(h->cusparseHandle, h->opA, h->opB, &alpha, h->descr_A,
                                                        h->descr_B, &beta, h->descr_C, h->scalarType, h->alg,
                                                        h->spgemmDescr));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetPointerMode(h->cusparseHandle, oldPtrMode));
  handle->set_call_numeric();
}

#elif (CUDA_VERSION >= 11000)
// 11.0-11.3 supports only the generic API, but not reuse.
template <typename KernelHandle, typename lno_t, typename ConstRowMapType, typename ConstEntriesType,
          typename ConstValuesType, typename EntriesType, typename ValuesType>
void spgemm_numeric_cusparse(KernelHandle *handle, lno_t /*m*/, lno_t /*n*/, lno_t /*k*/,
                             const ConstRowMapType &row_mapA, const ConstEntriesType &entriesA,
                             const ConstValuesType &valuesA, const ConstRowMapType &row_mapB,
                             const ConstEntriesType &entriesB, const ConstValuesType &valuesB,
                             const ConstRowMapType &row_mapC, const EntriesType &entriesC, const ValuesType &valuesC) {
  using scalar_type = typename KernelHandle::nnz_scalar_t;
  auto h            = handle->get_cusparse_spgemm_handle();
  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseCsrSetPointers(h->descr_A, (void *)row_mapA.data(), (void *)entriesA.data(), (void *)valuesA.data()));
  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseCsrSetPointers(h->descr_B, (void *)row_mapB.data(), (void *)entriesB.data(), (void *)valuesB.data()));
  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseCsrSetPointers(h->descr_C, (void *)row_mapC.data(), (void *)entriesC.data(), (void *)valuesC.data()));
  const auto alpha = Kokkos::ArithTraits<scalar_type>::one();
  const auto beta  = Kokkos::ArithTraits<scalar_type>::zero();
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_compute(h->cusparseHandle, h->opA, h->opB, &alpha, h->descr_A, h->descr_B,
                                                   &beta, h->descr_C, h->scalarType, CUSPARSE_SPGEMM_DEFAULT,
                                                   h->spgemmDescr, &h->bufferSize4, h->buffer4));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_copy(h->cusparseHandle, h->opA, h->opB, &alpha, h->descr_A, h->descr_B,
                                                &beta, h->descr_C, h->scalarType, CUSPARSE_SPGEMM_DEFAULT,
                                                h->spgemmDescr));
  handle->set_computed_entries();
  handle->set_call_numeric();
}

#else

// Generic (using overloads) wrapper for cusparseXcsrgemm (where X is S, D, C,
// or Z). Accepts Kokkos types (e.g. Kokkos::complex<float>) for Scalar and
// handles casting to cuSparse types internally.

#define CUSPARSE_XCSRGEMM_SPEC(KokkosType, CusparseType, Abbreviation)                                                 \
  inline cusparseStatus_t cusparseXcsrgemm(                                                                            \
      cusparseHandle_t handle, cusparseOperation_t transA, cusparseOperation_t transB, int m, int n, int k,            \
      const cusparseMatDescr_t descrA, const int nnzA, const KokkosType *csrSortedValA, const int *csrSortedRowPtrA,   \
      const int *csrSortedColIndA, const cusparseMatDescr_t descrB, const int nnzB, const KokkosType *csrSortedValB,   \
      const int *csrSortedRowPtrB, const int *csrSortedColIndB, const cusparseMatDescr_t descrC,                       \
      KokkosType *csrSortedValC, const int *csrSortedRowPtrC, int *csrSortedColIndC) {                                 \
    return cusparse##Abbreviation##csrgemm(                                                                            \
        handle, transA, transB, m, n, k, descrA, nnzA, reinterpret_cast<const CusparseType *>(csrSortedValA),          \
        csrSortedRowPtrA, csrSortedColIndA, descrB, nnzB, reinterpret_cast<const CusparseType *>(csrSortedValB),       \
        csrSortedRowPtrB, csrSortedColIndB, descrC, reinterpret_cast<CusparseType *>(csrSortedValC), csrSortedRowPtrC, \
        csrSortedColIndC);                                                                                             \
  }

CUSPARSE_XCSRGEMM_SPEC(float, float, S)
CUSPARSE_XCSRGEMM_SPEC(double, double, D)
CUSPARSE_XCSRGEMM_SPEC(Kokkos::complex<float>, cuComplex, C)
CUSPARSE_XCSRGEMM_SPEC(Kokkos::complex<double>, cuDoubleComplex, Z)

#undef CUSPARSE_XCSRGEMM_SPEC

// 10.x supports the pre-generic interface.
template <typename KernelHandle, typename lno_t, typename ConstRowMapType, typename ConstEntriesType,
          typename ConstValuesType, typename EntriesType, typename ValuesType>
void spgemm_numeric_cusparse(KernelHandle *handle, lno_t m, lno_t n, lno_t k, const ConstRowMapType &row_mapA,
                             const ConstEntriesType &entriesA, const ConstValuesType &valuesA,
                             const ConstRowMapType &row_mapB, const ConstEntriesType &entriesB,
                             const ConstValuesType &valuesB, const ConstRowMapType &row_mapC,
                             const EntriesType &entriesC, const ValuesType &valuesC) {
  auto h = handle->get_cusparse_spgemm_handle();

  int nnzA = entriesA.extent(0);
  int nnzB = entriesB.extent(0);

  // Only call numeric if C actually has entries
  if (handle->get_c_nnz()) {
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseXcsrgemm(
        h->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, m, k, n, h->generalDescr,
        nnzA, valuesA.data(), row_mapA.data(), entriesA.data(), h->generalDescr, nnzB, valuesB.data(), row_mapB.data(),
        entriesB.data(), h->generalDescr, valuesC.data(), row_mapC.data(), entriesC.data()));
  }
  handle->set_computed_entries();
  handle->set_call_numeric();
}

#endif

#define SPGEMM_NUMERIC_DECL_CUSPARSE(SCALAR, MEMSPACE, TPL_AVAIL)                                                    \
  template <>                                                                                                        \
  struct SPGEMM_NUMERIC<KokkosKernels::Experimental::KokkosKernelsHandle<const int, const int, const SCALAR,         \
                                                                         Kokkos::Cuda, MEMSPACE, MEMSPACE>,          \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                        Kokkos::View<const SCALAR *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,         \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                        Kokkos::View<const SCALAR *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,         \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                        Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,                  \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                        Kokkos::View<SCALAR *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,               \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                       \
                        true, TPL_AVAIL> {                                                                           \
    using KernelHandle    = KokkosKernels::Experimental::KokkosKernelsHandle<const int, const int, const SCALAR,     \
                                                                          Kokkos::Cuda, MEMSPACE, MEMSPACE>;      \
    using c_int_view_t    = Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,        \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                   \
    using int_view_t      = Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,              \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                   \
    using c_scalar_view_t = Kokkos::View<const SCALAR *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,     \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                   \
    using scalar_view_t   = Kokkos::View<SCALAR *, default_layout, Kokkos::Device<Kokkos::Cuda, MEMSPACE>,           \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                   \
    static void spgemm_numeric(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,                             \
                               typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,               \
                               c_int_view_t row_mapA, c_int_view_t entriesA, c_scalar_view_t valuesA, bool,          \
                               c_int_view_t row_mapB, c_int_view_t entriesB, c_scalar_view_t valuesB, bool,          \
                               c_int_view_t row_mapC, int_view_t entriesC, scalar_view_t valuesC) {                  \
      std::string label = "KokkosSparse::spgemm_numeric[TPL_CUSPARSE," + Kokkos::ArithTraits<SCALAR>::name() + "]";  \
      Kokkos::Profiling::pushRegion(label);                                                                          \
      spgemm_numeric_cusparse(handle->get_spgemm_handle(), m, n, k, row_mapA, entriesA, valuesA, row_mapB, entriesB, \
                              valuesB, row_mapC, entriesC, valuesC);                                                 \
      Kokkos::Profiling::popRegion();                                                                                \
    }                                                                                                                \
  };

#define SPGEMM_NUMERIC_DECL_CUSPARSE_S(SCALAR, TPL_AVAIL)            \
  SPGEMM_NUMERIC_DECL_CUSPARSE(SCALAR, Kokkos::CudaSpace, TPL_AVAIL) \
  SPGEMM_NUMERIC_DECL_CUSPARSE(SCALAR, Kokkos::CudaUVMSpace, TPL_AVAIL)

SPGEMM_NUMERIC_DECL_CUSPARSE_S(float, true)
SPGEMM_NUMERIC_DECL_CUSPARSE_S(double, true)
SPGEMM_NUMERIC_DECL_CUSPARSE_S(Kokkos::complex<float>, true)
SPGEMM_NUMERIC_DECL_CUSPARSE_S(Kokkos::complex<double>, true)

SPGEMM_NUMERIC_DECL_CUSPARSE_S(float, false)
SPGEMM_NUMERIC_DECL_CUSPARSE_S(double, false)
SPGEMM_NUMERIC_DECL_CUSPARSE_S(Kokkos::complex<float>, false)
SPGEMM_NUMERIC_DECL_CUSPARSE_S(Kokkos::complex<double>, false)
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
//=============================================================================
// Overload rocsparse_Xcsrgemm_numeric() over scalar types
#define ROCSPARSE_XCSRGEMM_NUMERIC_SPEC(scalar_type, TOKEN)                                                            \
  inline rocsparse_status rocsparse_Xcsrgemm_numeric(                                                                  \
      rocsparse_handle handle, rocsparse_operation trans_A, rocsparse_operation trans_B, rocsparse_int m,              \
      rocsparse_int n, rocsparse_int k, const scalar_type *alpha, const rocsparse_mat_descr descr_A,                   \
      rocsparse_int nnz_A, const scalar_type *csr_val_A, const rocsparse_int *csr_row_ptr_A,                           \
      const rocsparse_int *csr_col_ind_A, const rocsparse_mat_descr descr_B, rocsparse_int nnz_B,                      \
      const scalar_type *csr_val_B, const rocsparse_int *csr_row_ptr_B, const rocsparse_int *csr_col_ind_B,            \
      const scalar_type *beta, const rocsparse_mat_descr descr_D, rocsparse_int nnz_D, const scalar_type *csr_val_D,   \
      const rocsparse_int *csr_row_ptr_D, const rocsparse_int *csr_col_ind_D, const rocsparse_mat_descr descr_C,       \
      rocsparse_int nnz_C, scalar_type *csr_val_C, const rocsparse_int *csr_row_ptr_C,                                 \
      const rocsparse_int *csr_col_ind_C, const rocsparse_mat_info info_C, void *buffer) {                             \
    return rocsparse_##TOKEN##csrgemm_numeric(                                                                         \
        handle, trans_A, trans_B, m, n, k, alpha, descr_A, nnz_A, csr_val_A, csr_row_ptr_A, csr_col_ind_A, descr_B,    \
        nnz_B, csr_val_B, csr_row_ptr_B, csr_col_ind_B, beta, descr_D, nnz_D, csr_val_D, csr_row_ptr_D, csr_col_ind_D, \
        descr_C, nnz_C, csr_val_C, csr_row_ptr_C, csr_col_ind_C, info_C, buffer);                                      \
  }

ROCSPARSE_XCSRGEMM_NUMERIC_SPEC(float, s)
ROCSPARSE_XCSRGEMM_NUMERIC_SPEC(double, d)
ROCSPARSE_XCSRGEMM_NUMERIC_SPEC(rocsparse_float_complex, c)
ROCSPARSE_XCSRGEMM_NUMERIC_SPEC(rocsparse_double_complex, z)

template <typename KernelHandle, typename ain_row_index_view_type, typename ain_nonzero_index_view_type,
          typename ain_nonzero_value_view_type, typename bin_row_index_view_type, typename bin_nonzero_index_view_type,
          typename bin_nonzero_value_view_type, typename cin_row_index_view_type, typename cin_nonzero_index_view_type,
          typename cin_nonzero_value_view_type>
void spgemm_numeric_rocsparse(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,
                              typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,
                              ain_row_index_view_type rowptrA, ain_nonzero_index_view_type colidxA,
                              ain_nonzero_value_view_type valuesA, bin_row_index_view_type rowptrB,
                              bin_nonzero_index_view_type colidxB, bin_nonzero_value_view_type valuesB,
                              cin_row_index_view_type rowptrC, cin_nonzero_index_view_type colidxC,
                              cin_nonzero_value_view_type valuesC) {
  using scalar_type           = typename KernelHandle::nnz_scalar_t;
  using rocsparse_scalar_type = typename kokkos_to_rocsparse_type<scalar_type>::type;

  typename KernelHandle::rocSparseSpgemmHandleType *h = handle->get_rocsparse_spgemm_handle();

  const auto alpha = Kokkos::ArithTraits<scalar_type>::one();
  const auto beta  = Kokkos::ArithTraits<scalar_type>::zero();
  rocsparse_pointer_mode oldPtrMode;

  auto nnz_A = colidxA.extent(0);
  auto nnz_B = colidxB.extent(0);
  auto nnz_C = colidxC.extent(0);

  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_get_pointer_mode(h->rocsparseHandle, &oldPtrMode));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_set_pointer_mode(h->rocsparseHandle, rocsparse_pointer_mode_host));

  if (!handle->are_entries_computed()) {
    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_csrgemm_symbolic(
        h->rocsparseHandle, h->opA, h->opB, m, k, n, h->descr_A, nnz_A, rowptrA.data(), colidxA.data(), h->descr_B,
        nnz_B, rowptrB.data(), colidxB.data(), h->descr_D, 0, nullptr, nullptr, h->descr_C, nnz_C, rowptrC.data(),
        colidxC.data(), h->info_C, h->buffer));
    handle->set_computed_entries();
  }

  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_Xcsrgemm_numeric(
      h->rocsparseHandle, h->opA, h->opB, m, k, n, reinterpret_cast<const rocsparse_scalar_type *>(&alpha), h->descr_A,
      nnz_A, reinterpret_cast<const rocsparse_scalar_type *>(valuesA.data()), rowptrA.data(), colidxA.data(),
      h->descr_B, nnz_B, reinterpret_cast<const rocsparse_scalar_type *>(valuesB.data()), rowptrB.data(),
      colidxB.data(), reinterpret_cast<const rocsparse_scalar_type *>(&beta), h->descr_D, 0, nullptr, nullptr, nullptr,
      h->descr_C, nnz_C, reinterpret_cast<rocsparse_scalar_type *>(valuesC.data()), rowptrC.data(), colidxC.data(),
      h->info_C, h->buffer));
  // Restore old pointer mode
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_set_pointer_mode(h->rocsparseHandle, oldPtrMode));
  handle->set_call_numeric();
}

#define SPGEMM_NUMERIC_DECL_ROCSPARSE(SCALAR, TPL_AVAIL)                                                              \
  template <>                                                                                                         \
  struct SPGEMM_NUMERIC<KokkosKernels::Experimental::KokkosKernelsHandle<                                             \
                            const int, const int, const SCALAR, Kokkos::HIP, Kokkos::HIPSpace, Kokkos::HIPSpace>,     \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,      \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,      \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                        Kokkos::View<const SCALAR *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,   \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,      \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,      \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                        Kokkos::View<const SCALAR *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,   \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                        Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,      \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                        Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                        Kokkos::View<SCALAR *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,         \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                        \
                        true, TPL_AVAIL> {                                                                            \
    using KernelHandle =                                                                                              \
        KokkosKernels::Experimental::KokkosKernelsHandle<const int, const int, const SCALAR, Kokkos::HIP,             \
                                                         Kokkos::HIPSpace, Kokkos::HIPSpace>;                         \
    using c_int_view_t = Kokkos::View<const int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,     \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                       \
    using int_view_t   = Kokkos::View<int *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,           \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                       \
    using c_scalar_view_t =                                                                                           \
        Kokkos::View<const SCALAR *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,                   \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                        \
    using scalar_view_t = Kokkos::View<SCALAR *, default_layout, Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,       \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                      \
    static void spgemm_numeric(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,                              \
                               typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,                \
                               c_int_view_t row_mapA, c_int_view_t entriesA, c_scalar_view_t valuesA, bool,           \
                               c_int_view_t row_mapB, c_int_view_t entriesB, c_scalar_view_t valuesB, bool,           \
                               c_int_view_t row_mapC, int_view_t entriesC, scalar_view_t valuesC) {                   \
      std::string label = "KokkosSparse::spgemm_numeric[TPL_ROCSPARSE," + Kokkos::ArithTraits<SCALAR>::name() + "]";  \
      Kokkos::Profiling::pushRegion(label);                                                                           \
      spgemm_numeric_rocsparse(handle->get_spgemm_handle(), m, n, k, row_mapA, entriesA, valuesA, row_mapB, entriesB, \
                               valuesB, row_mapC, entriesC, valuesC);                                                 \
      Kokkos::Profiling::popRegion();                                                                                 \
    }                                                                                                                 \
  };

SPGEMM_NUMERIC_DECL_ROCSPARSE(float, true)
SPGEMM_NUMERIC_DECL_ROCSPARSE(double, true)
SPGEMM_NUMERIC_DECL_ROCSPARSE(Kokkos::complex<float>, true)
SPGEMM_NUMERIC_DECL_ROCSPARSE(Kokkos::complex<double>, true)

SPGEMM_NUMERIC_DECL_ROCSPARSE(float, false)
SPGEMM_NUMERIC_DECL_ROCSPARSE(double, false)
SPGEMM_NUMERIC_DECL_ROCSPARSE(Kokkos::complex<float>, false)
SPGEMM_NUMERIC_DECL_ROCSPARSE(Kokkos::complex<double>, false)
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
template <typename KernelHandle, typename ain_row_index_view_type, typename ain_nonzero_index_view_type,
          typename ain_nonzero_value_view_type, typename bin_row_index_view_type, typename bin_nonzero_index_view_type,
          typename bin_nonzero_value_view_type, typename cin_row_index_view_type, typename cin_nonzero_index_view_type,
          typename cin_nonzero_value_view_type>
void spgemm_numeric_mkl(KernelHandle *handle, typename KernelHandle::nnz_lno_t m, typename KernelHandle::nnz_lno_t n,
                        typename KernelHandle::nnz_lno_t k, ain_row_index_view_type rowptrA,
                        ain_nonzero_index_view_type colidxA, ain_nonzero_value_view_type valuesA,
                        bin_row_index_view_type rowptrB, bin_nonzero_index_view_type colidxB,
                        bin_nonzero_value_view_type valuesB, cin_row_index_view_type rowptrC,
                        cin_nonzero_index_view_type colidxC, cin_nonzero_value_view_type valuesC) {
  using ExecSpace   = typename KernelHandle::HandleExecSpace;
  using index_type  = typename KernelHandle::nnz_lno_t;
  using size_type   = typename KernelHandle::size_type;
  using scalar_type = typename KernelHandle::nnz_scalar_t;
  using MKLMatrix   = MKLSparseMatrix<scalar_type>;
  size_type c_nnz   = handle->get_c_nnz();
  if (c_nnz == size_type(0)) {
    handle->set_computed_entries();
    handle->set_call_numeric();
    return;
  }
  MKLMatrix A(m, n, const_cast<size_type *>(rowptrA.data()), const_cast<index_type *>(colidxA.data()),
              const_cast<scalar_type *>(valuesA.data()));
  MKLMatrix B(n, k, const_cast<size_type *>(rowptrB.data()), const_cast<index_type *>(colidxB.data()),
              const_cast<scalar_type *>(valuesB.data()));
  auto mklSpgemmHandle = handle->get_mkl_spgemm_handle();
  matrix_descr generalDescr;
  generalDescr.type = SPARSE_MATRIX_TYPE_GENERAL;
  generalDescr.mode = SPARSE_FILL_MODE_FULL;
  generalDescr.diag = SPARSE_DIAG_NON_UNIT;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, generalDescr, A,
                                              SPARSE_OPERATION_NON_TRANSPOSE, generalDescr, B,
                                              SPARSE_STAGE_FINALIZE_MULT_NO_VAL, &mklSpgemmHandle->C));
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, generalDescr, A,
                                              SPARSE_OPERATION_NON_TRANSPOSE, generalDescr, B,
                                              SPARSE_STAGE_FINALIZE_MULT, &mklSpgemmHandle->C));
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_order(mklSpgemmHandle->C));
  MKLMatrix wrappedC(mklSpgemmHandle->C);
  MKL_INT nrows = 0, ncols = 0;
  MKL_INT *rowptrRaw     = nullptr;
  MKL_INT *colidxRaw     = nullptr;
  scalar_type *valuesRaw = nullptr;
  wrappedC.export_data(nrows, ncols, rowptrRaw, colidxRaw, valuesRaw);
  Kokkos::View<index_type *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> colidxRawView(colidxRaw,
                                                                                                       c_nnz);
  Kokkos::View<scalar_type *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> valuesRawView(valuesRaw,
                                                                                                        c_nnz);
  Kokkos::deep_copy(ExecSpace(), colidxC, colidxRawView);
  Kokkos::deep_copy(ExecSpace(), valuesC, valuesRawView);
  handle->set_call_numeric();
  handle->set_computed_entries();
}

#define SPGEMM_NUMERIC_DECL_MKL(SCALAR, EXEC, TPL_AVAIL)                                                                  \
  template <>                                                                                                             \
  struct SPGEMM_NUMERIC<KokkosKernels::Experimental::KokkosKernelsHandle<const MKL_INT, const MKL_INT, const SCALAR,      \
                                                                         EXEC, Kokkos::HostSpace, Kokkos::HostSpace>,     \
                        Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
                        Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
                        Kokkos::View<const SCALAR *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,             \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
                        Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
                        Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
                        Kokkos::View<const SCALAR *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,             \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
                        Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,            \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
                        Kokkos::View<MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,                  \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
                        Kokkos::View<SCALAR *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,                   \
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                            \
                        true, TPL_AVAIL> {                                                                                \
    using KernelHandle    = KokkosKernels::Experimental::KokkosKernelsHandle<const MKL_INT, const MKL_INT, const SCALAR,  \
                                                                          EXEC, Kokkos::HostSpace, Kokkos::HostSpace>; \
    using c_int_view_t    = Kokkos::View<const MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,        \
                                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    using int_view_t      = Kokkos::View<MKL_INT *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,              \
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    using c_scalar_view_t = Kokkos::View<const SCALAR *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,         \
                                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    using scalar_view_t   = Kokkos::View<SCALAR *, default_layout, Kokkos::Device<EXEC, Kokkos::HostSpace>,               \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                        \
    static void spgemm_numeric(KernelHandle *handle, typename KernelHandle::nnz_lno_t m,                                  \
                               typename KernelHandle::nnz_lno_t n, typename KernelHandle::nnz_lno_t k,                    \
                               c_int_view_t row_mapA, c_int_view_t entriesA, c_scalar_view_t valuesA, bool,               \
                               c_int_view_t row_mapB, c_int_view_t entriesB, c_scalar_view_t valuesB, bool,               \
                               c_int_view_t row_mapC, int_view_t entriesC, scalar_view_t valuesC) {                       \
      std::string label = "KokkosSparse::spgemm_numeric[TPL_MKL," + Kokkos::ArithTraits<SCALAR>::name() + "]";            \
      Kokkos::Profiling::pushRegion(label);                                                                               \
      spgemm_numeric_mkl(handle->get_spgemm_handle(), m, n, k, row_mapA, entriesA, valuesA, row_mapB, entriesB,           \
                         valuesB, row_mapC, entriesC, valuesC);                                                           \
      Kokkos::Profiling::popRegion();                                                                                     \
    }                                                                                                                     \
  };

#define SPGEMM_NUMERIC_DECL_MKL_SE(SCALAR, EXEC) \
  SPGEMM_NUMERIC_DECL_MKL(SCALAR, EXEC, true)    \
  SPGEMM_NUMERIC_DECL_MKL(SCALAR, EXEC, false)

#define SPGEMM_NUMERIC_DECL_MKL_E(EXEC)                    \
  SPGEMM_NUMERIC_DECL_MKL_SE(float, EXEC)                  \
  SPGEMM_NUMERIC_DECL_MKL_SE(double, EXEC)                 \
  SPGEMM_NUMERIC_DECL_MKL_SE(Kokkos::complex<float>, EXEC) \
  SPGEMM_NUMERIC_DECL_MKL_SE(Kokkos::complex<double>, EXEC)

#ifdef KOKKOS_ENABLE_SERIAL
SPGEMM_NUMERIC_DECL_MKL_E(Kokkos::Serial)
#endif
#ifdef KOKKOS_ENABLE_OPENMP
SPGEMM_NUMERIC_DECL_MKL_E(Kokkos::OpenMP)
#endif
#endif

}  // namespace Impl
}  // namespace KokkosSparse

#endif
