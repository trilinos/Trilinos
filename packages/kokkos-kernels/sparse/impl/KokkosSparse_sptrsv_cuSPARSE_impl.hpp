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

#ifndef _KOKKOSSPTRSVCUSPARSE_HPP
#define _KOKKOSSPTRSVCUSPARSE_HPP

#include "KokkosSparse_Utils_cusparse.hpp"

namespace KokkosSparse {
namespace Impl {

template <typename KernelHandle, typename ain_row_index_view_type,
          typename ain_nonzero_index_view_type,
          typename ain_values_scalar_view_type>
void sptrsvcuSPARSE_symbolic(KernelHandle* sptrsv_handle,
                             typename KernelHandle::nnz_lno_t nrows,
                             ain_row_index_view_type row_map,
                             ain_nonzero_index_view_type entries,
                             ain_values_scalar_view_type values, bool trans) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#if (CUDA_VERSION >= 11030)
  typedef typename KernelHandle::nnz_lno_t idx_type;
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::scalar_t scalar_type;
  typedef typename KernelHandle::memory_space memory_space;
  typedef typename KernelHandle::nnz_scalar_view_t nnz_scalar_view_t;

  const bool is_cuda_space =
      std::is_same<memory_space, Kokkos::CudaSpace>::value ||
      std::is_same<memory_space, Kokkos::CudaUVMSpace>::value ||
      std::is_same<memory_space, Kokkos::CudaHostPinnedSpace>::value;

  const bool is_idx_type_supported = std::is_same<idx_type, int>::value ||
                                     std::is_same<idx_type, int64_t>::value;

  if (!is_cuda_space) {
    throw std::runtime_error(
        "KokkosKernels sptrsvcuSPARSE_symbolic: MEMORY IS NOT ALLOCATED IN GPU "
        "DEVICE for CUSPARSE\n");
  } else if (!is_idx_type_supported) {
    throw std::runtime_error(
        "CUSPARSE requires local ordinals to be integer (32 bits or 64 "
        "bits).\n");
  } else {
    bool is_lower = sptrsv_handle->is_lower_tri();
    sptrsv_handle->create_cuSPARSE_Handle(trans, is_lower);

    typename KernelHandle::SPTRSVcuSparseHandleType* h =
        sptrsv_handle->get_cuSparseHandle();

    int64_t nnz = static_cast<int64_t>(entries.extent(0));
    size_t pBufferSize;
    void* rm;
    // NOTE (Oct-29-2022):
    // cusparseCreateCsr only supports the same sizes (either 32 bits or 64
    // bits) for row_map_type and entries_type
    if (std::is_same<idx_type, int>::value) {
      if (!std::is_same<size_type, int>::value) {
        sptrsv_handle->allocate_tmp_int_rowmap(row_map.extent(0));
        rm = (void*)sptrsv_handle->get_int_rowmap_ptr_copy(row_map);
      } else {
        rm = (void*)row_map.data();
      }
    } else {  // idx_type has 64 bits
      if (!std::is_same<size_type, int64_t>::value) {
        sptrsv_handle->allocate_tmp_int64_rowmap(row_map.extent(0));
        rm = (void*)sptrsv_handle->get_int64_rowmap_ptr_copy(row_map);
      } else {
        rm = (void*)row_map.data();
      }
    }
    const scalar_type alpha = scalar_type(1.0);

    cusparseIndexType_t cudaCsrRowMapType =
        cusparse_index_type_t_from<idx_type>();
    cusparseIndexType_t cudaCsrColIndType =
        cusparse_index_type_t_from<idx_type>();
    cudaDataType cudaValueType = cuda_data_type_from<scalar_type>();

    // Create sparse matrix in CSR format
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(
        &(h->matDescr), static_cast<int64_t>(nrows),
        static_cast<int64_t>(nrows), nnz, rm, (void*)entries.data(),
        (void*)values.data(), cudaCsrRowMapType, cudaCsrColIndType,
        CUSPARSE_INDEX_BASE_ZERO, cudaValueType));

    // Create dummy dense vector B (RHS)
    nnz_scalar_view_t b_dummy("b_dummy", nrows);
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCreateDnVec(&(h->vecBDescr_dummy), static_cast<int64_t>(nrows),
                            b_dummy.data(), cudaValueType));

    // Create dummy dense vector X (LHS)
    nnz_scalar_view_t x_dummy("x_dummy", nrows);
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCreateDnVec(&(h->vecXDescr_dummy), static_cast<int64_t>(nrows),
                            x_dummy.data(), cudaValueType));

    // Specify Lower|Upper fill mode
    if (is_lower) {
      cusparseFillMode_t fillmode = CUSPARSE_FILL_MODE_LOWER;
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMatSetAttribute(
          h->matDescr, CUSPARSE_SPMAT_FILL_MODE, &fillmode, sizeof(fillmode)));
    } else {
      cusparseFillMode_t fillmode = CUSPARSE_FILL_MODE_UPPER;
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMatSetAttribute(
          h->matDescr, CUSPARSE_SPMAT_FILL_MODE, &fillmode, sizeof(fillmode)));
    }

    // Specify Unit|Non-Unit diagonal type.
    cusparseDiagType_t diagtype = CUSPARSE_DIAG_TYPE_NON_UNIT;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMatSetAttribute(
        h->matDescr, CUSPARSE_SPMAT_DIAG_TYPE, &diagtype, sizeof(diagtype)));

    // Allocate an external buffer for analysis
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpSV_bufferSize(
        h->handle, h->transpose, &alpha, h->matDescr, h->vecBDescr_dummy,
        h->vecXDescr_dummy, cudaValueType, CUSPARSE_SPSV_ALG_DEFAULT,
        h->spsvDescr, &pBufferSize));

    // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void**)&(h->pBuffer), pBufferSize));

    // Run analysis
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpSV_analysis(
        h->handle, h->transpose, &alpha, h->matDescr, h->vecBDescr_dummy,
        h->vecXDescr_dummy, cudaValueType, CUSPARSE_SPSV_ALG_DEFAULT,
        h->spsvDescr, h->pBuffer));

    // Destroy dummy dense vector descriptors
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(h->vecBDescr_dummy));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(h->vecXDescr_dummy));
  }
#else  // CUDA_VERSION < 11030
  typedef typename KernelHandle::nnz_lno_t idx_type;
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::scalar_t scalar_type;
  typedef typename KernelHandle::memory_space memory_space;

  const bool is_cuda_space =
      std::is_same<memory_space, Kokkos::CudaSpace>::value ||
      std::is_same<memory_space, Kokkos::CudaUVMSpace>::value ||
      std::is_same<memory_space, Kokkos::CudaHostPinnedSpace>::value;

  if (!is_cuda_space) {
    throw std::runtime_error(
        "KokkosKernels sptrsvcuSPARSE_symbolic: MEMORY IS NOT ALLOCATED IN GPU "
        "DEVICE for CUSPARSE\n");
  } else if (std::is_same<idx_type, int>::value) {
    bool is_lower = sptrsv_handle->is_lower_tri();
    sptrsv_handle->create_cuSPARSE_Handle(trans, is_lower);

    typename KernelHandle::SPTRSVcuSparseHandleType* h =
        sptrsv_handle->get_cuSparseHandle();

    cusparseStatus_t status;
    status = cusparseCreateCsrsv2Info(&(h->info));
    if (CUSPARSE_STATUS_SUCCESS != status)
      std::cout << "csrsv2info create status error name " << (status)
                << std::endl;

    // query how much memory used in csrsv2, and allocate the buffer
    int nnz = entries.extent_int(0);
    int pBufferSize;

    if (!std::is_same<size_type, int>::value)
      sptrsv_handle->allocate_tmp_int_rowmap(row_map.extent(0));
    const int* rm = !std::is_same<size_type, int>::value
                        ? sptrsv_handle->get_int_rowmap_ptr_copy(row_map)
                        : (const int*)row_map.data();
    const int* ent          = (const int*)entries.data();
    const scalar_type* vals = values.data();

    if (std::is_same<scalar_type, double>::value) {
      cusparseDcsrsv2_bufferSize(h->handle, h->transpose, nrows, nnz, h->descr,
                                 (double*)vals, (int*)rm, (int*)ent, h->info,
                                 &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void**)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name "
                  << cudaGetErrorString(my_error) << std::endl;

      status = cusparseDcsrsv2_analysis(
          h->handle, h->transpose, nrows, nnz, h->descr, (double*)vals,
          (int*)rm, (int*)ent, h->info, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, float>::value) {
      cusparseScsrsv2_bufferSize(h->handle, h->transpose, nrows, nnz, h->descr,
                                 (float*)vals, (int*)rm, (int*)ent, h->info,
                                 &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void**)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name "
                  << cudaGetErrorString(my_error) << std::endl;

      status = cusparseScsrsv2_analysis(
          h->handle, h->transpose, nrows, nnz, h->descr, (float*)vals, (int*)rm,
          (int*)ent, h->info, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, Kokkos::complex<double> >::value) {
      cusparseZcsrsv2_bufferSize(h->handle, h->transpose, nrows, nnz, h->descr,
                                 (cuDoubleComplex*)vals, (int*)rm, (int*)ent,
                                 h->info, &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void**)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name "
                  << cudaGetErrorString(my_error) << std::endl;

      status = cusparseZcsrsv2_analysis(
          h->handle, h->transpose, nrows, nnz, h->descr, (cuDoubleComplex*)vals,
          (int*)rm, (int*)ent, h->info, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, Kokkos::complex<float> >::value) {
      cusparseCcsrsv2_bufferSize(h->handle, h->transpose, nrows, nnz, h->descr,
                                 (cuComplex*)vals, (int*)rm, (int*)ent, h->info,
                                 &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void**)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name "
                  << cudaGetErrorString(my_error) << std::endl;

      status = cusparseCcsrsv2_analysis(
          h->handle, h->transpose, nrows, nnz, h->descr, (cuComplex*)vals,
          (int*)rm, (int*)ent, h->info, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "analysis status error name " << (status) << std::endl;
    } else {
      throw std::runtime_error("CUSPARSE wrapper error: unsupported type.\n");
    }
  } else {
    throw std::runtime_error(
        "CUSPARSE requires local ordinals to be integer.\n");
  }
#endif
#else
  (void)sptrsv_handle;
  (void)nrows;
  (void)row_map;
  (void)entries;
  (void)values;
  (void)trans;
  throw std::runtime_error("CUSPARSE IS NOT DEFINED\n");
  // return;
#endif
}

template <
    typename KernelHandle, typename ain_row_index_view_type,
    typename ain_nonzero_index_view_type, typename ain_values_scalar_view_type,
    typename b_values_scalar_view_type, typename x_values_scalar_view_type>
void sptrsvcuSPARSE_solve(KernelHandle* sptrsv_handle,
                          typename KernelHandle::nnz_lno_t nrows,
                          ain_row_index_view_type row_map,
                          ain_nonzero_index_view_type entries,
                          ain_values_scalar_view_type values,
                          b_values_scalar_view_type rhs,
                          x_values_scalar_view_type lhs, bool /*trans*/
) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#if (CUDA_VERSION >= 11030)
  typedef typename KernelHandle::nnz_lno_t idx_type;
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::scalar_t scalar_type;
  typedef typename KernelHandle::memory_space memory_space;

  const bool is_cuda_space =
      std::is_same<memory_space, Kokkos::CudaSpace>::value ||
      std::is_same<memory_space, Kokkos::CudaUVMSpace>::value ||
      std::is_same<memory_space, Kokkos::CudaHostPinnedSpace>::value;

  const bool is_idx_type_supported = std::is_same<idx_type, int>::value ||
                                     std::is_same<idx_type, int64_t>::value;

  if (!is_cuda_space) {
    throw std::runtime_error(
        "KokkosKernels sptrsvcuSPARSE_solve: MEMORY IS NOT ALLOCATED IN GPU "
        "DEVICE for CUSPARSE\n");
  } else if (!is_idx_type_supported) {
    throw std::runtime_error(
        "CUSPARSE requires local ordinals to be integer (32 bits or 64 "
        "bits).\n");
  } else {
    typename KernelHandle::SPTRSVcuSparseHandleType* h =
        sptrsv_handle->get_cuSparseHandle();

    const scalar_type alpha = scalar_type(1.0);

    cudaDataType cudaValueType = cuda_data_type_from<scalar_type>();

    // Create dense vector B (RHS)
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCreateDnVec(&(h->vecBDescr), static_cast<int64_t>(nrows),
                            (void*)rhs.data(), cudaValueType));

    // Create dense vector X (LHS)
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCreateDnVec(&(h->vecXDescr), static_cast<int64_t>(nrows),
                            (void*)lhs.data(), cudaValueType));

    // Solve
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpSV_solve(
        h->handle, h->transpose, &alpha, h->matDescr, h->vecBDescr,
        h->vecXDescr, cudaValueType, CUSPARSE_SPSV_ALG_DEFAULT, h->spsvDescr));

    // Destroy dense vector descriptors
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(h->vecBDescr));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(h->vecXDescr));
  }
#else  // CUDA_VERSION < 11030
  typedef typename KernelHandle::nnz_lno_t idx_type;
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::scalar_t scalar_type;

  if (std::is_same<idx_type, int>::value) {
    cusparseStatus_t status;

    typename KernelHandle::SPTRSVcuSparseHandleType* h =
        sptrsv_handle->get_cuSparseHandle();

    int nnz = entries.extent_int(0);

    const int* rm = !std::is_same<size_type, int>::value
                        ? sptrsv_handle->get_int_rowmap_ptr()
                        : (const int*)row_map.data();
    const int* ent          = (const int*)entries.data();
    const scalar_type* vals = values.data();
    const scalar_type* bv   = rhs.data();
    scalar_type* xv         = lhs.data();

    if (std::is_same<scalar_type, double>::value) {
      if (h->pBuffer == nullptr) {
        std::cout << "  pBuffer invalid" << std::endl;
      }
      const double alpha = double(1);

      status = cusparseDcsrsv2_solve(h->handle, h->transpose, nrows, nnz,
                                     &alpha, h->descr, (double*)vals, (int*)rm,
                                     (int*)ent, h->info, (double*)bv,
                                     (double*)xv, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, float>::value) {
      if (h->pBuffer == nullptr) {
        std::cout << "  pBuffer invalid" << std::endl;
      }
      const float alpha = float(1);

      status = cusparseScsrsv2_solve(h->handle, h->transpose, nrows, nnz,
                                     &alpha, h->descr, (float*)vals, (int*)rm,
                                     (int*)ent, h->info, (float*)bv, (float*)xv,
                                     h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, Kokkos::complex<double> >::value) {
      cuDoubleComplex cualpha;
      cualpha.x = 1.0;
      cualpha.y = 0.0;
      status    = cusparseZcsrsv2_solve(
          h->handle, h->transpose, nrows, nnz, &cualpha, h->descr,
          (cuDoubleComplex*)vals, (int*)rm, (int*)ent, h->info,
          (cuDoubleComplex*)bv, (cuDoubleComplex*)xv, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, Kokkos::complex<float> >::value) {
      cuComplex cualpha;
      cualpha.x = 1.0;
      cualpha.y = 0.0;
      status    = cusparseCcsrsv2_solve(
          h->handle, h->transpose, nrows, nnz, &cualpha, h->descr,
          (cuComplex*)vals, (int*)rm, (int*)ent, h->info, (cuComplex*)bv,
          (cuComplex*)xv, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status)
        std::cout << "solve status error name " << (status) << std::endl;
    } else {
      throw std::runtime_error("CUSPARSE wrapper error: unsupported type.\n");
    }

  } else {
    throw std::runtime_error(
        "CUSPARSE requires local ordinals to be integer.\n");
  }
#endif
#else
  (void)sptrsv_handle;
  (void)nrows;
  (void)row_map;
  (void)entries;
  (void)values;
  (void)rhs;
  (void)lhs;
  throw std::runtime_error("CUSPARSE IS NOT DEFINED\n");
#endif
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif
