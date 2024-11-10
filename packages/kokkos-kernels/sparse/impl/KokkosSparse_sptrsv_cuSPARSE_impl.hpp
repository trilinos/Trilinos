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

template <typename ExecutionSpace, typename KernelHandle, typename ain_row_index_view_type,
          typename ain_nonzero_index_view_type, typename ain_values_scalar_view_type>
void sptrsvcuSPARSE_symbolic(ExecutionSpace &space, KernelHandle *sptrsv_handle, typename KernelHandle::nnz_lno_t nrows,
                             ain_row_index_view_type row_map, ain_nonzero_index_view_type entries,
                             ain_values_scalar_view_type values, bool trans) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#if (CUDA_VERSION >= 11030)
  typedef typename KernelHandle::nnz_lno_t idx_type;
  typedef typename KernelHandle::size_type size_type;
  typedef typename KernelHandle::scalar_t scalar_type;
  typedef typename KernelHandle::memory_space memory_space;
  typedef typename KernelHandle::nnz_scalar_view_t nnz_scalar_view_t;

  const bool is_cuda_space = std::is_same<memory_space, Kokkos::CudaSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaUVMSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaHostPinnedSpace>::value;

  const bool is_idx_type_supported = std::is_same<idx_type, int>::value || std::is_same<idx_type, int64_t>::value;

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

    typename KernelHandle::SPTRSVcuSparseHandleType *h = sptrsv_handle->get_cuSparseHandle();

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetStream(h->handle, space.cuda_stream()));

    int64_t nnz = static_cast<int64_t>(entries.extent(0));
    size_t pBufferSize;
    void *rm;
    // NOTE (Oct-29-2022):
    // cusparseCreateCsr only supports the same sizes (either 32 bits or 64
    // bits) for row_map_type and entries_type
    if (std::is_same<idx_type, int>::value) {
      if (!std::is_same<size_type, int>::value) {
        sptrsv_handle->allocate_tmp_int_rowmap(row_map.extent(0));
        rm = (void *)sptrsv_handle->get_int_rowmap_ptr_copy(row_map);
      } else {
        rm = (void *)row_map.data();
      }
    } else {  // idx_type has 64 bits
      if (!std::is_same<size_type, int64_t>::value) {
        sptrsv_handle->allocate_tmp_int64_rowmap(row_map.extent(0));
        rm = (void *)sptrsv_handle->get_int64_rowmap_ptr_copy(row_map);
      } else {
        rm = (void *)row_map.data();
      }
    }
    const scalar_type alpha = scalar_type(1.0);

    cusparseIndexType_t cudaCsrRowMapType = cusparse_index_type_t_from<idx_type>();
    cusparseIndexType_t cudaCsrColIndType = cusparse_index_type_t_from<idx_type>();
    cudaDataType cudaValueType            = cuda_data_type_from<scalar_type>();

    // Create sparse matrix in CSR format
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(
        &(h->matDescr), static_cast<int64_t>(nrows), static_cast<int64_t>(nrows), nnz, rm, (void *)entries.data(),
        (void *)values.data(), cudaCsrRowMapType, cudaCsrColIndType, CUSPARSE_INDEX_BASE_ZERO, cudaValueType));

    // Create dummy dense vector B (RHS)
    nnz_scalar_view_t b_dummy(Kokkos::view_alloc(space, "b_dummy"), nrows);
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCreateDnVec(&(h->vecBDescr_dummy), static_cast<int64_t>(nrows), b_dummy.data(), cudaValueType));

    // Create dummy dense vector X (LHS)
    nnz_scalar_view_t x_dummy(Kokkos::view_alloc(space, "x_dummy"), nrows);
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCreateDnVec(&(h->vecXDescr_dummy), static_cast<int64_t>(nrows), x_dummy.data(), cudaValueType));

    // Specify Lower|Upper fill mode
    if (is_lower) {
      cusparseFillMode_t fillmode = CUSPARSE_FILL_MODE_LOWER;
      KOKKOS_CUSPARSE_SAFE_CALL(
          cusparseSpMatSetAttribute(h->matDescr, CUSPARSE_SPMAT_FILL_MODE, &fillmode, sizeof(fillmode)));
    } else {
      cusparseFillMode_t fillmode = CUSPARSE_FILL_MODE_UPPER;
      KOKKOS_CUSPARSE_SAFE_CALL(
          cusparseSpMatSetAttribute(h->matDescr, CUSPARSE_SPMAT_FILL_MODE, &fillmode, sizeof(fillmode)));
    }

    // Specify Unit|Non-Unit diagonal type.
    cusparseDiagType_t diagtype = CUSPARSE_DIAG_TYPE_NON_UNIT;
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseSpMatSetAttribute(h->matDescr, CUSPARSE_SPMAT_DIAG_TYPE, &diagtype, sizeof(diagtype)));

    // Allocate an external buffer for analysis
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpSV_bufferSize(h->handle, h->transpose, &alpha, h->matDescr, h->vecBDescr_dummy,
                                                      h->vecXDescr_dummy, cudaValueType, CUSPARSE_SPSV_ALG_DEFAULT,
                                                      h->spsvDescr, &pBufferSize));

    // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc((void **)&(h->pBuffer), pBufferSize));

    // Run analysis
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpSV_analysis(h->handle, h->transpose, &alpha, h->matDescr, h->vecBDescr_dummy,
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

  const bool is_cuda_space = std::is_same<memory_space, Kokkos::CudaSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaUVMSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaHostPinnedSpace>::value;

  if constexpr (!is_cuda_space) {
    throw std::runtime_error(
        "KokkosKernels sptrsvcuSPARSE_symbolic: MEMORY IS NOT ALLOCATED IN GPU "
        "DEVICE for CUSPARSE\n");
  } else if constexpr (std::is_same<idx_type, int>::value) {
    bool is_lower = sptrsv_handle->is_lower_tri();
    sptrsv_handle->create_cuSPARSE_Handle(trans, is_lower);

    typename KernelHandle::SPTRSVcuSparseHandleType *h = sptrsv_handle->get_cuSparseHandle();

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetStream(h->handle, space.cuda_stream()));

    cusparseStatus_t status;
    status = cusparseCreateCsrsv2Info(&(h->info));
    if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "csrsv2info create status error name " << (status) << std::endl;

    // query how much memory used in csrsv2, and allocate the buffer
    int nnz = entries.extent_int(0);
    int pBufferSize;

    if (!std::is_same<size_type, int>::value) sptrsv_handle->allocate_tmp_int_rowmap(row_map.extent(0));
    const int *rm           = !std::is_same<size_type, int>::value ? sptrsv_handle->get_int_rowmap_ptr_copy(row_map)
                                                                   : (const int *)row_map.data();
    const int *ent          = (const int *)entries.data();
    const scalar_type *vals = values.data();

    if (std::is_same<scalar_type, double>::value) {
      cusparseDcsrsv2_bufferSize(h->handle, h->transpose, nrows, nnz, h->descr, (double *)vals, (int *)rm, (int *)ent,
                                 h->info, &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void **)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name " << cudaGetErrorString(my_error) << std::endl;

      status = cusparseDcsrsv2_analysis(h->handle, h->transpose, nrows, nnz, h->descr, (double *)vals, (int *)rm,
                                        (int *)ent, h->info, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "analysis status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, float>::value) {
      cusparseScsrsv2_bufferSize(h->handle, h->transpose, nrows, nnz, h->descr, (float *)vals, (int *)rm, (int *)ent,
                                 h->info, &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void **)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name " << cudaGetErrorString(my_error) << std::endl;

      status = cusparseScsrsv2_analysis(h->handle, h->transpose, nrows, nnz, h->descr, (float *)vals, (int *)rm,
                                        (int *)ent, h->info, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "analysis status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, Kokkos::complex<double> >::value) {
      cusparseZcsrsv2_bufferSize(h->handle, h->transpose, nrows, nnz, h->descr, (cuDoubleComplex *)vals, (int *)rm,
                                 (int *)ent, h->info, &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void **)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name " << cudaGetErrorString(my_error) << std::endl;

      status = cusparseZcsrsv2_analysis(h->handle, h->transpose, nrows, nnz, h->descr, (cuDoubleComplex *)vals,
                                        (int *)rm, (int *)ent, h->info, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "analysis status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, Kokkos::complex<float> >::value) {
      cusparseCcsrsv2_bufferSize(h->handle, h->transpose, nrows, nnz, h->descr, (cuComplex *)vals, (int *)rm,
                                 (int *)ent, h->info, &pBufferSize);

      // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
      cudaError_t my_error;
      my_error = cudaMalloc((void **)&(h->pBuffer), pBufferSize);

      if (cudaSuccess != my_error)
        std::cout << "cudmalloc pBuffer error_t error name " << cudaGetErrorString(my_error) << std::endl;

      status = cusparseCcsrsv2_analysis(h->handle, h->transpose, nrows, nnz, h->descr, (cuComplex *)vals, (int *)rm,
                                        (int *)ent, h->info, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "analysis status error name " << (status) << std::endl;
    } else {
      throw std::runtime_error("CUSPARSE wrapper error: unsupported type.\n");
    }
  } else {
    throw std::runtime_error("CUSPARSE requires local ordinals to be integer.\n");
  }
#endif
#else
  (void)space;
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

template <typename ExecutionSpace, typename KernelHandle, typename ain_row_index_view_type,
          typename ain_nonzero_index_view_type, typename ain_values_scalar_view_type,
          typename b_values_scalar_view_type, typename x_values_scalar_view_type>
void sptrsvcuSPARSE_solve(ExecutionSpace &space, KernelHandle *sptrsv_handle, typename KernelHandle::nnz_lno_t nrows,
                          ain_row_index_view_type row_map, ain_nonzero_index_view_type entries,
                          ain_values_scalar_view_type values, b_values_scalar_view_type rhs,
                          x_values_scalar_view_type lhs, bool /*trans*/
) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#if (CUDA_VERSION >= 11030)
  typedef typename KernelHandle::nnz_lno_t idx_type;
  typedef typename KernelHandle::scalar_t scalar_type;
  typedef typename KernelHandle::memory_space memory_space;

  (void)row_map;
  (void)entries;
  (void)values;

  const bool is_cuda_space = std::is_same<memory_space, Kokkos::CudaSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaUVMSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaHostPinnedSpace>::value;

  const bool is_idx_type_supported = std::is_same<idx_type, int>::value || std::is_same<idx_type, int64_t>::value;

  if (!is_cuda_space) {
    throw std::runtime_error(
        "KokkosKernels sptrsvcuSPARSE_solve: MEMORY IS NOT ALLOCATED IN GPU "
        "DEVICE for CUSPARSE\n");
  } else if (!is_idx_type_supported) {
    throw std::runtime_error(
        "CUSPARSE requires local ordinals to be integer (32 bits or 64 "
        "bits).\n");
  } else {
    typename KernelHandle::SPTRSVcuSparseHandleType *h = sptrsv_handle->get_cuSparseHandle();

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetStream(h->handle, space.cuda_stream()));

    const scalar_type alpha = scalar_type(1.0);

    cudaDataType cudaValueType = cuda_data_type_from<scalar_type>();

    // Create dense vector B (RHS)
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCreateDnVec(&(h->vecBDescr), static_cast<int64_t>(nrows), (void *)rhs.data(), cudaValueType));

    // Create dense vector X (LHS)
    KOKKOS_CUSPARSE_SAFE_CALL(
        cusparseCreateDnVec(&(h->vecXDescr), static_cast<int64_t>(nrows), (void *)lhs.data(), cudaValueType));

    // Solve
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpSV_solve(h->handle, h->transpose, &alpha, h->matDescr, h->vecBDescr,
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

    typename KernelHandle::SPTRSVcuSparseHandleType *h = sptrsv_handle->get_cuSparseHandle();

    if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetStream(h->handle, space.cuda_stream()));
    }

    int nnz = entries.extent_int(0);

    const int *rm =
        !std::is_same<size_type, int>::value ? sptrsv_handle->get_int_rowmap_ptr() : (const int *)row_map.data();
    const int *ent          = (const int *)entries.data();
    const scalar_type *vals = values.data();
    const scalar_type *bv   = rhs.data();
    scalar_type *xv         = lhs.data();

    if (std::is_same<scalar_type, double>::value) {
      if (h->pBuffer == nullptr) {
        std::cout << "  pBuffer invalid" << std::endl;
      }
      const double alpha = double(1);

      status = cusparseDcsrsv2_solve(h->handle, h->transpose, nrows, nnz, &alpha, h->descr, (double *)vals, (int *)rm,
                                     (int *)ent, h->info, (double *)bv, (double *)xv, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "solve status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, float>::value) {
      if (h->pBuffer == nullptr) {
        std::cout << "  pBuffer invalid" << std::endl;
      }
      const float alpha = float(1);

      status = cusparseScsrsv2_solve(h->handle, h->transpose, nrows, nnz, &alpha, h->descr, (float *)vals, (int *)rm,
                                     (int *)ent, h->info, (float *)bv, (float *)xv, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "solve status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, Kokkos::complex<double> >::value) {
      cuDoubleComplex cualpha;
      cualpha.x = 1.0;
      cualpha.y = 0.0;
      status = cusparseZcsrsv2_solve(h->handle, h->transpose, nrows, nnz, &cualpha, h->descr, (cuDoubleComplex *)vals,
                                     (int *)rm, (int *)ent, h->info, (cuDoubleComplex *)bv, (cuDoubleComplex *)xv,
                                     h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "solve status error name " << (status) << std::endl;
    } else if (std::is_same<scalar_type, Kokkos::complex<float> >::value) {
      cuComplex cualpha;
      cualpha.x = 1.0;
      cualpha.y = 0.0;
      status =
          cusparseCcsrsv2_solve(h->handle, h->transpose, nrows, nnz, &cualpha, h->descr, (cuComplex *)vals, (int *)rm,
                                (int *)ent, h->info, (cuComplex *)bv, (cuComplex *)xv, h->policy, h->pBuffer);

      if (CUSPARSE_STATUS_SUCCESS != status) std::cout << "solve status error name " << (status) << std::endl;
    } else {
      throw std::runtime_error("CUSPARSE wrapper error: unsupported type.\n");
    }

  } else {
    throw std::runtime_error("CUSPARSE requires local ordinals to be integer.\n");
  }
#endif
#else
  (void)space;
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

// --------------------------------
// Stream interface
// --------------------------------

template <class ExecutionSpace, class KernelHandle, class ain_row_index_view_type, class ain_nonzero_index_view_type,
          class ain_values_scalar_view_type, class b_values_scalar_view_type, class x_values_scalar_view_type>
void sptrsvcuSPARSE_solve_streams(const std::vector<ExecutionSpace> &execspace_v, std::vector<KernelHandle> &handle_v,
                                  const std::vector<ain_row_index_view_type> &row_map_v,
                                  const std::vector<ain_nonzero_index_view_type> &entries_v,
                                  const std::vector<ain_values_scalar_view_type> &values_v,
                                  const std::vector<b_values_scalar_view_type> &rhs_v,
                                  std::vector<x_values_scalar_view_type> &lhs_v, bool /*trans*/
) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  using idx_type                 = typename KernelHandle::nnz_lno_t;
  using scalar_type              = typename KernelHandle::nnz_scalar_t;
  using memory_space             = typename KernelHandle::HandlePersistentMemorySpace;
  using sptrsvHandleType         = typename KernelHandle::SPTRSVHandleType;
  using sptrsvCuSparseHandleType = typename sptrsvHandleType::SPTRSVcuSparseHandleType;

  int nstreams = execspace_v.size();
#if (CUDA_VERSION >= 11030)
  (void)row_map_v;
  (void)entries_v;
  (void)values_v;

  const bool is_cuda_space = std::is_same<memory_space, Kokkos::CudaSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaUVMSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaHostPinnedSpace>::value;

  const bool is_idx_type_supported = std::is_same<idx_type, int>::value || std::is_same<idx_type, int64_t>::value;

  if constexpr (!is_cuda_space) {
    throw std::runtime_error(
        "KokkosKernels sptrsvcuSPARSE_solve_streams: MEMORY IS NOT ALLOCATED "
        "IN GPU DEVICE for CUSPARSE\n");
  } else if constexpr (!is_idx_type_supported) {
    throw std::runtime_error(
        "CUSPARSE requires local ordinals to be integer (32 bits or 64 "
        "bits).\n");
  } else {
    const scalar_type alpha = scalar_type(1.0);

    cudaDataType cudaValueType = cuda_data_type_from<scalar_type>();

    std::vector<sptrsvCuSparseHandleType *> h_v(nstreams);

    for (int i = 0; i < nstreams; i++) {
      sptrsvHandleType *sptrsv_handle = handle_v[i].get_sptrsv_handle();
      h_v[i]                          = sptrsv_handle->get_cuSparseHandle();

      // Bind cuspare handle to a stream
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetStream(h_v[i]->handle, execspace_v[i].cuda_stream()));

      int64_t nrows = static_cast<int64_t>(sptrsv_handle->get_nrows());

      // Create dense vector B (RHS)
      KOKKOS_CUSPARSE_SAFE_CALL(
          cusparseCreateDnVec(&(h_v[i]->vecBDescr), nrows, (void *)rhs_v[i].data(), cudaValueType));

      // Create dense vector X (LHS)
      KOKKOS_CUSPARSE_SAFE_CALL(
          cusparseCreateDnVec(&(h_v[i]->vecXDescr), nrows, (void *)lhs_v[i].data(), cudaValueType));
    }

    // Solve
    for (int i = 0; i < nstreams; i++) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpSV_solve(h_v[i]->handle, h_v[i]->transpose, &alpha, h_v[i]->matDescr,
                                                   h_v[i]->vecBDescr, h_v[i]->vecXDescr, cudaValueType,
                                                   CUSPARSE_SPSV_ALG_DEFAULT, h_v[i]->spsvDescr));
    }

    // Destroy dense vector descriptors
    for (int i = 0; i < nstreams; i++) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(h_v[i]->vecBDescr));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(h_v[i]->vecXDescr));
    }
  }
#else  // CUDA_VERSION < 11030
  using size_type = typename KernelHandle::size_type;

  const bool is_cuda_space = std::is_same<memory_space, Kokkos::CudaSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaUVMSpace>::value ||
                             std::is_same<memory_space, Kokkos::CudaHostPinnedSpace>::value;

  if constexpr (!is_cuda_space) {
    throw std::runtime_error(
        "KokkosKernels sptrsvcuSPARSE_solve_streams: MEMORY IS NOT ALLOCATED "
        "IN GPU DEVICE for CUSPARSE\n");
  } else if constexpr (!std::is_same<idx_type, int>::value) {
    throw std::runtime_error("CUSPARSE requires local ordinals to be integer.\n");
  } else {
    const scalar_type alpha = scalar_type(1.0);
    std::vector<sptrsvHandleType *> sptrsv_handle_v(nstreams);
    std::vector<sptrsvCuSparseHandleType *> h_v(nstreams);
    std::vector<const int *> rm_v(nstreams);
    std::vector<const int *> ent_v(nstreams);
    std::vector<const scalar_type *> vals_v(nstreams);
    std::vector<const scalar_type *> bv_v(nstreams);
    std::vector<scalar_type *> xv_v(nstreams);

    for (int i = 0; i < nstreams; i++) {
      sptrsv_handle_v[i] = handle_v[i].get_sptrsv_handle();
      h_v[i]             = sptrsv_handle_v[i]->get_cuSparseHandle();

      // Bind cuspare handle to a stream
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetStream(h_v[i]->handle, execspace_v[i].cuda_stream()));

      if (h_v[i]->pBuffer == nullptr) {
        std::cout << "  pBuffer invalid on stream " << i << std::endl;
      }
      rm_v[i]   = !std::is_same<size_type, int>::value ? sptrsv_handle_v[i]->get_int_rowmap_ptr()
                                                       : reinterpret_cast<const int *>(row_map_v[i].data());
      ent_v[i]  = reinterpret_cast<const int *>(entries_v[i].data());
      vals_v[i] = values_v[i].data();
      bv_v[i]   = rhs_v[i].data();
      xv_v[i]   = lhs_v[i].data();
    }

    for (int i = 0; i < nstreams; i++) {
      int nnz   = entries_v[i].extent_int(0);
      int nrows = static_cast<int>(sptrsv_handle_v[i]->get_nrows());
      if (std::is_same<scalar_type, double>::value) {
        KOKKOS_CUSPARSE_SAFE_CALL(cusparseDcsrsv2_solve(
            h_v[i]->handle, h_v[i]->transpose, nrows, nnz, reinterpret_cast<const double *>(&alpha), h_v[i]->descr,
            reinterpret_cast<const double *>(vals_v[i]), reinterpret_cast<const int *>(rm_v[i]),
            reinterpret_cast<const int *>(ent_v[i]), h_v[i]->info, reinterpret_cast<const double *>(bv_v[i]),
            reinterpret_cast<double *>(xv_v[i]), h_v[i]->policy, h_v[i]->pBuffer));
      } else if (std::is_same<scalar_type, float>::value) {
        KOKKOS_CUSPARSE_SAFE_CALL(cusparseScsrsv2_solve(
            h_v[i]->handle, h_v[i]->transpose, nrows, nnz, reinterpret_cast<const float *>(&alpha), h_v[i]->descr,
            reinterpret_cast<const float *>(vals_v[i]), reinterpret_cast<const int *>(rm_v[i]),
            reinterpret_cast<const int *>(ent_v[i]), h_v[i]->info, reinterpret_cast<const float *>(bv_v[i]),
            reinterpret_cast<float *>(xv_v[i]), h_v[i]->policy, h_v[i]->pBuffer));
      } else if (std::is_same<scalar_type, Kokkos::complex<double> >::value) {
        KOKKOS_CUSPARSE_SAFE_CALL(cusparseZcsrsv2_solve(
            h_v[i]->handle, h_v[i]->transpose, nrows, nnz, reinterpret_cast<const cuDoubleComplex *>(&alpha),
            h_v[i]->descr, reinterpret_cast<const cuDoubleComplex *>(vals_v[i]), reinterpret_cast<const int *>(rm_v[i]),
            reinterpret_cast<const int *>(ent_v[i]), h_v[i]->info, reinterpret_cast<const cuDoubleComplex *>(bv_v[i]),
            reinterpret_cast<cuDoubleComplex *>(xv_v[i]), h_v[i]->policy, h_v[i]->pBuffer));
      } else if (std::is_same<scalar_type, Kokkos::complex<float> >::value) {
        KOKKOS_CUSPARSE_SAFE_CALL(cusparseCcsrsv2_solve(
            h_v[i]->handle, h_v[i]->transpose, nrows, nnz, reinterpret_cast<const cuComplex *>(&alpha), h_v[i]->descr,
            reinterpret_cast<const cuComplex *>(vals_v[i]), reinterpret_cast<const int *>(rm_v[i]),
            reinterpret_cast<const int *>(ent_v[i]), h_v[i]->info, reinterpret_cast<const cuComplex *>(bv_v[i]),
            reinterpret_cast<cuComplex *>(xv_v[i]), h_v[i]->policy, h_v[i]->pBuffer));
      } else {
        throw std::runtime_error("CUSPARSE wrapper error: unsupported type.\n");
      }
    }
  }
#endif
#else
  (void)execspace_v;
  (void)handle_v;
  (void)row_map_v;
  (void)entries_v;
  (void)values_v;
  (void)rhs_v;
  (void)lhs_v;
  throw std::runtime_error("CUSPARSE IS NOT DEFINED\n");
#endif
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif
