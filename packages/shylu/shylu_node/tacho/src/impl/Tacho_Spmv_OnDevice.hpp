// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_SPMV_ON_DEVICE_HPP__
#define __TACHO_SPMV_ON_DEVICE_HPP__

/// \file  Tacho_Spmv_OnDevice.hpp
/// \brief Sparse-matrix dense-vector multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_SupernodeInfo.hpp"
#include "Tacho_TeamFunctor_ExtractCRS.hpp"

#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE))
  // SpMV flag
  #if (CUSPARSE_VERSION >= 11400)
    #define TACHO_CUSPARSE_SPMV_ALG CUSPARSE_SPMV_ALG_DEFAULT
  #else
    #define TACHO_CUSPARSE_SPMV_ALG CUSPARSE_MV_ALG_DEFAULT
  #endif
  // SpMM flag
  #if (CUSPARSE_VERSION >= 11000)
    #define TACHO_CUSPARSE_SPMM_ALG CUSPARSE_SPMM_ALG_DEFAULT
  #else
    #define TACHO_CUSPARSE_SPMM_ALG CUSPARSE_MM_ALG_DEFAULT
  #endif
#elif defined(KOKKOS_ENABLE_HIP)
  #if (ROCM_VERSION >= 60000)
    #define tacho_rocsparse_spmv rocsparse_spmv
  #elif (ROCM_VERSION >= 50400)
    #define tacho_rocsparse_spmv rocsparse_spmv_ex
  #else
    #define tacho_rocsparse_spmv rocsparse_spmv
  #endif
#endif

namespace Tacho {

template<typename supernode_info_type>
struct SpMV {

private:
  using device_type = typename supernode_info_type::device_type;
  using exec_space  = typename device_type::execution_space;

  using host_device_type = typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;
  using host_memory_space = typename host_device_type::memory_space;

  using value_type       = typename supernode_info_type::value_type;
  using value_type_array = typename supernode_info_type::value_type_array;
  using int_type_array   = typename supernode_info_type::int_type_array;

  using rowptr_view = Kokkos::View<int *, device_type>;
  using colind_view = Kokkos::View<int *, device_type>;
  using nzvals_view = Kokkos::View<value_type *, device_type>;

  bool _keep_zeros;
  bool _is_spmv_extracted;
  ordinal_type _nlvls;

  rowptr_view rowptrU;
  colind_view colindU;
  nzvals_view nzvalsU;

  rowptr_view rowptrL;
  colind_view colindL;
  nzvals_view nzvalsL;

  nzvals_view nzvalsD;
  value_type_array  buffer_L;
  value_type_array  buffer_U;
#if defined(KOKKOS_ENABLE_CUDA)
  #if defined(TACHO_HAVE_CUSPARSE)
  // workspace for SpMV
  // (separte for U and L, so that we can "destroy" without waiting for the other)
  cusparseDnMatDescr_t matL, matU, matW;
  cusparseDnVecDescr_t vecL, vecU, vecW;
  cusparseHandle_t sparseHandle;
  #endif
#elif defined(KOKKOS_ENABLE_HIP)
  // workspace for SpMV
  rocsparse_dnmat_descr matL, matU, matW;
  rocsparse_dnvec_descr vecL, vecU, vecW;
  rocsparse_handle sparseHandle;
#endif

  int _status;
  inline void check_OnDeviceSPMV_Status(const char *func, const char *lib) {
    if (_status != 0) {
      printf("Error: %s, %s returns non-zero status %d\n", lib, func, _status);
      std::runtime_error("checkStatus failed");
    }
  }

public:
  SpMV() {
    _keep_zeros = false;
    _is_spmv_extracted = false;
    _nlvls = 0;
  }

  SpMV(bool keep_zeros) {
    _keep_zeros = keep_zeros;
    _is_spmv_extracted = false;
    _nlvls = 0;
  }

  template <typename supernode_type, typename multi_vectors_type>
  int ApplyL_OnDevice(const ordinal_type lvl,
                            supernode_type &s0,
                      const multi_vectors_type &t,
                            multi_vectors_type &w) {
    int r_val(0);
    const ordinal_type m = t.extent(0);
    const ordinal_type nrhs = t.extent(1);
    const ordinal_type old_nrhs = w.extent(1);

    if (old_nrhs != nrhs) {
      // expand workspace
      Kokkos::resize(w, m, nrhs);
    }
#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE)) || \
     defined(KOKKOS_ENABLE_HIP)
    const ordinal_type ldt = t.stride(1);

#if defined(KOKKOS_ENABLE_CUDA)
    cudaDataType computeType = CUDA_R_64F;
    if (std::is_same<value_type, float>::value) {
      computeType = CUDA_R_32F;
    } else if (!std::is_same<value_type, double>::value) {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_datatype rocsparse_compute_type = rocsparse_datatype_f64_r;
    if (std::is_same<value_type, float>::value) {
      rocsparse_compute_type = rocsparse_datatype_f32_r;
    } else if (!std::is_same<value_type, double>::value) {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }
#endif
    // compute t = L^{-1}*w
    const value_type alpha (1);
    const value_type beta  (0);
    if (old_nrhs != nrhs) {
      // attach to Cusparse/Rocsparse data struct
      int ldw = w.stride(1);
#if defined(KOKKOS_ENABLE_CUDA)
      // destroy previous
      cusparseDestroyDnMat(matW);
      cusparseDestroyDnVec(vecW);
      // create new
      cusparseCreateDnMat(&matW, m, nrhs, ldw, (void*)(w.data()), computeType, CUSPARSE_ORDER_COL);
      cusparseCreateDnVec(&vecW, m, (void*)(w.data()), computeType);
#elif defined(KOKKOS_ENABLE_HIP)
      // destroy previous
      rocsparse_destroy_dnmat_descr(matW);
      rocsparse_destroy_dnvec_descr(vecW);
      // create new
      rocsparse_create_dnmat_descr(&matW, m, nrhs, ldw, (void*)(w.data()), rocsparse_compute_type, rocsparse_order_column);
      rocsparse_create_dnvec_descr(&vecW, m, (void*)(w.data()), rocsparse_compute_type);
#endif
    }
#if defined(KOKKOS_ENABLE_CUDA)
    // Desctory old CSR
    cusparseDestroySpMat(s0.L_cusparse);
    // Re-create CuSparse CSR
    if (s0.spmv_explicit_transpose) {
      cusparseCreateCsr(&s0.L_cusparse, m, m, s0.nnzL,
                        s0.rowptrL, s0.colindL, s0.nzvalsL,
                        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                        CUSPARSE_INDEX_BASE_ZERO, computeType);
    } else {
      cusparseCreateCsr(&s0.L_cusparse, m, m, s0.nnzU,
                        s0.rowptrU, s0.colindU, s0.nzvalsU,
                        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                        CUSPARSE_INDEX_BASE_ZERO, computeType);
    }
    // Call SpMV/SPMM
    cusparseStatus_t status;
    cusparseOperation_t opL = (s0.spmv_explicit_transpose ? CUSPARSE_OPERATION_NON_TRANSPOSE : CUSPARSE_OPERATION_TRANSPOSE);
    if (nrhs > 1) {
      if (lvl == _nlvls-1) {
        // start : destroy previous
        cusparseDestroyDnMat(matL);
        // start : create DnMat for T
        cusparseCreateDnMat(&matL, m, nrhs, ldt, (void*)(t.data()), computeType, CUSPARSE_ORDER_COL);
      }
      // create vectors
      auto matX = ((_nlvls-1-lvl)%2 == 0 ? matL : matW);
      auto matY = ((_nlvls-1-lvl)%2 == 0 ? matW : matL);
      // SpMM
      status = cusparseSpMM(sparseHandle, opL, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha, s0.L_cusparse,
                                    matX,
                            &beta,  matY,
                            computeType, TACHO_CUSPARSE_SPMM_ALG, (void*)buffer_L.data());
    } else {
      if (lvl == _nlvls-1) {
        // start : destroy previous
        cusparseDestroyDnVec(vecL);
        // start : create DnMat for T
        cusparseCreateDnVec(&vecL, m, (void*)(t.data()), computeType);
      }
      // create vectors
      auto vecX = ((_nlvls-1-lvl)%2 == 0 ? vecL : vecW);
      auto vecY = ((_nlvls-1-lvl)%2 == 0 ? vecW : vecL);
      // SpMV
      status = cusparseSpMV(sparseHandle, opL,
                            &alpha, s0.L_cusparse,
                                    vecX,
                            &beta,  vecY,
                            computeType, TACHO_CUSPARSE_SPMV_ALG, (void*)buffer_L.data());
    }
    if (CUSPARSE_STATUS_SUCCESS != status) {
      printf( " Failed cusparseSpMV for SpMV (lower)\n" );
    }
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_status status;
    if (nrhs > 1) {
      if (lvl == _nlvls-1) {
        // start : destroy previous
        rocsparse_destroy_dnmat_descr(matL);
        // start : create DnMat for T
        rocsparse_create_dnmat_descr(&matL, m, nrhs, ldt, (void*)(t.data()), rocsparse_compute_type, rocsparse_order_column);
      }
      // create vectors
      auto vecX = ((_nlvls-1-lvl)%2 == 0 ? matL : matW);
      auto vecY = ((_nlvls-1-lvl)%2 == 0 ? matW : matL);
      if (s0.spmv_explicit_transpose) {
        size_t buffer_size = buffer_L.extent(0);
        status = rocsparse_spmm(sparseHandle, rocsparse_operation_none, rocsparse_operation_none,
                                &alpha, s0.descrL, vecX, &beta, vecY,
                                rocsparse_compute_type, rocsparse_spmm_alg_default,
                                rocsparse_spmm_stage_compute,
                                &buffer_size, (void*)buffer_L.data());
      } else {
        size_t buffer_size = buffer_L.extent(0);
        status = rocsparse_spmm(sparseHandle, rocsparse_operation_transpose, rocsparse_operation_none,
                                &alpha, s0.descrL, vecX, &beta, vecY, // dscrL stores the same ptrs as descrU, but optimized for trans
                                rocsparse_compute_type, rocsparse_spmm_alg_default,
                                rocsparse_spmm_stage_compute,
                                &buffer_size, (void*)buffer_L.data());
      }
    } else {
      if (lvl == _nlvls-1) {
        // start : destroy previous
        rocsparse_destroy_dnvec_descr(vecL);
        // start : create DnVec for T
        rocsparse_create_dnvec_descr(&vecL, m, (void*)(t.data()), rocsparse_compute_type);
      }
      size_t buffer_size = buffer_L.extent(0);
      auto vecX = ((_nlvls-1-lvl)%2 == 0 ? vecL : vecW);
      auto vecY = ((_nlvls-1-lvl)%2 == 0 ? vecW : vecL);
      if (s0.spmv_explicit_transpose) {
        status = tacho_rocsparse_spmv
          (sparseHandle, rocsparse_operation_none,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_compute,
           #endif
           &buffer_size, (void*)buffer_L.data());
      } else {
        status = tacho_rocsparse_spmv
          (sparseHandle, rocsparse_operation_transpose,
           &alpha, s0.descrL, vecX, &beta, vecY, // dscrL stores the same ptrs as descrU, but optimized for trans
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_compute,
           #endif
           &buffer_size, (void*)buffer_L.data());
      }
    }
    if (rocsparse_status_success != status) {
      printf( " Failed rocsparse_spmv for L\n" );
    }
#endif
#else
    const value_type zero(0);
    auto h_w = Kokkos::create_mirror_view_and_copy(host_memory_space(), ((_nlvls-1-lvl)%2 == 0 ? t : w));
    auto h_t = Kokkos::create_mirror_view(host_memory_space(), ((_nlvls-1-lvl)%2 == 0 ? w : t));
    Kokkos::deep_copy(h_t, zero);

    if (s0.spmv_explicit_transpose) {
      UnmanagedViewType<int_type_array>    d_rowptrL(s0.rowptrL, m+1);
      UnmanagedViewType<int_type_array>    d_colindL(s0.colindL, s0.nnzL);
      UnmanagedViewType<value_type_array>  d_nzvalsL(s0.nzvalsL, s0.nnzL);
      auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrL);
      auto h_colind = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_colindL);
      auto h_nzvals = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nzvalsL);
      for (ordinal_type i = 0; i < m ; i++) {
        for (int k = h_rowptr(i); k < h_rowptr(i+1); k++) {
          for (int j = 0; j < nrhs; j++) {
            h_t(i, j) += h_nzvals(k) * h_w(h_colind(k), j);
          }
        }
      }
    } else {
      UnmanagedViewType<int_type_array>    d_rowptrU(s0.rowptrU, m+1);
      UnmanagedViewType<int_type_array>    d_colindU(s0.colindU, s0.nnzU);
      UnmanagedViewType<value_type_array>  d_nzvalsU(s0.nzvalsU, s0.nnzU);
      auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrU);
      auto h_colind = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_colindU);
      auto h_nzvals = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nzvalsU);
      for (ordinal_type i = 0; i < m ; i++) {
        for (int k = h_rowptr(i); k < h_rowptr(i+1); k++) {
          for (int j = 0; j < nrhs; j++) {
            h_t(h_colind(k), j) += h_nzvals(k) * h_w(i, j);
          }
        }
      }
    }
    if ((_nlvls-1-lvl)%2 == 0) {
      Kokkos::deep_copy(w, h_t);
    } else {
      Kokkos::deep_copy(t, h_t);
    }
#endif
    if (lvl == 0) {
      // end : copy to output
      if ((_nlvls-1)%2 == 0) {
        Kokkos::deep_copy(t, w);
      }
    }
    return r_val;
  }

  template <typename supernode_type,
            typename multi_vectors_type>
  int ApplyU_OnDevice(const ordinal_type lvl,
                            supernode_type &s0,
                      const multi_vectors_type &t,
                            multi_vectors_type &w) {
    int r_val(0);
    const ordinal_type m = t.extent(0);
    const ordinal_type nrhs = t.extent(1);

#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE)) || \
     defined(KOKKOS_ENABLE_HIP)
    // x = t & y = w (lvl = 0,2,4)
    // compute t = L^{-1}*w
    const value_type alpha (1);
    const value_type beta  (0);
    const ordinal_type ldt = t.stride(1);
#if defined(KOKKOS_ENABLE_CUDA)
    cudaDataType computeType = CUDA_R_64F;
    if (std::is_same<value_type, float>::value) {
      computeType = CUDA_R_32F;
    } else if (!std::is_same<value_type, double>::value) {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }

    cusparseStatus_t status;
    // Desctory old CSR
    cusparseDestroySpMat(s0.U_cusparse);
    // Re-create CuSparse CSR
    cusparseCreateCsr(&s0.U_cusparse, m, m, s0.nnzU,
                      s0.rowptrU, s0.colindU, s0.nzvalsU,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                      CUSPARSE_INDEX_BASE_ZERO, computeType);

    // Call SpMV/SPMM
    if (nrhs > 1) {
      if (lvl == 0) {
        // start : destroy previous
        cusparseDestroyDnMat(matU);
        // start : create DnMat for T
        cusparseCreateDnMat(&matU, m, nrhs, ldt, (void*)(t.data()), computeType, CUSPARSE_ORDER_COL);
      }
      auto vecX = (lvl%2 == 0 ? matU : matW);
      auto vecY = (lvl%2 == 0 ? matW : matU);
      // SpMM
      status = cusparseSpMM(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha, s0.U_cusparse,
                                    vecX,
                            &beta,  vecY,
                            computeType, TACHO_CUSPARSE_SPMM_ALG, (void*)buffer_U.data());
    } else {
      if (lvl == 0) {
        // start : destroy previous
        cusparseDestroyDnVec(vecU);
        // start : create DnMat for T
        cusparseCreateDnVec(&vecU, m, (void*)(t.data()), computeType);
      }
      auto vecX = (lvl%2 == 0 ? vecU : vecW);
      auto vecY = (lvl%2 == 0 ? vecW : vecU);
      // SpMV
      status = cusparseSpMV(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha, s0.U_cusparse,
                                    vecX,
                            &beta,  vecY,
                            computeType, TACHO_CUSPARSE_SPMV_ALG, (void*)buffer_U.data());
    }
    if (CUSPARSE_STATUS_SUCCESS != status) {
       printf( " Failed cusparseSpMV for SpMV (upper)\n" );
    }
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_datatype rocsparse_compute_type = rocsparse_datatype_f64_r;
    if (std::is_same<value_type, float>::value) {
      rocsparse_compute_type = rocsparse_datatype_f32_r;
    } else if (!std::is_same<value_type, double>::value) {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }
    size_t buffer_size = buffer_U.extent(0);
    rocsparse_status status;
    if (nrhs > 1) {
      if (lvl == 0) {
        // start : create DnMat for T
        rocsparse_destroy_dnmat_descr(matU);
        rocsparse_create_dnmat_descr(&matU, m, nrhs, ldt, (void*)(t.data()), rocsparse_compute_type, rocsparse_order_column);
      }
      auto vecX = (lvl%2 == 0 ? matU : matW);
      auto vecY = (lvl%2 == 0 ? matW : matU);
      status = rocsparse_spmm(sparseHandle, rocsparse_operation_none, rocsparse_operation_none,
                              &alpha, s0.descrU, vecX, &beta, vecY,
                              rocsparse_compute_type, rocsparse_spmm_alg_default,
                              rocsparse_spmm_stage_compute,
                              &buffer_size, (void*)buffer_U.data());
    } else {
      if (lvl == 0) {
        // start : create DnVec for T
        rocsparse_destroy_dnvec_descr(vecU);
        rocsparse_create_dnvec_descr(&vecU, m, (void*)(t.data()), rocsparse_compute_type);
      }
      auto vecX = (lvl%2 == 0 ? vecU : vecW);
      auto vecY = (lvl%2 == 0 ? vecW : vecU);
      status = tacho_rocsparse_spmv
          (sparseHandle, rocsparse_operation_none,
           &alpha, s0.descrU, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_compute,
           #endif
           &buffer_size, (void*)buffer_U.data());
    }
    if (rocsparse_status_success != status) {
      printf( " Failed rocsparse_spmv for U\n" );
    }
#endif
#else
    const value_type zero(0);
    auto h_w = Kokkos::create_mirror_view_and_copy(host_memory_space(), (lvl%2 == 0 ? t : w));
    auto h_t = Kokkos::create_mirror_view(host_memory_space(), (lvl%2 == 0 ? w : t));
    Kokkos::deep_copy(h_t, zero);

    UnmanagedViewType<int_type_array>    d_rowptrU(s0.rowptrU, m+1);
    UnmanagedViewType<int_type_array>    d_colindU(s0.colindU, s0.nnzU);
    UnmanagedViewType<value_type_array>  d_nzvalsU(s0.nzvalsU, s0.nnzU);
    auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrU);
    auto h_colind = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_colindU);
    auto h_nzvals = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nzvalsU);

    for (ordinal_type i = 0; i < m ; i++) {
      for (int k = h_rowptr(i); k < h_rowptr(i+1); k++) {
        for (int j = 0; j < nrhs; j++) {
          h_t(i, j) += h_nzvals(k) * h_w(h_colind(k), j);
        }
      }
    }
    if (lvl%2 == 0) {
      Kokkos::deep_copy(w, h_t);
    } else {
      Kokkos::deep_copy(t, h_t);
    }
#endif
    if (lvl == _nlvls-1) {
      // end : copy to output
      if (lvl%2 == 0) {
        Kokkos::deep_copy(t, w);
      }
    }
    return r_val;
  }

  template <typename size_type_array_host, typename ordinal_type_array_host, typename supernode_type_array_host>
  int Release(const bool release_all, const bool verbose,
              const size_type_array_host &h_level_ptr,
              const ordinal_type_array_host &h_level_sids,
              const supernode_type_array_host &h_supernodes) {
    int r_val(0);
    if (verbose) {
      printf("LevelSetTools:releaseCRS\n");
      printf("========================\n");
      printf(" Have%s been extracted\n", (_is_spmv_extracted == 0 ?  " not" : ""));
      if (release_all) printf(" Release all\n");
      printf("\n"); fflush(stdout);
    }
    if(_is_spmv_extracted != 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      Kokkos::fence();
      if (release_all) {
        for (ordinal_type lvl = 0; lvl < _nlvls; ++lvl) {
          const ordinal_type pbeg = h_level_ptr(lvl);
          // the first supernode in this lvl (where the CRS matrix is stored)
          auto &s0 = h_supernodes(h_level_sids(pbeg));
#if defined(KOKKOS_ENABLE_CUDA)
          cusparseDestroySpMat(s0.U_cusparse);
          cusparseDestroySpMat(s0.L_cusparse);
#elif defined(KOKKOS_ENABLE_HIP)
          rocsparse_destroy_spmat_descr(s0.descrU);
          rocsparse_destroy_spmat_descr(s0.descrL);
#endif
        }
      }
#if defined(TACHO_HAVE_CUSPARSE) && defined(KOKKOS_ENABLE_CUDA)
      cusparseDestroy(sparseHandle);
      cusparseDestroyDnMat(matL);
      cusparseDestroyDnVec(vecL);
      cusparseDestroyDnMat(matU);
      cusparseDestroyDnVec(vecU);
      cusparseDestroyDnMat(matW);
      cusparseDestroyDnVec(vecW);
#elif defined(KOKKOS_ENABLE_HIP)
      rocsparse_destroy_handle(sparseHandle);
      rocsparse_destroy_dnmat_descr(matL);
      rocsparse_destroy_dnvec_descr(vecL);
      rocsparse_destroy_dnmat_descr(matU);
      rocsparse_destroy_dnvec_descr(vecU);
      rocsparse_destroy_dnmat_descr(matW);
      rocsparse_destroy_dnvec_descr(vecW);
#endif
#endif
      _is_spmv_extracted = 0;
    }
    return r_val;
  }


  template <typename size_type_array_host, typename ordinal_type_array_host, typename ordinal_type_array, 
            typename supernode_type_array_host,
            typename stream_type,
            typename multi_vectors_type>
  int Setup(const bool store_transpose, const bool verbose,
            const int method, const ordinal_type m, const ordinal_type nlvls,
            const size_type_array_host &h_level_ptr,
            const ordinal_type_array &level_sids,
            const ordinal_type_array_host &h_level_sids,
            const ordinal_type_array_host &h_solve_mode,
            const supernode_info_type &supernode_info,
            const supernode_type_array_host &h_supernodes,
            const ordinal_type_array &solve_mode,
            ordinal_type_array &piv,
	    const stream_type &stream_0,
            multi_vectors_type &_w_vec) {
    using host_device_type = typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;
    using host_memory_space = typename host_device_type::memory_space;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    int r_val(0);
    if (verbose) {
      printf("LevelSetTools:setupCRS(method=%d)\n",method);
      printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    }
    const bool lu = (method == 3);
    const bool ldl = (method == 2);
    const bool ldl_nopiv = (method == 0);
    const ordinal_type nrhs = 1;

    // ========================
    // free CRS,
    // if it has been extracted
    _nlvls = nlvls; // store # of levels
    Release(true, verbose,
            h_level_ptr, h_level_sids, h_supernodes);

    // ========================
    // workspace
    Kokkos::resize(_w_vec, m, nrhs);
#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE)) || \
     defined(KOKKOS_ENABLE_HIP)
    const value_type zero(0);

    int ldw = _w_vec.stride(1);
#if defined(KOKKOS_ENABLE_CUDA)
    cudaDataType computeType;
    if (std::is_same<value_type, double>::value) {
      computeType = CUDA_R_64F;
    } else if (std::is_same<value_type, float>::value) {
      computeType = CUDA_R_32F;
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }
    // create cusparse handle
    cusparseCreate(&sparseHandle);
    // attach handle to stream-0
    cusparseSetStream(sparseHandle, stream_0);
    // attach to Cusparse data struct
    cusparseCreateDnMat(&matW, m, nrhs, ldw, (void*)(_w_vec.data()), computeType, CUSPARSE_ORDER_COL);
    cusparseCreateDnVec(&vecW, m, (void*)(_w_vec.data()), computeType);
    // also to T, to be destroyed before each SpMV call
    cusparseCreateDnMat(&matL, m, nrhs, ldw, (void*)(_w_vec.data()), computeType, CUSPARSE_ORDER_COL);
    cusparseCreateDnVec(&vecL, m, (void*)(_w_vec.data()), computeType);
    cusparseCreateDnMat(&matU, m, nrhs, ldw, (void*)(_w_vec.data()), computeType, CUSPARSE_ORDER_COL);
    cusparseCreateDnVec(&vecU, m, (void*)(_w_vec.data()), computeType);
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_datatype rocsparse_compute_type = rocsparse_datatype_f64_r;
    if (std::is_same<value_type, float>::value) {
      rocsparse_compute_type = rocsparse_datatype_f32_r;
    }
    // create rocsparse handle
    _status = rocsparse_create_handle(&sparseHandle);
    check_OnDeviceSPMV_Status("rocsparse_create_handle", "rocsparse");
    // attach handle to stream-0
    rocsparse_set_stream(sparseHandle, stream_0);
    check_OnDeviceSPMV_Status("rocsparse_create_handle", "rocsparse");
    // attach to Rocsparse data struct
    _status = rocsparse_create_dnmat_descr(&matW, m, nrhs, ldw, (void*)(_w_vec.data()), rocsparse_compute_type, rocsparse_order_column);
    check_OnDeviceSPMV_Status("rocsparse_dnmat_descr", "rocsparse");
    _status = rocsparse_create_dnvec_descr(&vecW, m, (void*)(_w_vec.data()), rocsparse_compute_type);
    check_OnDeviceSPMV_Status("rocsparse_dnmat_descr", "rocsparse");
    // also to T, to be destroyed before each SpMV call
    _status = rocsparse_create_dnmat_descr(&matL, m, nrhs, ldw, (void*)(_w_vec.data()), rocsparse_compute_type, rocsparse_order_column);
    check_OnDeviceSPMV_Status("rocsparse_dnmat_descr", "rocsparse");
    _status = rocsparse_create_dnvec_descr(&vecL, m, (void*)(_w_vec.data()), rocsparse_compute_type);
    check_OnDeviceSPMV_Status("rocsparse_dnvec_descr", "rocsparse");
    _status = rocsparse_create_dnmat_descr(&matU, m, nrhs, ldw, (void*)(_w_vec.data()), rocsparse_compute_type, rocsparse_order_column);
    check_OnDeviceSPMV_Status("rocsparse_dnmat_descr", "rocsparse");
    _status = rocsparse_create_dnvec_descr(&vecU, m, (void*)(_w_vec.data()), rocsparse_compute_type);
    check_OnDeviceSPMV_Status("rocsparse_dnvec_descr", "rocsparse");
#endif
#endif

    // allocate rowptrs
    Kokkos::resize(rowptrU, nlvls*(1+m));
    Kokkos::resize(rowptrL, nlvls*(1+m));
    Kokkos::deep_copy(rowptrL, 0);
    // counting nnz, first, so that we can allocate in NumericalTool
    size_t ptr = 0;
    size_t nnzU = 0;
    size_t nnzL = 0;
    typedef TeamFunctor_ExtractCrs<supernode_info_type> functor_type;
    for (ordinal_type lvl = 0; lvl < nlvls; ++lvl) {
      const ordinal_type pbeg = h_level_ptr(lvl), pend = h_level_ptr(lvl + 1);

      // the first supernode in this lvl (where the CRS matrix is stored)
      auto &s0 = h_supernodes(h_level_sids(pbeg));
      s0.spmv_explicit_transpose = store_transpose;

      // ========================
      // count nnz / row
      auto d_rowptrU = Kokkos::subview(rowptrU, range_type(ptr, ptr+m+1));
      s0.rowptrU = d_rowptrU.data();

      functor_type extractor_crs(_keep_zeros, supernode_info, solve_mode, level_sids);
      extractor_crs.setGlobalSize(m);
      extractor_crs.setRange(pbeg, pend);
      extractor_crs.setRowPtr(s0.rowptrU);
      if (ldl) {
        // integrate pivoting into U
        extractor_crs.integratePivots(true, piv);
      }
      {
        using team_policy_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                                    typename functor_type::ExtractPtrTag>;
        team_policy_type team_policy((pend-pbeg)+1, Kokkos::AUTO());

        Kokkos::parallel_for("extract rowptr", team_policy, extractor_crs);
        exec_space().fence();
      }

      // ========================
      // shift to generate rowptr
      {
        using range_policy_type = Kokkos::RangePolicy<exec_space>;
        Kokkos::parallel_scan("shiftRowptr", range_policy_type(0, m+1), rowptr_sum(s0.rowptrU));
        exec_space().fence();
        // get nnz
        auto d_nnz = Kokkos::subview(d_rowptrU, range_type(m, m+1));
        auto h_nnz = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nnz);
        s0.nnzU = h_nnz(0);
        nnzU += s0.nnzU;
      }

      if (lu) {
        // get nnz per row (L is stored by column)
        auto d_rowptrL = Kokkos::subview(rowptrL, range_type(ptr, ptr+m+1));
        s0.rowptrL = d_rowptrL.data();
        {
          using team_policy_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                                      typename functor_type::ExtractPtrColTag>;
          team_policy_type team_policy((pend-pbeg)+1, Kokkos::AUTO());

          extractor_crs.setRowPtr(s0.rowptrL);
          Kokkos::parallel_for("extract rowptr L", team_policy, extractor_crs);
          exec_space().fence();
        }
        {
          // convert to offset
          using range_policy_type = Kokkos::RangePolicy<exec_space>;
          Kokkos::parallel_scan("shiftRowptr L", range_policy_type(0, m+1), rowptr_sum(s0.rowptrL));
          exec_space().fence();
          // get nnz (on CPU for now)
          auto d_nnz = Kokkos::subview(d_rowptrL, range_type(m, m+1));
          auto h_nnz = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nnz);
          s0.nnzL = h_nnz(0);
          nnzL += s0.nnzL;
        }
        s0.spmv_explicit_transpose = true;
      } else if (s0.spmv_explicit_transpose) {
        // ========================
        // explicitly form transpose
        s0.nnzL = s0.nnzU;
        auto d_rowptrL = Kokkos::subview(rowptrL, range_type(ptr, ptr+m+1));
        s0.rowptrL = d_rowptrL.data();
        nnzL += s0.nnzL;
      }
      ptr += (1+m);
    }
    // allocate (TODO: move to symbolic)
    Kokkos::resize(colindU, nnzU);
    Kokkos::resize(nzvalsU, nnzU);
    Kokkos::resize(colindL, nnzL);
    Kokkos::resize(nzvalsL, nnzL);
    if (ldl_nopiv) {
      // Initialized to be I, then updated with diagonal values during numeric
      //  The location of the updated values should be the same / symbolic
      const value_type one(1);
      Kokkos::resize(nzvalsD, nlvls*m);
      Kokkos::deep_copy(nzvalsD, one);
    }
    _is_spmv_extracted = true;

    return r_val;
  }

  template <typename size_type_array_host, typename ordinal_type_array_host, typename ordinal_type_array,
            typename supernode_type_array_host,
            typename stream_type,
            typename multi_vectors_type>
  int Load(const bool store_transpose, const bool verbose,
           const int method, const ordinal_type m,
           const size_type_array_host &h_level_ptr,
           const ordinal_type_array &level_sids,
           const ordinal_type_array_host &h_level_sids,
           const ordinal_type_array_host &h_solve_mode,
           const supernode_info_type &supernode_info,
           const supernode_type_array_host &h_supernodes,
           const ordinal_type_array &solve_mode,
           ordinal_type_array &piv,
	   const stream_type &stream_0,
           multi_vectors_type &_w_vec) {

    using exec_space = typename multi_vectors_type::execution_space;
    using host_device_type = typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type;
    using host_memory_space = typename host_device_type::memory_space;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    const bool lu = (method == 3);
    const bool ldl = (method == 2);
    const bool ldl_nopiv = (method == 0);

    int r_val(0);
    // ========================
    // workspace
#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE)) || \
     defined(KOKKOS_ENABLE_HIP)
#if defined(KOKKOS_ENABLE_CUDA)
    // vectors used for preprocessing
    cudaDataType computeType;
    if (std::is_same<value_type, double>::value) {
      computeType = CUDA_R_64F;
    } else if (std::is_same<value_type, float>::value) {
      computeType = CUDA_R_32F;
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }
#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
    const ordinal_type nrhs = 1;
    cusparseDnMatDescr_t vecX, vecY;
    const ordinal_type ldx = _w_vec.stride(1);
    cusparseCreateDnMat(&vecX, m, nrhs, ldx, _w_vec.data(), computeType, CUSPARSE_ORDER_COL);
    cusparseCreateDnMat(&vecY, m, nrhs, ldx, _w_vec.data(), computeType, CUSPARSE_ORDER_COL);
#else
    cusparseDnVecDescr_t vecX, vecY;
    cusparseCreateDnVec(&vecX, m, _w_vec.data(), computeType);
    cusparseCreateDnVec(&vecY, m, _w_vec.data(), computeType);
#endif
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_datatype rocsparse_compute_type = rocsparse_datatype_f64_r;
    if (std::is_same<value_type, float>::value) {
      rocsparse_compute_type = rocsparse_datatype_f32_r;
    }
    // vectors used for preprocessing
    rocsparse_dnvec_descr vecX, vecY;
    rocsparse_create_dnvec_descr(&vecX, m, (void*)_w_vec.data(), rocsparse_compute_type);
    rocsparse_create_dnvec_descr(&vecY, m, (void*)_w_vec.data(), rocsparse_compute_type);
#endif
#endif
    size_t ptr = 0;
    size_t nnzU = 0;
    size_t nnzL = 0;
    typedef TeamFunctor_ExtractCrs<supernode_info_type> functor_type;

    // load nonzero val/ind
    for (ordinal_type lvl = 0; lvl < _nlvls; ++lvl) {
      const ordinal_type pbeg = h_level_ptr(lvl), pend = h_level_ptr(lvl + 1);

      // the first supernode in this lvl (where the CRS matrix is stored)
      auto &s0 = h_supernodes(h_level_sids(pbeg));

      // ========================
      // assign memory
      auto d_rowptrU = Kokkos::subview(rowptrU, range_type(ptr, ptr+m+1));
      auto d_colindU = Kokkos::subview(colindU, range_type(nnzU, nnzU+s0.nnzU));
      auto d_nzvalsU = Kokkos::subview(nzvalsU, range_type(nnzU, nnzU+s0.nnzU));
      s0.colindU = d_colindU.data();
      s0.nzvalsU = d_nzvalsU.data();
      nnzU += s0.nnzU;

      // ========================
      // extract nonzero element
      functor_type extractor_crs(_keep_zeros, supernode_info, solve_mode, level_sids);
      extractor_crs.setGlobalSize(m);
      extractor_crs.setRange(pbeg, pend);
      extractor_crs.setRowPtr(s0.rowptrU);
      extractor_crs.setCrsView(s0.colindU, s0.nzvalsU);
      if (ldl) {
        // integrate pivoting into U
        extractor_crs.integratePivots(true, piv);
      } else if (ldl_nopiv) {
        auto d_nzvalsD = Kokkos::subview(nzvalsD, range_type(lvl*m, (lvl+1)*m));
        s0.nzvalsD = d_nzvalsD.data();
        extractor_crs.withUnitDiag(s0.nzvalsD);
      }
      {
        using team_policy_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                                    typename functor_type::ExtractValTag>;
        team_policy_type team_policy((pend-pbeg)+1, Kokkos::AUTO());

        // >> launch functor to extract nonzero entries
        Kokkos::parallel_for("extract nzvals", team_policy, extractor_crs);
        exec_space().fence();
      }

      // ========================
      // shift back (TODO: shift first to avoid this)
      {
        //  copy to CPU, for now
        auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrU);
        for (ordinal_type i = m; i > 0 ; i--) h_rowptr(i) = h_rowptr(i-1);
        h_rowptr(0) = 0;
        Kokkos::deep_copy(d_rowptrU, h_rowptr);
      }

      if (lu) {
        auto d_rowptrL = Kokkos::subview(rowptrL, range_type(ptr, ptr+m+1));
        auto d_colindL = Kokkos::subview(colindL, range_type(nnzL, nnzL+s0.nnzL));
        auto d_nzvalsL = Kokkos::subview(nzvalsL, range_type(nnzL, nnzL+s0.nnzL));
        s0.colindL = d_colindL.data();
        s0.nzvalsL = d_nzvalsL.data();
        nnzL += s0.nnzL;

        // ========================
        // insert nonzeros
        extractor_crs.setRowPtr(s0.rowptrL);
        extractor_crs.setCrsView(s0.colindL, s0.nzvalsL);
        extractor_crs.setPivPtr(piv);
        {
          using team_policy_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                                      typename functor_type::ExtractValColTag>;
          team_policy_type team_policy((pend-pbeg)+1, Kokkos::AUTO());

          // >> launch functor to extract nonzero entries
          Kokkos::parallel_for("extract nzvals L", team_policy, extractor_crs);
          exec_space().fence();
        }
        // ========================
        // shift back
        // (TODO: shift first to avoid this)
        {
          //  copy to CPU, for now
          auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrL);
          for (ordinal_type i = m; i > 0 ; i--) h_rowptr(i) = h_rowptr(i-1);
          h_rowptr(0) = 0;
          Kokkos::deep_copy(d_rowptrL, h_rowptr);
        }
      } else if (s0.spmv_explicit_transpose) {
        // ========================
        // transpose
        // >> generate rowptr
        extractor_crs.setRowPtrT(s0.rowptrL);
        {
          // >> count nnz / row (transpose)
          using team_policy_type = Kokkos::RangePolicy<typename functor_type::TransPtrTag, exec_space>;
          team_policy_type team_policy(0, m);
          Kokkos::parallel_for("transpose pointer", team_policy, extractor_crs);
        }
        {
          // >> accumulate to generate rowptr (transpose)
          using range_policy_type = Kokkos::RangePolicy<exec_space>;
          Kokkos::parallel_scan("shiftRowptrT", range_policy_type(0, m+1), rowptr_sum(s0.rowptrL));
          exec_space().fence();
        }

        s0.nnzL = s0.nnzU;
        auto d_colindL = Kokkos::subview(colindL, range_type(nnzL, nnzL+s0.nnzL));
        auto d_nzvalsL = Kokkos::subview(nzvalsL, range_type(nnzL, nnzL+s0.nnzL));
        s0.colindL = d_colindL.data();
        s0.nzvalsL = d_nzvalsL.data();
        nnzL += s0.nnzL;

        // ========================
        // >> copy into transpose-matrix
        extractor_crs.setRowPtrT(s0.rowptrL);
        extractor_crs.setCrsViewT(s0.colindL, s0.nzvalsL);
        {
          using team_policy_type = Kokkos::RangePolicy<typename functor_type::TransMatTag, exec_space>;
          team_policy_type team_policy(0, m);
          Kokkos::parallel_for("transpose pointer", team_policy, extractor_crs);
          exec_space().fence();
        }
        // ========================
        // shift back
        // (TODO: shift first to avoid this)
        {
          // copy to CPU, for now
          auto d_rowptrL = Kokkos::subview(rowptrL, range_type(ptr, ptr+m+1));
          auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrL);
          for (ordinal_type i = m; i > 0 ; i--) h_rowptr(i) = h_rowptr(i-1);
          h_rowptr(0) = 0;
          Kokkos::deep_copy(d_rowptrL, h_rowptr);
        }
      }
      ptr += (1+m);
#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE)) || \
     defined(KOKKOS_ENABLE_HIP)
      const value_type one(1);
      const value_type zero(0);
      // ========================
      // create NVIDIA/AMD data structures for SpMV
      size_t buffer_size_L = 0;
      size_t buffer_size_U = 0;
      value_type alpha = one;
      value_type beta = zero;
#if defined(KOKKOS_ENABLE_CUDA)
      // create matrix
      cusparseCreateCsr(&s0.U_cusparse, m, m, s0.nnzU,
                        s0.rowptrU, s0.colindU, s0.nzvalsU,
                        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                        CUSPARSE_INDEX_BASE_ZERO, computeType);

#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
      cusparseSpMM_bufferSize(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                              &alpha, s0.U_cusparse, vecX, &beta, vecY,
                              computeType, TACHO_CUSPARSE_SPMM_ALG, &buffer_size_U);
#else
      cusparseSpMV_bufferSize(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, s0.U_cusparse, vecX, &beta, vecY,
                              computeType, TACHO_CUSPARSE_SPMV_ALG, &buffer_size_U);
#endif
      if (s0.spmv_explicit_transpose) {
        // create matrix (transpose(U) or L)
        cusparseCreateCsr(&s0.L_cusparse, m, m, s0.nnzL,
                          s0.rowptrL, s0.colindL, s0.nzvalsL,
                          CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                          CUSPARSE_INDEX_BASE_ZERO, computeType);
        // workspace size
#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
        cusparseSpMM_bufferSize(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, s0.L_cusparse, vecX, &beta, vecY,
                                computeType, TACHO_CUSPARSE_SPMM_ALG, &buffer_size_L);
#else
        cusparseSpMV_bufferSize(sparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, s0.L_cusparse, vecX, &beta, vecY,
                                computeType, TACHO_CUSPARSE_SPMV_ALG, &buffer_size_L);
#endif
      } else {
        // create matrix (L_cusparse stores the same ptrs as descrU, but optimized for trans)
        s0.nnzL = s0.nnzU;
        cusparseCreateCsr(&s0.L_cusparse, m, m, s0.nnzU,
                          s0.rowptrU, s0.colindU, s0.nzvalsU,
                          CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                          CUSPARSE_INDEX_BASE_ZERO, computeType);
        // workspace size for transpose SpMV
#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
        cusparseSpMM_bufferSize(sparseHandle, CUSPARSE_OPERATION_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, s0.L_cusparse, vecX, &beta, vecY,
                                computeType, TACHO_CUSPARSE_SPMM_ALG, &buffer_size_L);
#else
        cusparseSpMV_bufferSize(sparseHandle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, s0.L_cusparse, vecX, &beta, vecY,
                                computeType, TACHO_CUSPARSE_SPMV_ALG, &buffer_size_L);
#endif
      }
      // allocate workspace
      // allocate workspace
      if (buffer_size_U > buffer_U.extent(0)) {
        Kokkos::resize(buffer_U, buffer_size_U);
      }
      if (buffer_size_L > buffer_L.extent(0)) {
        Kokkos::resize(buffer_L, buffer_size_L);
      }
#elif defined(KOKKOS_ENABLE_HIP)
      // create matrix
      rocsparse_create_csr_descr(&(s0.descrU), m, m, s0.nnzU,
                                 s0.rowptrU, s0.colindU, s0.nzvalsU,
                                 rocsparse_indextype_i32, rocsparse_indextype_i32, rocsparse_index_base_zero, rocsparse_compute_type);
      // workspace
      tacho_rocsparse_spmv
          (sparseHandle, rocsparse_operation_none,
           &alpha, s0.descrU, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_buffer_size,
           #endif
           &buffer_size_U, nullptr);
      // allocate workspace
      if (buffer_size_U > buffer_U.extent(0)) {
        Kokkos::resize(buffer_U, buffer_size_U);
      }
      #if ROCM_VERSION >= 50400
      // preprocess
      buffer_size_U = buffer_U.extent(0);
      tacho_rocsparse_spmv
          (sparseHandle, rocsparse_operation_none,
           &alpha, s0.descrU, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           rocsparse_spmv_stage_preprocess,
           &buffer_size_U, (void*)buffer_U.data());
      #endif
      if (s0.spmv_explicit_transpose) {
        // create matrix (transpose)
        _status = rocsparse_create_csr_descr(&(s0.descrL), m, m, s0.nnzL,
                                             s0.rowptrL, s0.colindL, s0.nzvalsL,
                                             rocsparse_indextype_i32, rocsparse_indextype_i32,
                                             rocsparse_index_base_zero, rocsparse_compute_type);
        check_OnDeviceSPMV_Status("rocsparse_create_crs_descr", "rocsparse");
        // workspace
        _status = tacho_rocsparse_spmv
          (sparseHandle, rocsparse_operation_none,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_buffer_size,
           #endif
           &buffer_size_L, nullptr);
        check_OnDeviceSPMV_Status("rocsparse_create_spmv", "rocsparse");
        // allocate workspace
        if (buffer_size_L > buffer_L.extent(0)) {
          Kokkos::resize(buffer_L, buffer_size_L);
        }
        #if ROCM_VERSION >= 50400
        // preprocess
        buffer_size_L = buffer_L.extent(0);
        tacho_rocsparse_spmv
          (sparseHandle, rocsparse_operation_none,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
            rocsparse_spmv_stage_preprocess,
           &buffer_size_L, (void*)buffer_L.data());
        #endif
      } else {
        // create matrix, transpose (L_cusparse stores the same ptrs as descrU, but optimized for trans)
        _status = rocsparse_create_csr_descr(&(s0.descrL), m, m, s0.nnzU,
                                             s0.rowptrU, s0.colindU, s0.nzvalsU,
                                             rocsparse_indextype_i32, rocsparse_indextype_i32,
                                             rocsparse_index_base_zero, rocsparse_compute_type);
        check_OnDeviceSPMV_Status("rocsparse_create_crs_descr", "rocsparse");
        // workspace (transpose)
        _status = tacho_rocsparse_spmv
          (sparseHandle, rocsparse_operation_transpose,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_buffer_size,
           #endif
           &buffer_size_L, nullptr);
        check_OnDeviceSPMV_Status("rocsparse_create_spmv", "rocsparse");
        // allcate workspace
        if (buffer_size_L > buffer_L.extent(0)) {
          Kokkos::resize(buffer_L, buffer_size_L);
        }
        #if ROCM_VERSION >= 50400
        // preprocess
        buffer_size_L = buffer_L.extent(0);
        tacho_rocsparse_spmv
          (sparseHandle, rocsparse_operation_transpose,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
            rocsparse_spmv_stage_preprocess,
           &buffer_size_L, (void*)buffer_L.data());
        #endif
      }
#endif
#endif
    }
#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE)) || \
     defined(KOKKOS_ENABLE_HIP)
#if defined(KOKKOS_ENABLE_CUDA)
#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
    cusparseDestroyDnMat(vecX);
    cusparseDestroyDnMat(vecY);
#else
    cusparseDestroyDnVec(vecX);
    cusparseDestroyDnVec(vecY);
#endif
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_destroy_dnvec_descr(vecX);
    rocsparse_destroy_dnvec_descr(vecY);
#endif
#endif
    if (verbose) {
      printf("LevelSetTools:loadCRS(method = %d)\n",method);
      printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
      printf( " CRS with total nnzL = %ld and nnzU = %ld\n",colindL.extent(0),colindU.extent(0) );
      if (store_transpose) printf( "  > explicitly storing transpose\n" );
    }
    return r_val;
  }

 };
} // namespace Tacho
#endif
