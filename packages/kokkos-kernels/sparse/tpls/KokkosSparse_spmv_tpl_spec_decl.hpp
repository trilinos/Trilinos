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

#ifndef KOKKOSPARSE_SPMV_TPL_SPEC_DECL_HPP_
#define KOKKOSPARSE_SPMV_TPL_SPEC_DECL_HPP_

#include <sstream>

#include "KokkosKernels_Controls.hpp"

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#include "KokkosSparse_Utils_cusparse.hpp"

namespace KokkosSparse {
namespace Impl {

template <class AMatrix, class XVector, class YVector>
void spmv_cusparse(const Kokkos::Cuda& exec,
                   const KokkosKernels::Experimental::Controls& controls,
                   const char mode[],
                   typename YVector::non_const_value_type const& alpha,
                   const AMatrix& A, const XVector& x,
                   typename YVector::non_const_value_type const& beta,
                   const YVector& y) {
  using offset_type = typename AMatrix::non_const_size_type;
  using value_type  = typename AMatrix::non_const_value_type;

  /* initialize cusparse library */
  cusparseHandle_t cusparseHandle = controls.getCusparseHandle();
  /* Set cuSPARSE to use the given stream until this function exits */
  TemporarySetCusparseStream(cusparseHandle, exec);

  /* Set the operation mode */
  cusparseOperation_t myCusparseOperation;
  switch (toupper(mode[0])) {
    case 'N': myCusparseOperation = CUSPARSE_OPERATION_NON_TRANSPOSE; break;
    case 'T': myCusparseOperation = CUSPARSE_OPERATION_TRANSPOSE; break;
    case 'H':
      myCusparseOperation = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
      break;
    default: {
      std::ostringstream out;
      out << "Mode " << mode << " invalid for cuSPARSE SpMV.\n";
      throw std::invalid_argument(out.str());
    }
  }
  // cuSPARSE doesn't directly support mode H with real values, but this is
  // equivalent to mode T
  if (myCusparseOperation == CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE &&
      !Kokkos::ArithTraits<value_type>::isComplex)
    myCusparseOperation = CUSPARSE_OPERATION_TRANSPOSE;

#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)

  using entry_type = typename AMatrix::non_const_ordinal_type;
  /* Check that cusparse can handle the types of the input Kokkos::CrsMatrix */
  const cusparseIndexType_t myCusparseOffsetType =
      cusparse_index_type_t_from<offset_type>();
  const cusparseIndexType_t myCusparseEntryType =
      cusparse_index_type_t_from<entry_type>();

  cudaDataType myCudaDataType;
  if (std::is_same<value_type, float>::value)
    myCudaDataType = CUDA_R_32F;
  else if (std::is_same<value_type, double>::value)
    myCudaDataType = CUDA_R_64F;
  else if (std::is_same<value_type, Kokkos::complex<float>>::value)
    myCudaDataType = CUDA_C_32F;
  else if (std::is_same<value_type, Kokkos::complex<double>>::value)
    myCudaDataType = CUDA_C_64F;
  else
    throw std::logic_error(
        "Scalar (data) type of CrsMatrix isn't supported by cuSPARSE, yet TPL "
        "layer says it is");

  /* create matrix */
  cusparseSpMatDescr_t A_cusparse;
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(
      &A_cusparse, A.numRows(), A.numCols(), A.nnz(),
      (void*)A.graph.row_map.data(), (void*)A.graph.entries.data(),
      (void*)A.values.data(), myCusparseOffsetType, myCusparseEntryType,
      CUSPARSE_INDEX_BASE_ZERO, myCudaDataType));

  /* create lhs and rhs */
  cusparseDnVecDescr_t vecX, vecY;
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(
      &vecX, x.extent_int(0), (void*)x.data(), myCudaDataType));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(
      &vecY, y.extent_int(0), (void*)y.data(), myCudaDataType));

  size_t bufferSize = 0;
  void* dBuffer     = NULL;
#if CUSPARSE_VERSION >= 11301
  cusparseSpMVAlg_t alg = CUSPARSE_SPMV_ALG_DEFAULT;
#else
  cusparseSpMVAlg_t alg = CUSPARSE_MV_ALG_DEFAULT;
#endif
  if (controls.isParameter("algorithm")) {
    const std::string algName = controls.getParameter("algorithm");
    if (algName == "default")
#if CUSPARSE_VERSION >= 11301
      alg = CUSPARSE_SPMV_ALG_DEFAULT;
#else
      alg = CUSPARSE_MV_ALG_DEFAULT;
#endif
    else if (algName == "merge")
#if CUSPARSE_VERSION >= 11301
      alg = CUSPARSE_SPMV_CSR_ALG2;
#else
      alg = CUSPARSE_CSRMV_ALG2;
#endif
  }
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMV_bufferSize(
      cusparseHandle, myCusparseOperation, &alpha, A_cusparse, vecX, &beta,
      vecY, myCudaDataType, alg, &bufferSize));
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc(&dBuffer, bufferSize));

  /* perform SpMV */
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMV(cusparseHandle, myCusparseOperation,
                                         &alpha, A_cusparse, vecX, &beta, vecY,
                                         myCudaDataType, alg, dBuffer));

  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(dBuffer));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecX));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecY));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(A_cusparse));

#elif (9000 <= CUDA_VERSION)

  /* create and set the matrix descriptor */
  cusparseMatDescr_t descrA = 0;
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descrA));
  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
  KOKKOS_CUSPARSE_SAFE_CALL(
      cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));

  /* perform the actual SpMV operation */
  if (std::is_same<int, offset_type>::value) {
    if (std::is_same<value_type, float>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseScsrmv(
          cusparseHandle, myCusparseOperation, A.numRows(), A.numCols(),
          A.nnz(), reinterpret_cast<float const*>(&alpha), descrA,
          reinterpret_cast<float const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(),
          reinterpret_cast<float const*>(x.data()),
          reinterpret_cast<float const*>(&beta),
          reinterpret_cast<float*>(y.data())));

    } else if (std::is_same<value_type, double>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDcsrmv(
          cusparseHandle, myCusparseOperation, A.numRows(), A.numCols(),
          A.nnz(), reinterpret_cast<double const*>(&alpha), descrA,
          reinterpret_cast<double const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(),
          reinterpret_cast<double const*>(x.data()),
          reinterpret_cast<double const*>(&beta),
          reinterpret_cast<double*>(y.data())));
    } else if (std::is_same<value_type, Kokkos::complex<float>>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCcsrmv(
          cusparseHandle, myCusparseOperation, A.numRows(), A.numCols(),
          A.nnz(), reinterpret_cast<cuComplex const*>(&alpha), descrA,
          reinterpret_cast<cuComplex const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(),
          reinterpret_cast<cuComplex const*>(x.data()),
          reinterpret_cast<cuComplex const*>(&beta),
          reinterpret_cast<cuComplex*>(y.data())));
    } else if (std::is_same<value_type, Kokkos::complex<double>>::value) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseZcsrmv(
          cusparseHandle, myCusparseOperation, A.numRows(), A.numCols(),
          A.nnz(), reinterpret_cast<cuDoubleComplex const*>(&alpha), descrA,
          reinterpret_cast<cuDoubleComplex const*>(A.values.data()),
          A.graph.row_map.data(), A.graph.entries.data(),
          reinterpret_cast<cuDoubleComplex const*>(x.data()),
          reinterpret_cast<cuDoubleComplex const*>(&beta),
          reinterpret_cast<cuDoubleComplex*>(y.data())));
    } else {
      throw std::logic_error(
          "Trying to call cusparse SpMV with a scalar type not float/double, "
          "nor complex of either!");
    }
  } else {
    throw std::logic_error(
        "With cuSPARSE pre-10.0, offset type must be int. Something wrong with "
        "TPL avail logic.");
  }

  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyMatDescr(descrA));
#endif  // CUDA_VERSION
}

#define KOKKOSSPARSE_SPMV_CUSPARSE(SCALAR, ORDINAL, OFFSET, LAYOUT, SPACE,  \
                                   COMPILE_LIBRARY)                         \
  template <>                                                               \
  struct SPMV<                                                              \
      Kokkos::Cuda,                                                         \
      KokkosSparse::CrsMatrix<                                              \
          SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::Cuda, SPACE>, \
          Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const>,           \
      Kokkos::View<                                                         \
          SCALAR const*, LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>,       \
          Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,  \
      Kokkos::View<SCALAR*, LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>,    \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                \
      true, COMPILE_LIBRARY> {                                              \
    using device_type       = Kokkos::Device<Kokkos::Cuda, SPACE>;          \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;      \
    using AMatrix = CrsMatrix<SCALAR const, ORDINAL const, device_type,     \
                              memory_trait_type, OFFSET const>;             \
    using XVector = Kokkos::View<                                           \
        SCALAR const*, LAYOUT, device_type,                                 \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;    \
    using YVector =                                                         \
        Kokkos::View<SCALAR*, LAYOUT, device_type, memory_trait_type>;      \
    using Controls = KokkosKernels::Experimental::Controls;                 \
                                                                            \
    using coefficient_type = typename YVector::non_const_value_type;        \
                                                                            \
    static void spmv(const Kokkos::Cuda& exec, const Controls& controls,    \
                     const char mode[], const coefficient_type& alpha,      \
                     const AMatrix& A, const XVector& x,                    \
                     const coefficient_type& beta, const YVector& y) {      \
      std::string label = "KokkosSparse::spmv[TPL_CUSPARSE," +              \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";        \
      Kokkos::Profiling::pushRegion(label);                                 \
      spmv_cusparse(exec, controls, mode, alpha, A, x, beta, y);            \
      Kokkos::Profiling::popRegion();                                       \
    }                                                                       \
  };

// BMK: cuSPARSE that comes with CUDA 9 does not support tranpose or conjugate
// transpose modes. No version of cuSPARSE supports mode C (conjugate, non
// transpose). In those cases, fall back to KokkosKernels native spmv.

#if (9000 <= CUDA_VERSION)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutRight,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutRight,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int,
                           Kokkos::LayoutLeft, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int,
                           Kokkos::LayoutRight, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int,
                           Kokkos::LayoutRight, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int, int, Kokkos::LayoutRight,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int, int, Kokkos::LayoutRight,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int,
                           Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int, int,
                           Kokkos::LayoutRight, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int, int,
                           Kokkos::LayoutRight, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)

#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutLeft,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutRight,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutLeft,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutRight,
                           Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int64_t, size_t,
                           Kokkos::LayoutLeft, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int64_t, size_t,
                           Kokkos::LayoutRight, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int64_t, size_t,
                           Kokkos::LayoutLeft, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int64_t, size_t,
                           Kokkos::LayoutRight, Kokkos::CudaSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutLeft,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutRight,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutLeft,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutRight,
                           Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int64_t, size_t,
                           Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<double>, int64_t, size_t,
                           Kokkos::LayoutRight, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int64_t, size_t,
                           Kokkos::LayoutLeft, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_CUSPARSE(Kokkos::complex<float>, int64_t, size_t,
                           Kokkos::LayoutRight, Kokkos::CudaUVMSpace,
                           KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif  // defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)
#endif  // 9000 <= CUDA_VERSION

#undef KOKKOSSPARSE_SPMV_CUSPARSE

}  // namespace Impl
}  // namespace KokkosSparse
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// rocSPARSE
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)
#include "KokkosSparse_Utils_rocsparse.hpp"

namespace KokkosSparse {
namespace Impl {

template <class AMatrix, class XVector, class YVector>
void spmv_rocsparse(const Kokkos::HIP& exec,
                    const KokkosKernels::Experimental::Controls& controls,
                    const char mode[],
                    typename YVector::non_const_value_type const& alpha,
                    const AMatrix& A, const XVector& x,
                    typename YVector::non_const_value_type const& beta,
                    const YVector& y) {
  using offset_type = typename AMatrix::non_const_size_type;
  using entry_type  = typename AMatrix::non_const_ordinal_type;
  using value_type  = typename AMatrix::non_const_value_type;

  /* initialize rocsparse library */
  rocsparse_handle handle = controls.getRocsparseHandle();
  /* Set rocsparse to use the given stream until this function exits */
  TemporarySetRocsparseStream(handle, exec);

  /* Set the operation mode */
  rocsparse_operation myRocsparseOperation = mode_kk_to_rocsparse(mode);

  /* Set the index type */
  rocsparse_indextype offset_index_type = rocsparse_index_type<offset_type>();
  rocsparse_indextype entry_index_type  = rocsparse_index_type<entry_type>();

  /* Set the scalar type */
  rocsparse_datatype compute_type = rocsparse_compute_type<value_type>();

  /* Create the rocsparse mat and csr descr */
  rocsparse_mat_descr Amat;
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&Amat));
  rocsparse_spmat_descr Aspmat;
  // We need to do some casting to void*
  // Note that row_map is always a const view so const_cast is necessary,
  // however entries and values may not be const so we need to check first.
  void* csr_row_ptr =
      static_cast<void*>(const_cast<offset_type*>(A.graph.row_map.data()));
  void* csr_col_ind =
      static_cast<void*>(const_cast<entry_type*>(A.graph.entries.data()));
  void* csr_val = static_cast<void*>(const_cast<value_type*>(A.values.data()));

  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_csr_descr(
      &Aspmat, A.numRows(), A.numCols(), A.nnz(), csr_row_ptr, csr_col_ind,
      csr_val, offset_index_type, entry_index_type, rocsparse_index_base_zero,
      compute_type));

  /* Create rocsparse dense vectors for X and Y */
  rocsparse_dnvec_descr vecX, vecY;
  void* x_data = static_cast<void*>(
      const_cast<typename XVector::non_const_value_type*>(x.data()));
  void* y_data = static_cast<void*>(
      const_cast<typename YVector::non_const_value_type*>(y.data()));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_dnvec_descr(
      &vecX, x.extent_int(0), x_data,
      rocsparse_compute_type<typename XVector::non_const_value_type>()));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_dnvec_descr(
      &vecY, y.extent_int(0), y_data,
      rocsparse_compute_type<typename YVector::non_const_value_type>()));

  /* Actually perform the SpMV operation, first size buffer, then compute result
   */
  size_t buffer_size     = 0;
  void* tmp_buffer       = nullptr;
  rocsparse_spmv_alg alg = rocsparse_spmv_alg_default;
  // Note, Dec 6th 2021 - lbv:
  // rocSPARSE offers two diffrent algorithms for spmv
  // 1. ocsparse_spmv_alg_csr_adaptive
  // 2. rocsparse_spmv_alg_csr_stream
  // it is unclear which one is the default algorithm
  // or what both algorithms are intended for?
  if (controls.isParameter("algorithm")) {
    const std::string algName = controls.getParameter("algorithm");
    if (algName == "default")
      alg = rocsparse_spmv_alg_default;
    else if (algName == "merge")
      alg = rocsparse_spmv_alg_csr_stream;
  }

#if KOKKOSSPARSE_IMPL_ROCM_VERSION >= 60000
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(
      rocsparse_spmv(handle, myRocsparseOperation, &alpha, Aspmat, vecX, &beta,
                     vecY, compute_type, alg, rocsparse_spmv_stage_buffer_size,
                     &buffer_size, tmp_buffer));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(&tmp_buffer, buffer_size));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(
      rocsparse_spmv(handle, myRocsparseOperation, &alpha, Aspmat, vecX, &beta,
                     vecY, compute_type, alg, rocsparse_spmv_stage_compute,
                     &buffer_size, tmp_buffer));
#elif KOKKOSSPARSE_IMPL_ROCM_VERSION >= 50400
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_spmv_ex(
      handle, myRocsparseOperation, &alpha, Aspmat, vecX, &beta, vecY,
      compute_type, alg, rocsparse_spmv_stage_auto, &buffer_size, tmp_buffer));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(&tmp_buffer, buffer_size));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_spmv_ex(
      handle, myRocsparseOperation, &alpha, Aspmat, vecX, &beta, vecY,
      compute_type, alg, rocsparse_spmv_stage_auto, &buffer_size, tmp_buffer));
#else
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(
      rocsparse_spmv(handle, myRocsparseOperation, &alpha, Aspmat, vecX, &beta,
                     vecY, compute_type, alg, &buffer_size, tmp_buffer));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(&tmp_buffer, buffer_size));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(
      rocsparse_spmv(handle, myRocsparseOperation, &alpha, Aspmat, vecX, &beta,
                     vecY, compute_type, alg, &buffer_size, tmp_buffer));
#endif
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(tmp_buffer));

  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_dnvec_descr(vecY));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_dnvec_descr(vecX));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_spmat_descr(Aspmat));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_mat_descr(Amat));
}

#define KOKKOSSPARSE_SPMV_ROCSPARSE(SCALAR, LAYOUT, COMPILE_LIBRARY)          \
  template <>                                                                 \
  struct SPMV<                                                                \
      Kokkos::HIP,                                                            \
      KokkosSparse::CrsMatrix<SCALAR const, rocsparse_int const,              \
                              Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,  \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>,        \
                              rocsparse_int const>,                           \
      Kokkos::View<                                                           \
          SCALAR const*, LAYOUT,                                              \
          Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,                      \
          Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,    \
      Kokkos::View<SCALAR*, LAYOUT,                                           \
                   Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>,             \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                  \
      true, COMPILE_LIBRARY> {                                                \
    using device_type       = Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>;  \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;        \
    using AMatrix = CrsMatrix<SCALAR const, rocsparse_int const, device_type, \
                              memory_trait_type, rocsparse_int const>;        \
    using XVector = Kokkos::View<                                             \
        SCALAR const*, LAYOUT, device_type,                                   \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;      \
    using YVector =                                                           \
        Kokkos::View<SCALAR*, LAYOUT, device_type, memory_trait_type>;        \
    using Controls = KokkosKernels::Experimental::Controls;                   \
                                                                              \
    using coefficient_type = typename YVector::non_const_value_type;          \
                                                                              \
    static void spmv(const Kokkos::HIP& exec, const Controls& controls,       \
                     const char mode[], const coefficient_type& alpha,        \
                     const AMatrix& A, const XVector& x,                      \
                     const coefficient_type& beta, const YVector& y) {        \
      std::string label = "KokkosSparse::spmv[TPL_ROCSPARSE," +               \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";          \
      Kokkos::Profiling::pushRegion(label);                                   \
      spmv_rocsparse(exec, controls, mode, alpha, A, x, beta, y);             \
      Kokkos::Profiling::popRegion();                                         \
    }                                                                         \
  };

KOKKOSSPARSE_SPMV_ROCSPARSE(double, Kokkos::LayoutLeft,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ROCSPARSE(double, Kokkos::LayoutRight,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ROCSPARSE(float, Kokkos::LayoutLeft,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ROCSPARSE(float, Kokkos::LayoutRight,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<double>, Kokkos::LayoutLeft,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<double>, Kokkos::LayoutRight,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<float>, Kokkos::LayoutLeft,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ROCSPARSE(Kokkos::complex<float>, Kokkos::LayoutRight,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)

#undef KOKKOSSPARSE_SPMV_ROCSPARSE

}  // namespace Impl
}  // namespace KokkosSparse
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include <mkl.h>
#include "KokkosSparse_Utils_mkl.hpp"

namespace KokkosSparse {
namespace Impl {

#if (__INTEL_MKL__ > 2017)
// MKL 2018 and above: use new interface: sparse_matrix_t and mkl_sparse_?_mv()

inline void spmv_mkl(sparse_operation_t op, float alpha, float beta, MKL_INT m,
                     MKL_INT n, const MKL_INT* Arowptrs,
                     const MKL_INT* Aentries, const float* Avalues,
                     const float* x, float* y) {
  sparse_matrix_t A_mkl;
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_s_create_csr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, m, n, const_cast<MKL_INT*>(Arowptrs),
      const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries),
      const_cast<float*>(Avalues)));
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_s_mv(op, alpha, A_mkl, A_descr, x, beta, y));
}

inline void spmv_mkl(sparse_operation_t op, double alpha, double beta,
                     MKL_INT m, MKL_INT n, const MKL_INT* Arowptrs,
                     const MKL_INT* Aentries, const double* Avalues,
                     const double* x, double* y) {
  sparse_matrix_t A_mkl;
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_d_create_csr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, m, n, const_cast<MKL_INT*>(Arowptrs),
      const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries),
      const_cast<double*>(Avalues)));
  KOKKOSKERNELS_MKL_SAFE_CALL(
      mkl_sparse_d_mv(op, alpha, A_mkl, A_descr, x, beta, y));
}

inline void spmv_mkl(sparse_operation_t op, Kokkos::complex<float> alpha,
                     Kokkos::complex<float> beta, MKL_INT m, MKL_INT n,
                     const MKL_INT* Arowptrs, const MKL_INT* Aentries,
                     const Kokkos::complex<float>* Avalues,
                     const Kokkos::complex<float>* x,
                     Kokkos::complex<float>* y) {
  sparse_matrix_t A_mkl;
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_create_csr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, m, n, const_cast<MKL_INT*>(Arowptrs),
      const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries),
      (MKL_Complex8*)Avalues));
  MKL_Complex8 alpha_mkl{alpha.real(), alpha.imag()};
  MKL_Complex8 beta_mkl{beta.real(), beta.imag()};
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_c_mv(
      op, alpha_mkl, A_mkl, A_descr, reinterpret_cast<const MKL_Complex8*>(x),
      beta_mkl, reinterpret_cast<MKL_Complex8*>(y)));
}

inline void spmv_mkl(sparse_operation_t op, Kokkos::complex<double> alpha,
                     Kokkos::complex<double> beta, MKL_INT m, MKL_INT n,
                     const MKL_INT* Arowptrs, const MKL_INT* Aentries,
                     const Kokkos::complex<double>* Avalues,
                     const Kokkos::complex<double>* x,
                     Kokkos::complex<double>* y) {
  sparse_matrix_t A_mkl;
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_create_csr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, m, n, const_cast<MKL_INT*>(Arowptrs),
      const_cast<MKL_INT*>(Arowptrs + 1), const_cast<MKL_INT*>(Aentries),
      (MKL_Complex16*)Avalues));
  MKL_Complex16 alpha_mkl{alpha.real(), alpha.imag()};
  MKL_Complex16 beta_mkl{beta.real(), beta.imag()};
  KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_z_mv(
      op, alpha_mkl, A_mkl, A_descr, reinterpret_cast<const MKL_Complex16*>(x),
      beta_mkl, reinterpret_cast<MKL_Complex16*>(y)));
}

// Note: classic MKL runs on Serial/OpenMP but can't use our execution space
// instances
#define KOKKOSSPARSE_SPMV_MKL(SCALAR, EXECSPACE, COMPILE_LIBRARY)              \
  template <>                                                                  \
  struct SPMV<EXECSPACE,                                                       \
              KokkosSparse::CrsMatrix<                                         \
                  SCALAR const, MKL_INT const,                                 \
                  Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,                \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const>,     \
              Kokkos::View<SCALAR const*, Kokkos::LayoutLeft,                  \
                           Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,       \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged |            \
                                                Kokkos::RandomAccess>>,        \
              Kokkos::View<SCALAR*, Kokkos::LayoutLeft,                        \
                           Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,       \
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>>,           \
              true, COMPILE_LIBRARY> {                                         \
    using device_type = Kokkos::Device<EXECSPACE, Kokkos::HostSpace>;          \
    using AMatrix =                                                            \
        CrsMatrix<SCALAR const, MKL_INT const, device_type,                    \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, MKL_INT const>;     \
    using XVector = Kokkos::View<                                              \
        SCALAR const*, Kokkos::LayoutLeft, device_type,                        \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector = Kokkos::View<SCALAR*, Kokkos::LayoutLeft, device_type,     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using coefficient_type = typename YVector::non_const_value_type;           \
    using Controls         = KokkosKernels::Experimental::Controls;            \
                                                                               \
    static void spmv(const EXECSPACE&, const Controls&, const char mode[],     \
                     const coefficient_type& alpha, const AMatrix& A,          \
                     const XVector& x, const coefficient_type& beta,           \
                     const YVector& y) {                                       \
      std::string label = "KokkosSparse::spmv[TPL_MKL," +                      \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";           \
      Kokkos::Profiling::pushRegion(label);                                    \
      spmv_mkl(mode_kk_to_mkl(mode[0]), alpha, beta, A.numRows(), A.numCols(), \
               A.graph.row_map.data(), A.graph.entries.data(),                 \
               A.values.data(), x.data(), y.data());                           \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
  };

#ifdef KOKKOS_ENABLE_SERIAL
KOKKOSSPARSE_SPMV_MKL(float, Kokkos::Serial, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(double, Kokkos::Serial,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<float>, Kokkos::Serial,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<double>, Kokkos::Serial,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif

#ifdef KOKKOS_ENABLE_OPENMP
KOKKOSSPARSE_SPMV_MKL(float, Kokkos::OpenMP, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(double, Kokkos::OpenMP,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<float>, Kokkos::OpenMP,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_MKL(Kokkos::complex<double>, Kokkos::OpenMP,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif

#undef KOKKOSSPARSE_SPMV_MKL
#endif

#if defined(KOKKOS_ENABLE_SYCL) && \
    !defined(KOKKOSKERNELS_ENABLE_TPL_MKL_SYCL_OVERRIDE)
inline oneapi::mkl::transpose mode_kk_to_onemkl(char mode_kk) {
  switch (toupper(mode_kk)) {
    case 'N': return oneapi::mkl::transpose::nontrans;
    case 'T': return oneapi::mkl::transpose::trans;
    case 'H': return oneapi::mkl::transpose::conjtrans;
    default:;
  }
  throw std::invalid_argument(
      "Invalid mode for oneMKL (should be one of N, T, H)");
}

template <bool is_complex>
struct spmv_onemkl_wrapper {};

template <>
struct spmv_onemkl_wrapper<false> {
  template <class execution_space, class matrix_type, class xview_type,
            class yview_type>
  static void spmv(const execution_space& exec, oneapi::mkl::transpose mkl_mode,
                   typename matrix_type::non_const_value_type const alpha,
                   const matrix_type& A, const xview_type& x,
                   typename matrix_type::non_const_value_type const beta,
                   const yview_type& y) {
    using scalar_type  = typename matrix_type::non_const_value_type;
    using ordinal_type = typename matrix_type::non_const_ordinal_type;

    // oneAPI doesn't directly support mode H with real values, but this is
    // equivalent to mode T
    if (mkl_mode == oneapi::mkl::transpose::conjtrans &&
        !Kokkos::ArithTraits<scalar_type>::isComplex)
      mkl_mode = oneapi::mkl::transpose::trans;

    oneapi::mkl::sparse::matrix_handle_t handle = nullptr;
    oneapi::mkl::sparse::init_matrix_handle(&handle);
    auto ev_set = oneapi::mkl::sparse::set_csr_data(
        exec.sycl_queue(), handle, A.numRows(), A.numCols(),
        oneapi::mkl::index_base::zero,
        const_cast<ordinal_type*>(A.graph.row_map.data()),
        const_cast<ordinal_type*>(A.graph.entries.data()),
        const_cast<scalar_type*>(A.values.data()));
    auto ev_opt = oneapi::mkl::sparse::optimize_gemv(
        exec.sycl_queue(), mkl_mode, handle, {ev_set});
    auto ev_gemv =
        oneapi::mkl::sparse::gemv(exec.sycl_queue(), mkl_mode, alpha, handle,
                                  x.data(), beta, y.data(), {ev_opt});
    auto ev_release = oneapi::mkl::sparse::release_matrix_handle(
        exec.sycl_queue(), &handle, {ev_gemv});
    ev_release.wait();
  }
};

template <>
struct spmv_onemkl_wrapper<true> {
  template <class execution_space, class matrix_type, class xview_type,
            class yview_type>
  static void spmv(const execution_space& exec, oneapi::mkl::transpose mkl_mode,
                   typename matrix_type::non_const_value_type const alpha,
                   const matrix_type& A, const xview_type& x,
                   typename matrix_type::non_const_value_type const beta,
                   const yview_type& y) {
    using scalar_type  = typename matrix_type::non_const_value_type;
    using ordinal_type = typename matrix_type::non_const_ordinal_type;
    using mag_type     = typename Kokkos::ArithTraits<scalar_type>::mag_type;

    oneapi::mkl::sparse::matrix_handle_t handle = nullptr;
    oneapi::mkl::sparse::init_matrix_handle(&handle);
    auto ev_set = oneapi::mkl::sparse::set_csr_data(
        exec.sycl_queue(), handle, static_cast<ordinal_type>(A.numRows()),
        static_cast<ordinal_type>(A.numCols()), oneapi::mkl::index_base::zero,
        const_cast<ordinal_type*>(A.graph.row_map.data()),
        const_cast<ordinal_type*>(A.graph.entries.data()),
        reinterpret_cast<std::complex<mag_type>*>(
            const_cast<scalar_type*>(A.values.data())));
    auto ev_opt = oneapi::mkl::sparse::optimize_gemv(
        exec.sycl_queue(), mkl_mode, handle, {ev_set});
    auto ev_gemv = oneapi::mkl::sparse::gemv(
        exec.sycl_queue(), mkl_mode, alpha, handle,
        reinterpret_cast<std::complex<mag_type>*>(
            const_cast<scalar_type*>(x.data())),
        beta, reinterpret_cast<std::complex<mag_type>*>(y.data()), {ev_opt});
    auto ev_release = oneapi::mkl::sparse::release_matrix_handle(
        exec.sycl_queue(), &handle, {ev_gemv});
    ev_release.wait();
  }
};

#define KOKKOSSPARSE_SPMV_ONEMKL(SCALAR, ORDINAL, MEMSPACE, COMPILE_LIBRARY) \
  template <>                                                                \
  struct SPMV<                                                               \
      Kokkos::Experimental::SYCL,                                            \
      KokkosSparse::CrsMatrix<                                               \
          SCALAR const, ORDINAL const,                                       \
          Kokkos::Device<Kokkos::Experimental::SYCL, MEMSPACE>,              \
          Kokkos::MemoryTraits<Kokkos::Unmanaged>, ORDINAL const>,           \
      Kokkos::View<                                                          \
          SCALAR const*, Kokkos::LayoutLeft,                                 \
          Kokkos::Device<Kokkos::Experimental::SYCL, MEMSPACE>,              \
          Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,   \
      Kokkos::View<SCALAR*, Kokkos::LayoutLeft,                              \
                   Kokkos::Device<Kokkos::Experimental::SYCL, MEMSPACE>,     \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                 \
      true, COMPILE_LIBRARY> {                                               \
    using execution_space = Kokkos::Experimental::SYCL;                      \
    using device_type     = Kokkos::Device<execution_space, MEMSPACE>;       \
    using AMatrix =                                                          \
        CrsMatrix<SCALAR const, ORDINAL const, device_type,                  \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, ORDINAL const>;   \
    using XVector = Kokkos::View<                                            \
        SCALAR const*, Kokkos::LayoutLeft, device_type,                      \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;     \
    using YVector = Kokkos::View<SCALAR*, Kokkos::LayoutLeft, device_type,   \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;   \
    using coefficient_type = typename YVector::non_const_value_type;         \
    using Controls         = KokkosKernels::Experimental::Controls;          \
                                                                             \
    static void spmv(const execution_space& exec, const Controls&,           \
                     const char mode[], const coefficient_type& alpha,       \
                     const AMatrix& A, const XVector& x,                     \
                     const coefficient_type& beta, const YVector& y) {       \
      std::string label = "KokkosSparse::spmv[TPL_ONEMKL," +                 \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";         \
      Kokkos::Profiling::pushRegion(label);                                  \
      oneapi::mkl::transpose mkl_mode = mode_kk_to_onemkl(mode[0]);          \
      spmv_onemkl_wrapper<Kokkos::ArithTraits<SCALAR>::is_complex>::spmv(    \
          exec, mkl_mode, alpha, A, x, beta, y);                             \
      Kokkos::Profiling::popRegion();                                        \
    }                                                                        \
  };

KOKKOSSPARSE_SPMV_ONEMKL(float, std::int32_t,
                         Kokkos::Experimental::SYCLDeviceUSMSpace,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ONEMKL(double, std::int32_t,
                         Kokkos::Experimental::SYCLDeviceUSMSpace,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ONEMKL(Kokkos::complex<float>, std::int32_t,
                         Kokkos::Experimental::SYCLDeviceUSMSpace,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ONEMKL(Kokkos::complex<double>, std::int32_t,
                         Kokkos::Experimental::SYCLDeviceUSMSpace,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)

KOKKOSSPARSE_SPMV_ONEMKL(float, std::int64_t,
                         Kokkos::Experimental::SYCLDeviceUSMSpace,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ONEMKL(double, std::int64_t,
                         Kokkos::Experimental::SYCLDeviceUSMSpace,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ONEMKL(Kokkos::complex<float>, std::int64_t,
                         Kokkos::Experimental::SYCLDeviceUSMSpace,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
KOKKOSSPARSE_SPMV_ONEMKL(Kokkos::complex<double>, std::int64_t,
                         Kokkos::Experimental::SYCLDeviceUSMSpace,
                         KOKKOSKERNELS_IMPL_COMPILE_LIBRARY)
#endif
}  // namespace Impl
}  // namespace KokkosSparse
#endif

#endif  // KOKKOSPARSE_SPMV_TPL_SPEC_DECL_HPP_
