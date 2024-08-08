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

#ifndef KOKKOSPARSE_SPMV_MV_TPL_SPEC_DECL_HPP_
#define KOKKOSPARSE_SPMV_MV_TPL_SPEC_DECL_HPP_

#include <sstream>
#include "KokkosKernels_tpl_handles_decl.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

/* CUSPARSE_VERSION < 10301 either doesn't have cusparseSpMM
   or the non-tranpose version produces incorrect results.

   Version 11702 corresponds to CUDA 11.6.1, which also produces incorrect
   results. 11701 (CUDA 11.6.0) is OK.
*/
#if defined(CUSPARSE_VERSION) && (10301 <= CUSPARSE_VERSION) && (CUSPARSE_VERSION != 11702)
#include "cusparse.h"
#include "KokkosSparse_Utils_cusparse.hpp"

namespace KokkosSparse {
namespace Impl {

/* Derive a compute type for various operand types.
   cusparseSpMM does not always allow the same compute type as operand types
   This should be consistent with the allowed operand types for cusparseSpMM,
   as needed for TPL availability. Current definition does not comprehensively
   cover all cusparseSpMM options.

   cuSparse 11.5.1+ does not support uniform precision for FP16
   Otherwise, uniform precision is supported
*/
template <typename AScalar, typename XScalar = AScalar, typename YScalar = AScalar>
cudaDataType compute_type() {
  return cuda_data_type_from<AScalar>();
}
#if CUSPARSE_VERSION >= 11501
template <>
inline cudaDataType compute_type<Kokkos::Experimental::half_t>() {
  return CUDA_R_32F;
}
#else
template <>
inline cudaDataType compute_type<Kokkos::Experimental::half_t>() {
  return cuda_data_type_from<Kokkos::Experimental::half_t>();
}
#endif

/*! \brief convert a 2D view to a cusparseDnMatDescr_t

*/
template <typename ViewType, std::enable_if_t<ViewType::rank == 2, bool> = true>
cusparseDnMatDescr_t make_cusparse_dn_mat_descr_t(ViewType &view) {
  // If the view is LayoutRight, we still need to create descr as column-major
  // but it should be an implicit transpose, meaning dimensions and strides are
  // swapped
  bool transpose    = std::is_same_v<typename ViewType::array_layout, Kokkos::LayoutRight>;
  const size_t rows = transpose ? view.extent(1) : view.extent(0);
  const size_t cols = transpose ? view.extent(0) : view.extent(1);
  const size_t ld   = transpose ? view.stride(0) : view.stride(1);

  // cusparseCreateCsr notes it is safe to const_cast this away for input
  // pointers to a descriptor as long as that descriptor is not an output
  // parameter
  void *values = const_cast<typename ViewType::non_const_value_type *>(view.data());

  cudaDataType valueType = cuda_data_type_from<typename ViewType::non_const_value_type>();

  // col-major is the only supported order in 10301
  // ignore the layout of the provided view, and expect the caller to
  // fix with a transpose operation, if possible.
  // This should be revisited once cusparse supports row-major dense matrices
  const cusparseOrder_t order = CUSPARSE_ORDER_COL;

  cusparseDnMatDescr_t descr;
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateDnMat(&descr, static_cast<int64_t>(rows), static_cast<int64_t>(cols),
                                                static_cast<int64_t>(ld), values, valueType, order));

  return descr;
}

template <class Handle, class AMatrix, class XVector, class YVector>
void spmv_mv_cusparse(const Kokkos::Cuda &exec, Handle *handle, const char mode[],
                      typename YVector::const_value_type &alpha, const AMatrix &A, const XVector &x,
                      typename YVector::const_value_type &beta, const YVector &y) {
  static_assert(XVector::rank == 2, "should only be instantiated for multivector");
  static_assert(YVector::rank == 2, "should only be instantiated for multivector");

  using offset_type  = typename AMatrix::non_const_size_type;
  using entry_type   = typename AMatrix::non_const_ordinal_type;
  using value_type   = typename AMatrix::non_const_value_type;
  using x_value_type = typename XVector::non_const_value_type;
  using y_value_type = typename YVector::non_const_value_type;

  /* initialize cusparse library */
  cusparseHandle_t cusparseHandle = KokkosKernels::Impl::CusparseSingleton::singleton().cusparseHandle;
  /* Set cuSPARSE to use the given stream until this function exits */
  TemporarySetCusparseStream tscs(cusparseHandle, exec);

  /* Check that cusparse can handle the types of the input Kokkos::CrsMatrix */
  const cusparseIndexType_t myCusparseOffsetType = cusparse_index_type_t_from<offset_type>();
  const cusparseIndexType_t myCusparseEntryType  = cusparse_index_type_t_from<entry_type>();
  const cudaDataType aCusparseType               = cuda_data_type_from<value_type>();

  /* Set the operation mode */
  cusparseOperation_t opA;
  switch (toupper(mode[0])) {
    case 'N': opA = CUSPARSE_OPERATION_NON_TRANSPOSE; break;
    case 'T': opA = CUSPARSE_OPERATION_TRANSPOSE; break;
    case 'H': opA = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE; break;
    default: {
      std::ostringstream out;
      out << "Mode " << mode << " invalid for cuSPARSE SpMV MV.\n";
      throw std::invalid_argument(out.str());
    }
  }

  /* create lhs and rhs
     NOTE: The descriptions always say vecX and vecY are column-major cusparse
     order. For CUSPARSE_VERSION 10301 this is the only supported ordering. if X
     is not LayoutLeft, we can fix with a transpose. If cusparseSpMM ever
     supports row-major dense matrices, this logic will have to be reworked */
  constexpr bool xIsLL = std::is_same<typename XVector::array_layout, Kokkos::LayoutLeft>::value;
  constexpr bool xIsLR = std::is_same<typename XVector::array_layout, Kokkos::LayoutRight>::value;
  static_assert(xIsLL || xIsLR, "X multivector was not LL or LR (TPL error)");
  static_assert(std::is_same_v<typename YVector::array_layout, Kokkos::LayoutLeft>,
                "Y multivector was not LL (TPL error)");
  cusparseDnMatDescr_t vecX = make_cusparse_dn_mat_descr_t(x);
  cusparseDnMatDescr_t vecY = make_cusparse_dn_mat_descr_t(y);
  cusparseOperation_t opB   = xIsLL ? CUSPARSE_OPERATION_NON_TRANSPOSE : CUSPARSE_OPERATION_TRANSPOSE;

// CUSPARSE_MM_ALG_DEFAULT was deprecated in CUDA 11.0.1 / cuSPARSE 11.0.0 and
// removed in CUDA 12.0.0 / cuSPARSE 12.0.0
#if CUSPARSE_VERSION < 11000
  cusparseSpMMAlg_t algo = CUSPARSE_MM_ALG_DEFAULT;
#else
  cusparseSpMMAlg_t algo = CUSPARSE_SPMM_ALG_DEFAULT;
#endif

  // the precision of the SpMV
  const cudaDataType computeType = compute_type<value_type, x_value_type, y_value_type>();

  // cuSPARSE fails when conjugate_transpose is requested on R types
  // to avoid this problem we switch to transpose since the two are
  // equivalent in that case.
  if ((computeType == CUDA_R_32F) || (computeType == CUDA_R_64F)) {
    if (opA == CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE) {
      opA = CUSPARSE_OPERATION_TRANSPOSE;
    }
    if (opB == CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE) {
      opB = CUSPARSE_OPERATION_TRANSPOSE;
    }
  }

  KokkosSparse::Impl::CuSparse10_SpMV_Data *subhandle;
  if (handle->tpl_rank2) {
    subhandle = dynamic_cast<KokkosSparse::Impl::CuSparse10_SpMV_Data *>(handle->tpl_rank2);
    if (!subhandle) throw std::runtime_error("KokkosSparse::spmv: subhandle is not set up for cusparse");
    subhandle->set_exec_space(exec);
  } else {
    subhandle         = new KokkosSparse::Impl::CuSparse10_SpMV_Data(exec);
    handle->tpl_rank2 = subhandle;
    /* create matrix */
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&subhandle->mat, A.numRows(), A.numCols(), A.nnz(),
                                                (void *)A.graph.row_map.data(), (void *)A.graph.entries.data(),
                                                (void *)A.values.data(), myCusparseOffsetType, myCusparseEntryType,
                                                CUSPARSE_INDEX_BASE_ZERO, aCusparseType));

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMM_bufferSize(cusparseHandle, opA, opB, &alpha, subhandle->mat, vecX, &beta,
                                                      vecY, computeType, algo, &subhandle->bufferSize));

    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc(&subhandle->buffer, subhandle->bufferSize));
  }

  KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMM(cusparseHandle, opA, opB, &alpha, subhandle->mat, vecX, &beta, vecY,
                                         computeType, algo, subhandle->buffer));

  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnMat(vecX));
  KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnMat(vecY));
}

#define KOKKOSSPARSE_SPMV_MV_CUSPARSE(SCALAR, ORDINAL, OFFSET, XL, YL, SPACE)                                       \
  template <>                                                                                                       \
  struct SPMV_MV<                                                                                                   \
      Kokkos::Cuda, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, SPACE, SCALAR, OFFSET, ORDINAL>,               \
      KokkosSparse::CrsMatrix<SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::Cuda, SPACE>,                     \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const>,                               \
      Kokkos::View<SCALAR const **, XL, Kokkos::Device<Kokkos::Cuda, SPACE>,                                        \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                 \
      Kokkos::View<SCALAR **, YL, Kokkos::Device<Kokkos::Cuda, SPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,    \
      false, true> {                                                                                                \
    using device_type       = Kokkos::Device<Kokkos::Cuda, SPACE>;                                                  \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;                                              \
    using Handle            = KokkosSparse::Impl::SPMVHandleImpl<Kokkos::Cuda, SPACE, SCALAR, OFFSET, ORDINAL>;     \
    using AMatrix           = CrsMatrix<SCALAR const, ORDINAL const, device_type, memory_trait_type, OFFSET const>; \
    using XVector           = Kokkos::View<SCALAR const **, XL, device_type,                                        \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;         \
    using YVector           = Kokkos::View<SCALAR **, YL, device_type, memory_trait_type>;                          \
                                                                                                                    \
    using coefficient_type = typename YVector::non_const_value_type;                                                \
                                                                                                                    \
    static void spmv_mv(const Kokkos::Cuda &exec, Handle *handle, const char mode[], const coefficient_type &alpha, \
                        const AMatrix &A, const XVector &x, const coefficient_type &beta, const YVector &y) {       \
      std::string label = "KokkosSparse::spmv_mv[TPL_CUSPARSE," + Kokkos::ArithTraits<SCALAR>::name() + "]";        \
      Kokkos::Profiling::pushRegion(label);                                                                         \
      spmv_mv_cusparse(exec, handle, mode, alpha, A, x, beta, y);                                                   \
      Kokkos::Profiling::popRegion();                                                                               \
    }                                                                                                               \
  };

/* cusparseSpMM with following restrictions
 column-major ordering for Y
 col-major or row-major for X (see note below)
 32-bit indices for matrix A */
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                              Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                              Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                              Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                              Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                              Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                              Kokkos::CudaUVMSpace)

KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                              Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                              Kokkos::CudaUVMSpace)

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::Experimental::half_t, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                              Kokkos::CudaSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::Experimental::half_t, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                              Kokkos::CudaSpace)

KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::Experimental::half_t, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                              Kokkos::CudaUVMSpace)
KOKKOSSPARSE_SPMV_MV_CUSPARSE(Kokkos::Experimental::half_t, int, int, Kokkos::LayoutRight, Kokkos::LayoutLeft,
                              Kokkos::CudaUVMSpace)

#endif

#undef KOKKOSSPARSE_SPMV_MV_CUSPARSE

}  // namespace Impl
}  // namespace KokkosSparse
#endif  // defined(CUSPARSE_VERSION) && (10301 <= CUSPARSE_VERSION)
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// rocSPARSE
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)
#include "KokkosSparse_Utils_rocsparse.hpp"

namespace KokkosSparse {
namespace Impl {

template <class Handle, class AMatrix, class XVector, class YVector>
void spmv_mv_rocsparse(const Kokkos::HIP &exec, Handle *handle, const char mode[],
                       typename YVector::const_value_type &alpha, const AMatrix &A, const XVector &x,
                       typename YVector::const_value_type &beta, const YVector &y) {
  using offset_type = typename AMatrix::non_const_size_type;
  using entry_type  = typename AMatrix::non_const_ordinal_type;
  using value_type  = typename AMatrix::non_const_value_type;

  // initialize rocsparse library
  rocsparse_handle rocsparseHandle = KokkosKernels::Impl::RocsparseSingleton::singleton().rocsparseHandle;
  // Set rocsparse to use the given stream until this function exits
  TemporarySetRocsparseStream tsrs(rocsparseHandle, exec);

  rocsparse_operation rocsparseOperation = mode_kk_to_rocsparse(mode);
  rocsparse_indextype offset_index_type  = rocsparse_index_type<offset_type>();
  rocsparse_indextype entry_index_type   = rocsparse_index_type<entry_type>();
  rocsparse_datatype compute_type        = rocsparse_compute_type<value_type>();

  // Create rocsparse dense multivectors for X and Y
  void *x_data = static_cast<void *>(const_cast<typename XVector::non_const_value_type *>(x.data()));
  void *y_data = static_cast<void *>(const_cast<typename YVector::non_const_value_type *>(y.data()));

  size_t x_ld, y_ld;
  rocsparse_order x_order, y_order;
  if constexpr (std::is_same_v<typename XVector::array_layout, Kokkos::LayoutLeft>) {
    x_ld    = x.stride(1);
    x_order = rocsparse_order_column;
  } else {
    static_assert(std::is_same_v<typename XVector::array_layout, Kokkos::LayoutRight>,
                  "rocsparse_spmm internal logic error: x is neither LayoutLeft nor "
                  "LayoutRight");
    x_ld    = x.stride(0);
    x_order = rocsparse_order_row;
  }
  if constexpr (std::is_same_v<typename YVector::array_layout, Kokkos::LayoutLeft>) {
    y_ld    = y.stride(1);
    y_order = rocsparse_order_column;
  } else {
    static_assert(std::is_same_v<typename YVector::array_layout, Kokkos::LayoutRight>,
                  "rocsparse_spmm internal logic error: y is neither LayoutLeft nor "
                  "LayoutRight");
    y_ld    = y.stride(0);
    y_order = rocsparse_order_row;
  }

  rocsparse_dnmat_descr vecX, vecY;
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(
      rocsparse_create_dnmat_descr(&vecX, x.extent(0), x.extent(1), x_ld, x_data,
                                   rocsparse_compute_type<typename XVector::non_const_value_type>(), x_order));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(
      rocsparse_create_dnmat_descr(&vecY, y.extent(0), y.extent(1), y_ld, y_data,
                                   rocsparse_compute_type<typename YVector::non_const_value_type>(), y_order));

  rocsparse_spmm_alg alg = rocsparse_spmm_alg_default;

  KokkosSparse::Impl::RocSparse_CRS_SpMV_Data *subhandle;
  if (handle->tpl_rank2) {
    subhandle = dynamic_cast<KokkosSparse::Impl::RocSparse_CRS_SpMV_Data *>(handle->tpl_rank2);
    if (!subhandle)
      throw std::runtime_error(
          "KokkosSparse::spmv: rank-2 subhandle is not set up for rocsparse "
          "CRS");
    subhandle->set_exec_space(exec);
  } else {
    subhandle         = new KokkosSparse::Impl::RocSparse_CRS_SpMV_Data(exec);
    handle->tpl_rank2 = subhandle;
    // Create the rocsparse csr descr
    // We need to do some casting to void*
    // Note that row_map is always a const view so const_cast is necessary,
    // however entries and values may not be const so we need to check first.
    void *csr_row_ptr = static_cast<void *>(const_cast<offset_type *>(A.graph.row_map.data()));
    void *csr_col_ind = static_cast<void *>(const_cast<entry_type *>(A.graph.entries.data()));
    void *csr_val     = static_cast<void *>(const_cast<value_type *>(A.values.data()));

    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_csr_descr(
        &subhandle->mat, A.numRows(), A.numCols(), A.nnz(), csr_row_ptr, csr_col_ind, csr_val, offset_index_type,
        entry_index_type, rocsparse_index_base_zero, compute_type));

    // Size and allocate buffer, and analyze the matrix

    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_spmm(rocsparseHandle, rocsparseOperation, rocsparse_operation_none,
                                                   &alpha, subhandle->mat, vecX, &beta, vecY, compute_type, alg,
                                                   rocsparse_spmm_stage_buffer_size, &subhandle->bufferSize, nullptr));

    KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(&subhandle->buffer, subhandle->bufferSize));

    KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_spmm(
        rocsparseHandle, rocsparseOperation, rocsparse_operation_none, &alpha, subhandle->mat, vecX, &beta, vecY,
        compute_type, alg, rocsparse_spmm_stage_preprocess, &subhandle->bufferSize, subhandle->buffer));
  }

  // Perform the actual computation
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(
      rocsparse_spmm(rocsparseHandle, rocsparseOperation, rocsparse_operation_none, &alpha, subhandle->mat, vecX, &beta,
                     vecY, compute_type, alg, rocsparse_spmm_stage_compute, &subhandle->bufferSize, subhandle->buffer));

  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_dnmat_descr(vecY));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_dnmat_descr(vecX));
}

#define KOKKOSSPARSE_SPMV_MV_ROCSPARSE(SCALAR, XL, YL, MEMSPACE)                                                       \
  template <>                                                                                                          \
  struct SPMV_MV<                                                                                                      \
      Kokkos::HIP, KokkosSparse::Impl::SPMVHandleImpl<Kokkos::HIP, MEMSPACE, SCALAR, rocsparse_int, rocsparse_int>,    \
      KokkosSparse::CrsMatrix<SCALAR const, rocsparse_int const, Kokkos::Device<Kokkos::HIP, MEMSPACE>,                \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, rocsparse_int const>,                           \
      Kokkos::View<SCALAR const **, XL, Kokkos::Device<Kokkos::HIP, MEMSPACE>,                                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>,                                    \
      Kokkos::View<SCALAR **, YL, Kokkos::Device<Kokkos::HIP, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,     \
      false, true> {                                                                                                   \
    using device_type       = Kokkos::Device<Kokkos::HIP, MEMSPACE>;                                                   \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;                                                 \
    using Handle  = KokkosSparse::Impl::SPMVHandleImpl<Kokkos::HIP, MEMSPACE, SCALAR, rocsparse_int, rocsparse_int>;   \
    using AMatrix = CrsMatrix<SCALAR const, rocsparse_int const, device_type, memory_trait_type, rocsparse_int const>; \
    using XVector = Kokkos::View<SCALAR const **, XL, device_type,                                                     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;                      \
    using YVector = Kokkos::View<SCALAR **, YL, device_type, memory_trait_type>;                                       \
                                                                                                                       \
    using coefficient_type = typename YVector::non_const_value_type;                                                   \
                                                                                                                       \
    static void spmv_mv(const Kokkos::HIP &exec, Handle *handle, const char mode[], const coefficient_type &alpha,     \
                        const AMatrix &A, const XVector &x, const coefficient_type &beta, const YVector &y) {          \
      std::string label = "KokkosSparse::spmv_mv[TPL_ROCSPARSE," + Kokkos::ArithTraits<SCALAR>::name() + "]";          \
      Kokkos::Profiling::pushRegion(label);                                                                            \
      spmv_mv_rocsparse(exec, handle, mode, alpha, A, x, beta, y);                                                     \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define INST_ROCSPARSE_SCALAR_MEMSPACE(SCALAR, MEMSPACE)                                    \
  KOKKOSSPARSE_SPMV_MV_ROCSPARSE(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutLeft, MEMSPACE)  \
  KOKKOSSPARSE_SPMV_MV_ROCSPARSE(SCALAR, Kokkos::LayoutLeft, Kokkos::LayoutRight, MEMSPACE) \
  KOKKOSSPARSE_SPMV_MV_ROCSPARSE(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutLeft, MEMSPACE) \
  KOKKOSSPARSE_SPMV_MV_ROCSPARSE(SCALAR, Kokkos::LayoutRight, Kokkos::LayoutRight, MEMSPACE)

#define INST_ROCSPARSE_SCALAR(SCALAR)                      \
  INST_ROCSPARSE_SCALAR_MEMSPACE(SCALAR, Kokkos::HIPSpace) \
  INST_ROCSPARSE_SCALAR_MEMSPACE(SCALAR, Kokkos::HIPManagedSpace)

INST_ROCSPARSE_SCALAR(float)
INST_ROCSPARSE_SCALAR(double)
INST_ROCSPARSE_SCALAR(Kokkos::complex<float>)
INST_ROCSPARSE_SCALAR(Kokkos::complex<double>)

#undef INST_ROCSPARSE_SCALAR_MEMSPACE
#undef INST_ROCSPARSE_SCALAR
#undef KOKKOSSPARSE_SPMV_MV_ROCSPARSE

}  // namespace Impl
}  // namespace KokkosSparse
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#endif  // KOKKOSPARSE_SPMV_MV_TPL_SPEC_DECL_HPP_
