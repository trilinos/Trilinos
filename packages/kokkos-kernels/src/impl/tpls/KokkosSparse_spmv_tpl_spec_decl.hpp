/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSPARSE_SPMV_TPL_SPEC_DECL_HPP_
#define KOKKOSPARSE_SPMV_TPL_SPEC_DECL_HPP_

#include "KokkosKernels_Controls.hpp"

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#include "KokkosKernels_SparseUtils_cusparse.hpp"

namespace KokkosSparse {
namespace Impl {

template <class AMatrix, class XVector, class YVector>
void spmv_cusparse(const KokkosKernels::Experimental::Controls& controls,
                   const char mode[],
                   typename YVector::non_const_value_type const& alpha,
                   const AMatrix& A, const XVector& x,
                   typename YVector::non_const_value_type const& beta,
                   const YVector& y) {
  using offset_type = typename AMatrix::non_const_size_type;
  using entry_type  = typename AMatrix::non_const_ordinal_type;
  using value_type  = typename AMatrix::non_const_value_type;

  /* initialize cusparse library */
  cusparseHandle_t cusparseHandle = controls.getCusparseHandle();

  /* Set the operation mode */
  cusparseOperation_t myCusparseOperation;
  switch (toupper(mode[0])) {
    case 'N': myCusparseOperation = CUSPARSE_OPERATION_NON_TRANSPOSE; break;
    case 'T': myCusparseOperation = CUSPARSE_OPERATION_TRANSPOSE; break;
    case 'H':
      myCusparseOperation = CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE;
      break;
    default: {
      std::cerr << "Mode " << mode << " invalid for cuSPARSE SpMV.\n";
      throw std::invalid_argument("Invalid mode");
    }
  }

#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)

  /* Check that cusparse can handle the types of the input Kokkos::CrsMatrix */
  cusparseIndexType_t myCusparseOffsetType;
  if (std::is_same<offset_type, int>::value)
    myCusparseOffsetType = CUSPARSE_INDEX_32I;
  else if (std::is_same<offset_type, int64_t>::value ||
           std::is_same<offset_type, size_t>::value)
    myCusparseOffsetType = CUSPARSE_INDEX_64I;
  else
    throw std::logic_error(
        "Offset type of CrsMatrix isn't supported by cuSPARSE, yet TPL layer "
        "says it is");
  cusparseIndexType_t myCusparseEntryType;
  if (std::is_same<entry_type, int>::value)
    myCusparseEntryType = CUSPARSE_INDEX_32I;
  else if (std::is_same<entry_type, int64_t>::value)
    myCusparseEntryType = CUSPARSE_INDEX_64I;
  else
    throw std::logic_error(
        "Ordinal (entry) type of CrsMatrix isn't supported by cuSPARSE, yet "
        "TPL layer says it is");
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

  size_t bufferSize     = 0;
  void* dBuffer         = NULL;
  cusparseSpMVAlg_t alg = CUSPARSE_MV_ALG_DEFAULT;
  if (controls.isParameter("algorithm")) {
    const std::string algName = controls.getParameter("algorithm");
    if (algName == "default")
      alg = CUSPARSE_MV_ALG_DEFAULT;
    else if (algName == "merge")
      alg = CUSPARSE_CSRMV_ALG2;
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

#define KOKKOSSPARSE_SPMV_CUSPARSE(SCALAR, ORDINAL, OFFSET, LAYOUT, SPACE,     \
                                   COMPILE_LIBRARY)                            \
  template <>                                                                  \
  struct SPMV<                                                                 \
      SCALAR const, ORDINAL const, Kokkos::Device<Kokkos::Cuda, SPACE>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const, SCALAR const*,    \
      LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>,                             \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, SCALAR*, \
      LAYOUT, Kokkos::Device<Kokkos::Cuda, SPACE>,                             \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, COMPILE_LIBRARY> {        \
    using device_type       = Kokkos::Device<Kokkos::Cuda, SPACE>;             \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;         \
    using AMatrix = CrsMatrix<SCALAR const, ORDINAL const, device_type,        \
                              memory_trait_type, OFFSET const>;                \
    using XVector = Kokkos::View<                                              \
        SCALAR const*, LAYOUT, device_type,                                    \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector =                                                            \
        Kokkos::View<SCALAR*, LAYOUT, device_type, memory_trait_type>;         \
    using Controls = KokkosKernels::Experimental::Controls;                    \
                                                                               \
    using coefficient_type = typename YVector::non_const_value_type;           \
                                                                               \
    static void spmv(const Controls& controls, const char mode[],              \
                     const coefficient_type& alpha, const AMatrix& A,          \
                     const XVector& x, const coefficient_type& beta,           \
                     const YVector& y) {                                       \
      std::string label = "KokkosSparse::spmv[TPL_CUSPARSE," +                 \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";           \
      Kokkos::Profiling::pushRegion(label);                                    \
      spmv_cusparse(controls, mode, alpha, A, x, beta, y);                     \
      Kokkos::Profiling::popRegion();                                          \
    }                                                                          \
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
#endif
#endif

#undef KOKKOSSPARSE_SPMV_CUSPARSE

}  // namespace Impl
}  // namespace KokkosSparse
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// rocSPARSE
#if defined(KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE)
#include <rocsparse.h>
#include "KokkosKernels_SparseUtils_rocsparse.hpp"

namespace KokkosSparse {
namespace Impl {

template <class AMatrix, class XVector, class YVector>
void spmv_rocsparse(const KokkosKernels::Experimental::Controls& controls,
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
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(
      rocsparse_spmv(handle, myRocsparseOperation, &alpha, Aspmat, vecX, &beta,
                     vecY, compute_type, alg, &buffer_size, tmp_buffer));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(&tmp_buffer, buffer_size));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(
      rocsparse_spmv(handle, myRocsparseOperation, &alpha, Aspmat, vecX, &beta,
                     vecY, compute_type, alg, &buffer_size, tmp_buffer));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(tmp_buffer));

  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_dnvec_descr(vecY));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_dnvec_descr(vecX));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_spmat_descr(Aspmat));
  KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_destroy_mat_descr(Amat));
}

#define KOKKOSSPARSE_SPMV_ROCSPARSE(SCALAR, LAYOUT, COMPILE_LIBRARY)          \
  template <>                                                                 \
  struct SPMV<SCALAR const, rocsparse_int const,                              \
              Kokkos::Device<Kokkos::Experimental::HIP,                       \
                             Kokkos::Experimental::HIPSpace>,                 \
              Kokkos::MemoryTraits<Kokkos::Unmanaged>, rocsparse_int const,   \
              SCALAR const*, LAYOUT,                                          \
              Kokkos::Device<Kokkos::Experimental::HIP,                       \
                             Kokkos::Experimental::HIPSpace>,                 \
              Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, \
              SCALAR*, LAYOUT,                                                \
              Kokkos::Device<Kokkos::Experimental::HIP,                       \
                             Kokkos::Experimental::HIPSpace>,                 \
              Kokkos::MemoryTraits<Kokkos::Unmanaged>, true,                  \
              COMPILE_LIBRARY> {                                              \
    using device_type       = Kokkos::Device<Kokkos::Experimental::HIP,       \
                                       Kokkos::Experimental::HIPSpace>; \
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
    static void spmv(const Controls& controls, const char mode[],             \
                     const coefficient_type& alpha, const AMatrix& A,         \
                     const XVector& x, const coefficient_type& beta,          \
                     const YVector& y) {                                      \
      std::string label = "KokkosSparse::spmv[TPL_ROCSPARSE," +               \
                          Kokkos::ArithTraits<SCALAR>::name() + "]";          \
      Kokkos::Profiling::pushRegion(label);                                   \
      spmv_rocsparse(controls, mode, alpha, A, x, beta, y);                   \
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

}  // namespace Impl
}  // namespace KokkosSparse
#endif  // KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include <mkl.h>

namespace KokkosSparse {
namespace Impl {

#if (__INTEL_MKL__ > 2017)
// MKL 2018 and above: use new interface: sparse_matrix_t and mkl_sparse_?_mv()

// Note 12/03/21 - lbv:
// mkl_safe_call and mode_kk_to_mkl should
// be moved to some sparse or mkl utility
// header. It is likely that these will be
// reused for other kernels.
inline void mkl_safe_call(int errcode) {
  if (errcode != SPARSE_STATUS_SUCCESS)
    throw std::runtime_error("MKL returned non-success error code");
}

inline sparse_operation_t mode_kk_to_mkl(char mode_kk) {
  switch (toupper(mode_kk)) {
    case 'N': return SPARSE_OPERATION_NON_TRANSPOSE;
    case 'T': return SPARSE_OPERATION_TRANSPOSE;
    case 'H': return SPARSE_OPERATION_CONJUGATE_TRANSPOSE;
    default:;
  }
  throw std::invalid_argument(
      "Invalid mode for MKL (should be one of N, T, H)");
}

inline void spmv_mkl(sparse_operation_t op, float alpha, float beta, int m,
                     int n, const int* Arowptrs, const int* Aentries,
                     const float* Avalues, const float* x, float* y) {
  sparse_matrix_t A_mkl;
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  mkl_safe_call(mkl_sparse_s_create_csr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, m, n, const_cast<int*>(Arowptrs),
      const_cast<int*>(Arowptrs + 1), const_cast<int*>(Aentries),
      const_cast<float*>(Avalues)));
  mkl_safe_call(mkl_sparse_s_mv(op, alpha, A_mkl, A_descr, x, beta, y));
}

inline void spmv_mkl(sparse_operation_t op, double alpha, double beta, int m,
                     int n, const int* Arowptrs, const int* Aentries,
                     const double* Avalues, const double* x, double* y) {
  sparse_matrix_t A_mkl;
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  mkl_safe_call(mkl_sparse_d_create_csr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, m, n, const_cast<int*>(Arowptrs),
      const_cast<int*>(Arowptrs + 1), const_cast<int*>(Aentries),
      const_cast<double*>(Avalues)));
  mkl_safe_call(mkl_sparse_d_mv(op, alpha, A_mkl, A_descr, x, beta, y));
}

inline void spmv_mkl(sparse_operation_t op, Kokkos::complex<float> alpha,
                     Kokkos::complex<float> beta, int m, int n,
                     const int* Arowptrs, const int* Aentries,
                     const Kokkos::complex<float>* Avalues,
                     const Kokkos::complex<float>* x,
                     Kokkos::complex<float>* y) {
  sparse_matrix_t A_mkl;
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  mkl_safe_call(mkl_sparse_c_create_csr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, m, n, const_cast<int*>(Arowptrs),
      const_cast<int*>(Arowptrs + 1), const_cast<int*>(Aentries),
      (MKL_Complex8*)Avalues));
  MKL_Complex8& alpha_mkl = reinterpret_cast<MKL_Complex8&>(alpha);
  MKL_Complex8& beta_mkl  = reinterpret_cast<MKL_Complex8&>(beta);
  mkl_safe_call(mkl_sparse_c_mv(op, alpha_mkl, A_mkl, A_descr,
                                reinterpret_cast<const MKL_Complex8*>(x),
                                beta_mkl, reinterpret_cast<MKL_Complex8*>(y)));
}

inline void spmv_mkl(sparse_operation_t op, Kokkos::complex<double> alpha,
                     Kokkos::complex<double> beta, int m, int n,
                     const int* Arowptrs, const int* Aentries,
                     const Kokkos::complex<double>* Avalues,
                     const Kokkos::complex<double>* x,
                     Kokkos::complex<double>* y) {
  sparse_matrix_t A_mkl;
  matrix_descr A_descr;
  A_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
  A_descr.mode = SPARSE_FILL_MODE_FULL;
  A_descr.diag = SPARSE_DIAG_NON_UNIT;
  mkl_safe_call(mkl_sparse_z_create_csr(
      &A_mkl, SPARSE_INDEX_BASE_ZERO, m, n, const_cast<int*>(Arowptrs),
      const_cast<int*>(Arowptrs + 1), const_cast<int*>(Aentries),
      (MKL_Complex16*)Avalues));
  MKL_Complex16& alpha_mkl = reinterpret_cast<MKL_Complex16&>(alpha);
  MKL_Complex16& beta_mkl  = reinterpret_cast<MKL_Complex16&>(beta);
  mkl_safe_call(mkl_sparse_z_mv(op, alpha_mkl, A_mkl, A_descr,
                                reinterpret_cast<const MKL_Complex16*>(x),
                                beta_mkl, reinterpret_cast<MKL_Complex16*>(y)));
}

#define KOKKOSSPARSE_SPMV_MKL(SCALAR, EXECSPACE, COMPILE_LIBRARY)              \
  template <>                                                                  \
  struct SPMV<                                                                 \
      SCALAR const, int const, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,   \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, int const, SCALAR const*,       \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, SCALAR*, \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, COMPILE_LIBRARY> {        \
    using device_type = Kokkos::Device<EXECSPACE, Kokkos::HostSpace>;          \
    using AMatrix =                                                            \
        CrsMatrix<SCALAR const, int const, device_type,                        \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, int const>;         \
    using XVector = Kokkos::View<                                              \
        SCALAR const*, Kokkos::LayoutLeft, device_type,                        \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector = Kokkos::View<SCALAR*, Kokkos::LayoutLeft, device_type,     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using coefficient_type = typename YVector::non_const_value_type;           \
    using Controls         = KokkosKernels::Experimental::Controls;            \
                                                                               \
    static void spmv(const Controls&, const char mode[],                       \
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
#endif

#if (__INTEL_MKL__ == 2017)
// MKL 2017: use old interface: mkl_?csrmv
inline char mode_kk_to_mkl(char mode_kk) {
  switch (toupper(mode_kk)) {
    case 'N': return 'N';
    case 'T': return 'T';
    case 'H': return 'C';
    default:;
  }
  throw std::invalid_argument(
      "Invalid mode for MKL (should be one of N, T, H)");
}

inline void spmv_mkl(char mode, float alpha, float beta, int m, int n,
                     const int* Arowptrs, const int* Aentries,
                     const float* Avalues, const float* x, float* y) {
  mkl_scsrmv(&mode, &m, &n, &alpha, "G**C", Avalues, Aentries, Arowptrs,
             Arowptrs + 1, x, &beta, y);
}

inline void spmv_mkl(char mode, double alpha, double beta, int m, int n,
                     const int* Arowptrs, const int* Aentries,
                     const double* Avalues, const double* x, double* y) {
  mkl_dcsrmv(&mode, &m, &n, &alpha, "G**C", Avalues, Aentries, Arowptrs,
             Arowptrs + 1, x, &beta, y);
}

inline void spmv_mkl(char mode, Kokkos::complex<float> alpha,
                     Kokkos::complex<float> beta, int m, int n,
                     const int* Arowptrs, const int* Aentries,
                     const Kokkos::complex<float>* Avalues,
                     const Kokkos::complex<float>* x,
                     Kokkos::complex<float>* y) {
  const MKL_Complex8* alpha_mkl = reinterpret_cast<const MKL_Complex8*>(&alpha);
  const MKL_Complex8* beta_mkl  = reinterpret_cast<const MKL_Complex8*>(&beta);
  const MKL_Complex8* Avalues_mkl =
      reinterpret_cast<const MKL_Complex8*>(Avalues);
  const MKL_Complex8* x_mkl = reinterpret_cast<const MKL_Complex8*>(x);
  MKL_Complex8* y_mkl       = reinterpret_cast<MKL_Complex8*>(y);
  mkl_ccsrmv(&mode, &m, &n, alpha_mkl, "G**C", Avalues_mkl, Aentries, Arowptrs,
             Arowptrs + 1, x_mkl, beta_mkl, y_mkl);
}

inline void spmv_mkl(char mode, Kokkos::complex<double> alpha,
                     Kokkos::complex<double> beta, int m, int n,
                     const int* Arowptrs, const int* Aentries,
                     const Kokkos::complex<double>* Avalues,
                     const Kokkos::complex<double>* x,
                     Kokkos::complex<double>* y) {
  const MKL_Complex16* alpha_mkl =
      reinterpret_cast<const MKL_Complex16*>(&alpha);
  const MKL_Complex16* beta_mkl = reinterpret_cast<const MKL_Complex16*>(&beta);
  const MKL_Complex16* Avalues_mkl =
      reinterpret_cast<const MKL_Complex16*>(Avalues);
  const MKL_Complex16* x_mkl = reinterpret_cast<const MKL_Complex16*>(x);
  MKL_Complex16* y_mkl       = reinterpret_cast<MKL_Complex16*>(y);
  mkl_zcsrmv(&mode, &m, &n, alpha_mkl, "G**C", Avalues_mkl, Aentries, Arowptrs,
             Arowptrs + 1, x_mkl, beta_mkl, y_mkl);
}

#define KOKKOSSPARSE_SPMV_MKL(SCALAR, EXECSPACE, COMPILE_LIBRARY)              \
  template <>                                                                  \
  struct SPMV<                                                                 \
      SCALAR const, int const, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,   \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, int const, SCALAR const*,       \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, SCALAR*, \
      Kokkos::LayoutLeft, Kokkos::Device<EXECSPACE, Kokkos::HostSpace>,        \
      Kokkos::MemoryTraits<Kokkos::Unmanaged>, true, COMPILE_LIBRARY> {        \
    using device_type = Kokkos::Device<EXECSPACE, Kokkos::HostSpace>;          \
    using AMatrix =                                                            \
        CrsMatrix<SCALAR const, int const, device_type,                        \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, int const>;         \
    using XVector = Kokkos::View<                                              \
        SCALAR const*, Kokkos::LayoutLeft, device_type,                        \
        Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;       \
    using YVector = Kokkos::View<SCALAR*, Kokkos::LayoutLeft, device_type,     \
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;     \
    using coefficient_type = typename YVector::non_const_value_type;           \
    using Controls         = KokkosKernels::Experimental::Controls;            \
                                                                               \
    static void spmv(const Controls&, const char mode[],                       \
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
#endif

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
}  // namespace Impl
}  // namespace KokkosSparse
#endif

#endif  // KOKKOSPARSE_SPMV_TPL_SPEC_DECL_HPP_
