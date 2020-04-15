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

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#include "KokkosKernels_SparseUtils_cusparse.hpp"

namespace KokkosSparse {
namespace Impl {

  template <class AMatrix, class XVector, class YVector>
  void spmv_cusparse(const char mode[],
		     typename YVector::non_const_value_type const & alpha,
		     const AMatrix& A,
		     const XVector& x,
		     typename YVector::non_const_value_type const & beta,
		     const YVector& y) {
    using offset_type  = typename AMatrix::non_const_size_type;
    using value_type   = typename AMatrix::non_const_value_type;

    /* initialize cusparse library */
    cusparseHandle_t cusparseHandle = 0;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreate(&cusparseHandle));

    /* Set the operation mode */
    cusparseOperation_t myCusparseOperation = CUSPARSE_OPERATION_NON_TRANSPOSE;
    if(mode[0] == Transpose[0]) {myCusparseOperation = CUSPARSE_OPERATION_TRANSPOSE;}

#if defined(CUSPARSE_VERSION) && (10300 <= CUSPARSE_VERSION)

    /* Check that cusparse can handle the types of the input Kokkos::CrsMatrix */
    cusparseIndexType_t myCusparseIndexType;
    if(std::is_same<offset_type, int>::value)     {myCusparseIndexType = CUSPARSE_INDEX_32I;}
    if(std::is_same<offset_type, int64_t>::value) {myCusparseIndexType = CUSPARSE_INDEX_64I;}
    cudaDataType myCudaDataType;
    if(std::is_same<value_type, float>::value)  {myCudaDataType = CUDA_R_32F;}
    if(std::is_same<value_type, double>::value) {myCudaDataType = CUDA_R_64F;}

    /* create matrix */
    cusparseSpMatDescr_t A_cusparse;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateCsr(&A_cusparse, A.numRows(), A.numCols(), A.nnz(),
						const_cast<offset_type*>(A.graph.row_map.data()),
						const_cast<offset_type*>(A.graph.entries.data()),
						const_cast<value_type*>(A.values.data()),
						myCusparseIndexType,
						myCusparseIndexType,
						CUSPARSE_INDEX_BASE_ZERO,
						myCudaDataType));

    /* create lhs and rhs */
    cusparseDnVecDescr_t vecX, vecY;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&vecX, x.extent_int(0), const_cast<value_type*>(x.data()), myCudaDataType));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&vecY, y.extent_int(0), const_cast<value_type*>(y.data()), myCudaDataType));

    size_t bufferSize = 0;
    void*  dBuffer    = NULL;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMV_bufferSize(cusparseHandle, myCusparseOperation,
						      &alpha, A_cusparse, vecX, &beta, vecY, myCudaDataType,
						      CUSPARSE_CSRMV_ALG1, &bufferSize));
    CUDA_SAFE_CALL(cudaMalloc(&dBuffer, bufferSize));

    /* perform SpMV */
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpMV(cusparseHandle, myCusparseOperation,
					   &alpha, A_cusparse, vecX, &beta, vecY, myCudaDataType,
					   CUSPARSE_CSRMV_ALG1, dBuffer));

    CUDA_SAFE_CALL(cudaFree(dBuffer));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecX));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecY));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(A_cusparse));
#else

    /* create and set the matrix descriptor */
    cusparseMatDescr_t descrA = 0;
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descrA));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
    KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO));

    /* perform the actual SpMV operation */
    if(std::is_same<int, offset_type>::value) {
      if (std::is_same<value_type,float>::value) {
	KOKKOS_CUSPARSE_SAFE_CALL(cusparseScsrmv(cusparseHandle, myCusparseOperation,
						 A.numRows(), A.numCols(), A.nnz(),
						 reinterpret_cast<const float *>(&alpha), descrA,
						 reinterpret_cast<const float *>(A.values.data()),
						 A.graph.row_map.data(), A.graph.entries.data(),
						 reinterpret_cast<const float *>(x.data()),
						 reinterpret_cast<const float *>(&beta),
						 reinterpret_cast<float *>(y.data()) ));

      } else  if (std::is_same<value_type,double>::value) {
	KOKKOS_CUSPARSE_SAFE_CALL(cusparseDcsrmv(cusparseHandle, myCusparseOperation,
						 A.numRows(), A.numCols(), A.nnz(),
						 reinterpret_cast<double const *>(&alpha), descrA,
						 reinterpret_cast<double const *>(A.values.data()),
						 A.graph.row_map.data(), A.graph.entries.data(),
						 reinterpret_cast<double const *>(x.data()),
						 reinterpret_cast<double const *>(&beta),
						 reinterpret_cast<double *>(y.data()) ));
      } else {
	throw std::logic_error("Trying to call cusparse SpMV with a scalar type that is not float or double!");
      }
    } else {
      throw std::logic_error("Trying to call cusparse SpMV with an offset type that is not int!");
    }

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyMatDescr(descrA));
#endif // CUSPARSE_VERSION

    KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroy(cusparseHandle));
    cusparseHandle = 0;
  }

#define KOKKOSSPARSE_SPMV_CUSPARSE(SCALAR, OFFSET, LAYOUT, COMPILE_LIBRARY) \
  template<>								\
  struct SPMV<SCALAR const,  OFFSET const, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>, Kokkos::MemoryTraits<Kokkos::Unmanaged>, OFFSET const, \
	      SCALAR const*, LAYOUT,       Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>, Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
	      SCALAR*,       LAYOUT,       Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>, Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
	      true, COMPILE_LIBRARY> {					\
    using device_type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>; \
    using memory_trait_type = Kokkos::MemoryTraits<Kokkos::Unmanaged>;	\
    using AMatrix = CrsMatrix<SCALAR const, OFFSET const, device_type, memory_trait_type, OFFSET const>; \
    using XVector = Kokkos::View<SCALAR const*, LAYOUT,device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>>; \
    using YVector = Kokkos::View<SCALAR*, LAYOUT, device_type, memory_trait_type>; \
									\
    using coefficient_type = typename YVector::non_const_value_type;	\
									\
    static void spmv (const char mode[],				\
		      const coefficient_type& alpha,			\
		      const AMatrix& A,					\
		      const XVector& x,					\
		      const coefficient_type& beta,			\
		      const YVector& y) {				\
      std::string label = "KokkosSparse::spmv[TPL_CUSPARSE," + Kokkos::ArithTraits<SCALAR>::name() + "]"; \
      Kokkos::Profiling::pushRegion(label);				\
      spmv_cusparse(mode, alpha, A, x, beta, y);			\
      Kokkos::Profiling::popRegion();					\
    }									\
  };

  KOKKOSSPARSE_SPMV_CUSPARSE(double, int, Kokkos::LayoutLeft,  true)
  KOKKOSSPARSE_SPMV_CUSPARSE(double, int, Kokkos::LayoutLeft,  false)
  KOKKOSSPARSE_SPMV_CUSPARSE(double, int, Kokkos::LayoutRight, true)
  KOKKOSSPARSE_SPMV_CUSPARSE(double, int, Kokkos::LayoutRight, false)
  KOKKOSSPARSE_SPMV_CUSPARSE(float,  int, Kokkos::LayoutLeft,  true)
  KOKKOSSPARSE_SPMV_CUSPARSE(float,  int, Kokkos::LayoutLeft,  false)
  KOKKOSSPARSE_SPMV_CUSPARSE(float,  int, Kokkos::LayoutRight, true)
  KOKKOSSPARSE_SPMV_CUSPARSE(float,  int, Kokkos::LayoutRight, false)

#undef KOKKOSSPARSE_SPMV_CUSPARSE

} // namespace Impl
} // namespace KokkosSparse
#endif // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#endif // KOKKOSPARSE_SPMV_TPL_SPEC_DECL_HPP_
