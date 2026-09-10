// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_CUSOLVER_FUNCTIONMAP_HPP
#define AMESOS2_CUSOLVER_FUNCTIONMAP_HPP

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_cuSOLVER_TypeMap.hpp"

#include <cuda.h>
#include <cusolverSp.h>
#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>

#ifdef HAVE_TEUCHOS_COMPLEX
#include <cuComplex.h>
#endif

namespace Amesos2 {

  template <>
  struct FunctionMap<cuSOLVER,double>
  {
    static cusolverStatus_t bufferInfo(
                 cusolverDnHandle_t handle,
                 int n,
                 double * A,
                 int lda,
                 int * lwork)
    {
      return cusolverDnDgetrf_bufferSize(handle, n, n, A, lda, lwork);
    }

    static cusolverStatus_t sparseBufferInfo(
                 cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 cusparseMatDescr_t & desc,
                 const double * values,
                 const int * rowPtr,
                 const int * colIdx,
                 csrcholInfo_t & chol_info,
                 size_t * internalDataInBytes,
                 size_t * workspaceInBytes)
    {
      return cusolverSpDcsrcholBufferInfo(handle, size, nnz, desc, values,
        rowPtr, colIdx, chol_info, internalDataInBytes, workspaceInBytes);
    }

    static cusolverStatus_t numeric(
                 cusolverDnHandle_t handle,
                 int n,
                 double * A,
                 int lda,
                 double * work,
                 int * ipiv,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnDgetrf(
        handle, n, n, A, lda, work, ipiv, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cusolverStatus_t sparseNumeric(
                 cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 cusparseMatDescr_t & desc,
                 const double * values,
                 const int * rowPtr,
                 const int * colIdx,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      cusolverStatus_t status = cusolverSpDcsrcholFactor(
        handle, size, nnz, desc, values, rowPtr, colIdx, chol_info, buffer);
      cudaDeviceSynchronize();
      return status;
    }

    static cusolverStatus_t solveLU(
                 cusolverDnHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const double * LU,
                 int ldlu,
                 const int * ipiv,
                 double * B,
                 int ldb,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnDgetrs(
        handle, trans, n, nrhs, LU, ldlu, ipiv, B, ldb, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    // Solve X = op(A^{-1}) * B via dense matrix-matrix product.
    static cublasStatus_t solve(
                 cublasHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const double * inverse,
                 int ldinv,
                 const double * B,
                 int ldb,
                 double * X,
                 int ldx)
    {
      const double one = 1.0, zero = 0.0;
      return cublasDgemm(handle, trans, CUBLAS_OP_N,
        n, nrhs, n, &one, inverse, ldinv, B, ldb, &zero, X, ldx);
    }

    static cusolverStatus_t sparseSolve(
                 cusolverSpHandle_t handle,
                 int size,
                 const double * b,
                 double * x,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      cusolverStatus_t status = cusolverSpDcsrcholSolve(
        handle, size, b, x, chol_info, buffer);
      cudaDeviceSynchronize();
      return status;
    }
  };

  template <>
  struct FunctionMap<cuSOLVER,float>
  {
    static cusolverStatus_t bufferInfo(
                 cusolverDnHandle_t handle,
                 int n,
                 float * A,
                 int lda,
                 int * lwork)
    {
      return cusolverDnSgetrf_bufferSize(handle, n, n, A, lda, lwork);
    }

    static cusolverStatus_t sparseBufferInfo(
                 cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 cusparseMatDescr_t & desc,
                 const float * values,
                 const int * rowPtr,
                 const int * colIdx,
                 csrcholInfo_t & chol_info,
                 size_t * internalDataInBytes,
                 size_t * workspaceInBytes)
    {
      return cusolverSpScsrcholBufferInfo(handle, size, nnz, desc, values,
        rowPtr, colIdx, chol_info, internalDataInBytes, workspaceInBytes);
    }

    static cusolverStatus_t numeric(
                 cusolverDnHandle_t handle,
                 int n,
                 float * A,
                 int lda,
                 float * work,
                 int * ipiv,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnSgetrf(
        handle, n, n, A, lda, work, ipiv, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cusolverStatus_t sparseNumeric(
                 cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 cusparseMatDescr_t & desc,
                 const float * values,
                 const int * rowPtr,
                 const int * colIdx,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      cusolverStatus_t status = cusolverSpScsrcholFactor(
        handle, size, nnz, desc, values, rowPtr, colIdx, chol_info, buffer);
      cudaDeviceSynchronize();
      return status;
    }

    static cusolverStatus_t solveLU(
                 cusolverDnHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const float * LU,
                 int ldlu,
                 const int * ipiv,
                 float * B,
                 int ldb,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnSgetrs(
        handle, trans, n, nrhs, LU, ldlu, ipiv, B, ldb, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cublasStatus_t solve(
                 cublasHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const float * inverse,
                 int ldinv,
                 const float * B,
                 int ldb,
                 float * X,
                 int ldx)
    {
      const float one = 1.0f, zero = 0.0f;
      return cublasSgemm(handle, trans, CUBLAS_OP_N,
        n, nrhs, n, &one, inverse, ldinv, B, ldb, &zero, X, ldx);
    }

    static cusolverStatus_t sparseSolve(
                 cusolverSpHandle_t handle,
                 int size,
                 const float * b,
                 float * x,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      cusolverStatus_t status = cusolverSpScsrcholSolve(
        handle, size, b, x, chol_info, buffer);
      cudaDeviceSynchronize();
      return status;
    }
  };

#ifdef HAVE_TEUCHOS_COMPLEX
  template <>
  struct FunctionMap<cuSOLVER,Kokkos::complex<double>>
  {
    static cusolverStatus_t bufferInfo(
                 cusolverDnHandle_t handle,
                 int n,
                 void * A,
                 int lda,
                 int * lwork)
    {
      return cusolverDnZgetrf_bufferSize(
        handle, n, n, reinterpret_cast<cuDoubleComplex*>(A), lda, lwork);
    }

    static cusolverStatus_t sparseBufferInfo(
                 cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 cusparseMatDescr_t & desc,
                 const void * values,
                 const int * rowPtr,
                 const int * colIdx,
                 csrcholInfo_t & chol_info,
                 size_t * internalDataInBytes,
                 size_t * workspaceInBytes)
    {
      typedef cuDoubleComplex scalar_t;
      const scalar_t * cu_values = reinterpret_cast<const scalar_t *>(values);
      return cusolverSpZcsrcholBufferInfo(handle, size, nnz, desc,
        cu_values, rowPtr, colIdx, chol_info,
        internalDataInBytes, workspaceInBytes);
    }

    static cusolverStatus_t numeric(
                 cusolverDnHandle_t handle,
                 int n,
                 void * A,
                 int lda,
                 void * work,
                 int * ipiv,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnZgetrf(
        handle, n, n,
        reinterpret_cast<cuDoubleComplex*>(A), lda,
        reinterpret_cast<cuDoubleComplex*>(work),
        ipiv, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cusolverStatus_t sparseNumeric(
                 cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 cusparseMatDescr_t & desc,
                 const void * values,
                 const int * rowPtr,
                 const int * colIdx,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      typedef cuDoubleComplex scalar_t;
      const scalar_t * cu_values =
        reinterpret_cast<const scalar_t *>(values);
      cusolverStatus_t status = cusolverSpZcsrcholFactor(
        handle, size, nnz, desc, cu_values, rowPtr, colIdx, chol_info, buffer);
      cudaDeviceSynchronize();
      return status;
    }

    static cusolverStatus_t solveLU(
                 cusolverDnHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const void * LU,
                 int ldlu,
                 const int * ipiv,
                 void * B,
                 int ldb,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnZgetrs(
        handle, trans, n, nrhs,
        reinterpret_cast<const cuDoubleComplex*>(LU), ldlu, ipiv,
        reinterpret_cast<cuDoubleComplex*>(B), ldb, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cublasStatus_t solve(
                 cublasHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const void * inverse,
                 int ldinv,
                 const void * B,
                 int ldb,
                 void * X,
                 int ldx)
    {
      cuDoubleComplex one  = make_cuDoubleComplex(1, 0);
      cuDoubleComplex zero = make_cuDoubleComplex(0, 0);
      return cublasZgemm(handle, trans, CUBLAS_OP_N, n, nrhs, n,
        &one,  reinterpret_cast<const cuDoubleComplex*>(inverse), ldinv,
               reinterpret_cast<const cuDoubleComplex*>(B), ldb,
        &zero, reinterpret_cast<cuDoubleComplex*>(X), ldx);
    }

    static cusolverStatus_t sparseSolve(
                 cusolverSpHandle_t handle,
                 int size,
                 const void * b,
                 void * x,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      typedef cuDoubleComplex scalar_t;
      const scalar_t * cu_b = reinterpret_cast<const scalar_t *>(b);
      scalar_t * cu_x = reinterpret_cast<scalar_t *>(x);
      cusolverStatus_t status = cusolverSpZcsrcholSolve(
        handle, size, cu_b, cu_x, chol_info, buffer);
      cudaDeviceSynchronize();
      return status;
    }
  };

  template <>
  struct FunctionMap<cuSOLVER,Kokkos::complex<float>>
  {
    static cusolverStatus_t bufferInfo(
                 cusolverDnHandle_t handle,
                 int n,
                 void * A,
                 int lda,
                 int * lwork)
    {
      return cusolverDnCgetrf_bufferSize(
        handle, n, n, reinterpret_cast<cuFloatComplex*>(A), lda, lwork);
    }

    static cusolverStatus_t sparseBufferInfo(
                 cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 cusparseMatDescr_t & desc,
                 const void * values,
                 const int * rowPtr,
                 const int * colIdx,
                 csrcholInfo_t & chol_info,
                 size_t * internalDataInBytes,
                 size_t * workspaceInBytes)
    {
      typedef cuFloatComplex scalar_t;
      const scalar_t * cu_values = reinterpret_cast<const scalar_t *>(values);
      return cusolverSpCcsrcholBufferInfo(handle, size, nnz, desc,
        cu_values, rowPtr, colIdx, chol_info,
        internalDataInBytes, workspaceInBytes);
    }

    static cusolverStatus_t numeric(
                 cusolverDnHandle_t handle,
                 int n,
                 void * A,
                 int lda,
                 void * work,
                 int * ipiv,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnCgetrf(
        handle, n, n,
        reinterpret_cast<cuFloatComplex*>(A), lda,
        reinterpret_cast<cuFloatComplex*>(work),
        ipiv, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cusolverStatus_t sparseNumeric(
                 cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 cusparseMatDescr_t & desc,
                 const void * values,
                 const int * rowPtr,
                 const int * colIdx,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      typedef cuFloatComplex scalar_t;
      const scalar_t * cu_values = reinterpret_cast<const scalar_t *>(values);
      cusolverStatus_t status = cusolverSpCcsrcholFactor(
        handle, size, nnz, desc, cu_values, rowPtr, colIdx, chol_info, buffer);
      cudaDeviceSynchronize();
      return status;
    }

    static cusolverStatus_t solveLU(
                 cusolverDnHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const void * LU,
                 int ldlu,
                 const int * ipiv,
                 void * B,
                 int ldb,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnCgetrs(
        handle, trans, n, nrhs,
        reinterpret_cast<const cuFloatComplex*>(LU), ldlu, ipiv,
        reinterpret_cast<cuFloatComplex*>(B), ldb, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cublasStatus_t solve(
                 cublasHandle_t handle,
                 cublasOperation_t trans,
                 int n,
                 int nrhs,
                 const void * inverse,
                 int ldinv,
                 const void * B,
                 int ldb,
                 void * X,
                 int ldx)
    {
      cuFloatComplex one  = make_cuFloatComplex(1, 0);
      cuFloatComplex zero = make_cuFloatComplex(0, 0);
      return cublasCgemm(handle, trans, CUBLAS_OP_N, n, nrhs, n,
        &one,  reinterpret_cast<const cuFloatComplex*>(inverse), ldinv,
               reinterpret_cast<const cuFloatComplex*>(B), ldb,
        &zero, reinterpret_cast<cuFloatComplex*>(X), ldx);
    }

    static cusolverStatus_t sparseSolve(
                 cusolverSpHandle_t handle,
                 int size,
                 const void * b,
                 void * x,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      typedef cuFloatComplex scalar_t;
      const scalar_t * cu_b = reinterpret_cast<const scalar_t *>(b);
      scalar_t * cu_x = reinterpret_cast<scalar_t *>(x);
      cusolverStatus_t status = cusolverSpCcsrcholSolve(
        handle, size, cu_b, cu_x, chol_info, buffer);
      cudaDeviceSynchronize();
      return status;
    }
  };
#endif

} // end namespace Amesos2

#endif  // AMESOS2_CUSOLVER_FUNCTIONMAP_HPP
