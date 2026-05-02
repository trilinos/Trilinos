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
#include <cusolverDn.h>
#include <cublas_v2.h>

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

    // Compute explicit inverse: solves LU * inverse = I, giving inverse = A^{-1}.
    // inverse must be initialized to the identity matrix before this call.
    static cusolverStatus_t invert(
                 cusolverDnHandle_t handle,
                 int n,
                 const double * LU,
                 int ldlu,
                 const int * ipiv,
                 double * inverse,
                 int ldinv,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnDgetrs(
        handle, CUBLAS_OP_N, n, n, LU, ldlu, ipiv, inverse, ldinv, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    // Solve X = A^{-1} * B via dense matrix-matrix product.
    static cublasStatus_t solve(
                 cublasHandle_t handle,
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
      return cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        n, nrhs, n, &one, inverse, ldinv, B, ldb, &zero, X, ldx);
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

    static cusolverStatus_t invert(
                 cusolverDnHandle_t handle,
                 int n,
                 const float * LU,
                 int ldlu,
                 const int * ipiv,
                 float * inverse,
                 int ldinv,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnSgetrs(
        handle, CUBLAS_OP_N, n, n, LU, ldlu, ipiv, inverse, ldinv, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cublasStatus_t solve(
                 cublasHandle_t handle,
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
      return cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
        n, nrhs, n, &one, inverse, ldinv, B, ldb, &zero, X, ldx);
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

    static cusolverStatus_t invert(
                 cusolverDnHandle_t handle,
                 int n,
                 const void * LU,
                 int ldlu,
                 const int * ipiv,
                 void * inverse,
                 int ldinv,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnZgetrs(
        handle, CUBLAS_OP_N, n, n,
        reinterpret_cast<const cuDoubleComplex*>(LU), ldlu, ipiv,
        reinterpret_cast<cuDoubleComplex*>(inverse), ldinv, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cublasStatus_t solve(
                 cublasHandle_t handle,
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
      return cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, nrhs, n,
        &one,  reinterpret_cast<const cuDoubleComplex*>(inverse), ldinv,
               reinterpret_cast<const cuDoubleComplex*>(B), ldb,
        &zero, reinterpret_cast<cuDoubleComplex*>(X), ldx);
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

    static cusolverStatus_t invert(
                 cusolverDnHandle_t handle,
                 int n,
                 const void * LU,
                 int ldlu,
                 const int * ipiv,
                 void * inverse,
                 int ldinv,
                 int * devInfo)
    {
      cusolverStatus_t status = cusolverDnCgetrs(
        handle, CUBLAS_OP_N, n, n,
        reinterpret_cast<const cuFloatComplex*>(LU), ldlu, ipiv,
        reinterpret_cast<cuFloatComplex*>(inverse), ldinv, devInfo);
      cudaDeviceSynchronize();
      return status;
    }

    static cublasStatus_t solve(
                 cublasHandle_t handle,
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
      return cublasCgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, nrhs, n,
        &one,  reinterpret_cast<const cuFloatComplex*>(inverse), ldinv,
               reinterpret_cast<const cuFloatComplex*>(B), ldb,
        &zero, reinterpret_cast<cuFloatComplex*>(X), ldx);
    }
  };
#endif

} // end namespace Amesos2

#endif  // AMESOS2_CUSOLVER_FUNCTIONMAP_HPP
