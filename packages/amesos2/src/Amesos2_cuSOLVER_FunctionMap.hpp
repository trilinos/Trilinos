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
      cusolverStatus_t status =
        cusolverSpDcsrcholBufferInfo(handle, size, nnz, desc, values,
          rowPtr, colIdx, chol_info, internalDataInBytes, workspaceInBytes);
      return status;
    }

    static cusolverStatus_t numeric(
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
      return status;
    }

    static cusolverStatus_t solve(
                 cusolverSpHandle_t handle,
                 int size,
                 const double * b,
                 double * x,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      cusolverStatus_t status = cusolverSpDcsrcholSolve(
        handle, size, b, x, chol_info, buffer);
      return status;
    }
  };

  template <>
  struct FunctionMap<cuSOLVER,float>
  {
    static cusolverStatus_t bufferInfo(
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
      cusolverStatus_t status =
        cusolverSpScsrcholBufferInfo(handle, size, nnz, desc, values,
          rowPtr, colIdx, chol_info, internalDataInBytes, workspaceInBytes);
      return status;
    }

    static cusolverStatus_t numeric(
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
      return status;
    }

    static cusolverStatus_t solve(
                 cusolverSpHandle_t handle,
                 int size,
                 const float * b,
                 float * x,
                 csrcholInfo_t & chol_info,
                 void * buffer)
    {
      cusolverStatus_t status = cusolverSpScsrcholSolve(
        handle, size, b, x, chol_info, buffer);
      return status;
    }
  };

#ifdef HAVE_TEUCHOS_COMPLEX
  template <>
  struct FunctionMap<cuSOLVER,Kokkos::complex<double>>
  {
    static cusolverStatus_t bufferInfo(
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
      cusolverStatus_t status =
        cusolverSpZcsrcholBufferInfo(handle, size, nnz, desc,
          cu_values, rowPtr, colIdx, chol_info,
          internalDataInBytes, workspaceInBytes);
      return status;
    }

    static cusolverStatus_t numeric(
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
      return status;
    }

    static cusolverStatus_t solve(
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
      return status;
    }
  };

  template <>
  struct FunctionMap<cuSOLVER,Kokkos::complex<float>>
  {
    static cusolverStatus_t bufferInfo(
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
      cusolverStatus_t status =
        cusolverSpCcsrcholBufferInfo(handle, size, nnz, desc,
          cu_values, rowPtr, colIdx, chol_info,
          internalDataInBytes, workspaceInBytes);
      return status;
    }

    static cusolverStatus_t numeric(
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
      return status;
    }

    static cusolverStatus_t solve(
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
      return status;
    }
  };
#endif

} // end namespace Amesos2

#endif  // AMESOS2_CUSOLVER_FUNCTIONMAP_HPP
