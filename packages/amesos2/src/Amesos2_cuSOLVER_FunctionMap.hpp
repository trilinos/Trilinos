// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef AMESOS2_CUSOLVER_FUNCTIONMAP_HPP
#define AMESOS2_CUSOLVER_FUNCTIONMAP_HPP

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_cuSOLVER_TypeMap.hpp"

#ifdef HAVE_TEUCHOS_COMPLEX
  // For cudaDoubleComplex and cudaFloatComplex
  // Not namespaced because of extern C. TODO
#include <cuComplex.h>
#endif

namespace Amesos2 {

  namespace CUSOLVER {
#    include <cuda.h>
#    include <cusolverSp.h>
#    include <cusolverDn.h>
#    include <cusparse.h>
#    include <cusolverSp_LOWLEVEL_PREVIEW.h>
  }

  template <>
  struct FunctionMap<cuSOLVER,double>
  {
    static CUSOLVER::cusolverStatus_t bufferInfo(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 CUSOLVER::cusparseMatDescr_t & desc,
                 const double * values,
                 const int * rowPtr,
                 const int * colIdx,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 size_t * internalDataInBytes,
                 size_t * workspaceInBytes)
    {
      CUSOLVER::cusolverStatus_t status =
        CUSOLVER::cusolverSpDcsrcholBufferInfo(handle, size, nnz, desc, values,
          rowPtr, colIdx, chol_info, internalDataInBytes, workspaceInBytes);
      return status;
    }

    static CUSOLVER::cusolverStatus_t numeric(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 CUSOLVER::cusparseMatDescr_t & desc,
                 const double * values,
                 const int * rowPtr,
                 const int * colIdx,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 void * buffer)
    {
      CUSOLVER::cusolverStatus_t status = cusolverSpDcsrcholFactor(
        handle, size, nnz, desc, values, rowPtr, colIdx, chol_info, buffer);
      return status;
    }

    static CUSOLVER::cusolverStatus_t solve(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 const double * b,
                 double * x,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 void * buffer)
    {
      CUSOLVER::cusolverStatus_t status = CUSOLVER::cusolverSpDcsrcholSolve(
        handle, size, b, x, chol_info, buffer);
      return status;
    }
  };

  template <>
  struct FunctionMap<cuSOLVER,float>
  {
    static CUSOLVER::cusolverStatus_t bufferInfo(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 CUSOLVER::cusparseMatDescr_t & desc,
                 const float * values,
                 const int * rowPtr,
                 const int * colIdx,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 size_t * internalDataInBytes,
                 size_t * workspaceInBytes)
    {
      CUSOLVER::cusolverStatus_t status =
        CUSOLVER::cusolverSpScsrcholBufferInfo(handle, size, nnz, desc, values,
          rowPtr, colIdx, chol_info, internalDataInBytes, workspaceInBytes);
      return status;
    }

    static CUSOLVER::cusolverStatus_t numeric(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 CUSOLVER::cusparseMatDescr_t & desc,
                 const float * values,
                 const int * rowPtr,
                 const int * colIdx,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 void * buffer)
    {
      CUSOLVER::cusolverStatus_t status = cusolverSpScsrcholFactor(
        handle, size, nnz, desc, values, rowPtr, colIdx, chol_info, buffer);
      return status;
    }

    static CUSOLVER::cusolverStatus_t solve(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 const float * b,
                 float * x,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 void * buffer)
    {
      CUSOLVER::cusolverStatus_t status = CUSOLVER::cusolverSpScsrcholSolve(
        handle, size, b, x, chol_info, buffer);
      return status;
    }
  };

#ifdef HAVE_TEUCHOS_COMPLEX
  template <>
  struct FunctionMap<cuSOLVER,Kokkos::complex<double>>
  {
    static CUSOLVER::cusolverStatus_t bufferInfo(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 CUSOLVER::cusparseMatDescr_t & desc,
                 const void * values,
                 const int * rowPtr,
                 const int * colIdx,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 size_t * internalDataInBytes,
                 size_t * workspaceInBytes)
    {
      typedef cuDoubleComplex scalar_t;
      const scalar_t * cu_values = reinterpret_cast<const scalar_t *>(values);
      CUSOLVER::cusolverStatus_t status =
        CUSOLVER::cusolverSpZcsrcholBufferInfo(handle, size, nnz, desc,
          cu_values, rowPtr, colIdx, chol_info,
          internalDataInBytes, workspaceInBytes);
      return status;
    }

    static CUSOLVER::cusolverStatus_t numeric(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 CUSOLVER::cusparseMatDescr_t & desc,
                 const void * values,
                 const int * rowPtr,
                 const int * colIdx,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 void * buffer)
    {
      typedef cuDoubleComplex scalar_t;
      const scalar_t * cu_values =
        reinterpret_cast<const scalar_t *>(values);
      CUSOLVER::cusolverStatus_t status = cusolverSpZcsrcholFactor(
        handle, size, nnz, desc, cu_values, rowPtr, colIdx, chol_info, buffer);
      return status;
    }

    static CUSOLVER::cusolverStatus_t solve(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 const void * b,
                 void * x,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 void * buffer)
    {
      typedef cuDoubleComplex scalar_t;
      const scalar_t * cu_b = reinterpret_cast<const scalar_t *>(b);
      scalar_t * cu_x = reinterpret_cast<scalar_t *>(x);
      CUSOLVER::cusolverStatus_t status = CUSOLVER::cusolverSpZcsrcholSolve(
        handle, size, cu_b, cu_x, chol_info, buffer);
      return status;
    }
  };

  template <>
  struct FunctionMap<cuSOLVER,Kokkos::complex<float>>
  {
    static CUSOLVER::cusolverStatus_t bufferInfo(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 CUSOLVER::cusparseMatDescr_t & desc,
                 const void * values,
                 const int * rowPtr,
                 const int * colIdx,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 size_t * internalDataInBytes,
                 size_t * workspaceInBytes)
    {
      typedef cuFloatComplex scalar_t;
      const scalar_t * cu_values = reinterpret_cast<const scalar_t *>(values);
      CUSOLVER::cusolverStatus_t status =
        CUSOLVER::cusolverSpCcsrcholBufferInfo(handle, size, nnz, desc,
          cu_values, rowPtr, colIdx, chol_info,
          internalDataInBytes, workspaceInBytes);
      return status;
    }

    static CUSOLVER::cusolverStatus_t numeric(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 int nnz,
                 CUSOLVER::cusparseMatDescr_t & desc,
                 const void * values,
                 const int * rowPtr,
                 const int * colIdx,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 void * buffer)
    {
      typedef cuFloatComplex scalar_t;
      const scalar_t * cu_values = reinterpret_cast<const scalar_t *>(values);
      CUSOLVER::cusolverStatus_t status = cusolverSpCcsrcholFactor(
        handle, size, nnz, desc, cu_values, rowPtr, colIdx, chol_info, buffer);
      return status;
    }

    static CUSOLVER::cusolverStatus_t solve(
                 CUSOLVER::cusolverSpHandle_t handle,
                 int size,
                 const void * b,
                 void * x,
                 CUSOLVER::csrcholInfo_t & chol_info,
                 void * buffer)
    {
      typedef cuFloatComplex scalar_t;
      const scalar_t * cu_b = reinterpret_cast<const scalar_t *>(b);
      scalar_t * cu_x = reinterpret_cast<scalar_t *>(x);
      CUSOLVER::cusolverStatus_t status = CUSOLVER::cusolverSpCcsrcholSolve(
        handle, size, cu_b, cu_x, chol_info, buffer);
      return status;
    }
  };
#endif

} // end namespace Amesos2

#endif  // AMESOS2_CUSOLVER_FUNCTIONMAP_HPP
