// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 *  Copyright 2008-2009 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#pragma once

#include <cusp/format.h>
#include <cusp/csr_matrix.h>
#include <cusp/detail/functional.h>
#include <iostream>

// SpMM
#include <cusp/detail/device/spmm/csr_vector.h>
#include <cusp/detail/device/spmm/array2d.h>

#include "Teuchos_TimeMonitor.hpp"

namespace cusp
{
namespace detail
{
namespace device
{

////////////////////////////////////////
// Sparse Matrix-BlockVector Multiply //
////////////////////////////////////////
template <typename Matrix,
         typename Vector1,
         typename Vector2>
void MVmultiply(const Matrix&  A,
              const Vector1& B,
              Vector2& C,
              cusp::sparse_format,
              cusp::array2d_format,
              cusp::array2d_format
              )
{
  TEUCHOS_FUNC_TIME_MONITOR("CUSP Matrix block-apply");
  cusp::detail::device::spmm_csr_vector(A,B,C);
  cudaDeviceSynchronize();
}

////////////////////////////////////////
// Sparse Matrix-BlockVector Multiply //
////////////////////////////////////////
template <typename Matrix,
         typename Vector1,
         typename Vector2>
void OVmultiply(const Matrix&  A,
              const Vector1& B,
              Vector2& C,
              cusp::sparse_format,
              cusp::array2d_format,
              cusp::array2d_format
              )
{
  TEUCHOS_FUNC_TIME_MONITOR("CUSP Operator block-apply");
  cusp::detail::device::spmm_csr_vector(A,B,C);
  cudaDeviceSynchronize();
}

/////////////////////////////////////////
//// Dense Matrix-Matrix Multiplication where B is diagonal and entries stored in 1d array//
///////////////////////////////////////////
template <typename Matrix1,
         typename Matrix2,
         typename Matrix3>
void MVmultiply(const Matrix1& A,
              const Matrix2& B,
              Matrix3& C,
              cusp::array2d_format,
              cusp::array1d_format,
              cusp::array2d_format)
{
  TEUCHOS_FUNC_TIME_MONITOR("CUSP Dense-diag");
  cusp::detail::device::spmm_dense_diag(A,B,C);
  cudaDeviceSynchronize();
}

/////////////////////////////////////////
// Dot Product: Computes C[i] = A[i] ^ T * B[i] where A[i] and B[i] are the i-th columns of A and B
/////////////////////////////////////////////


template <typename MV,
         typename MV1,
         typename MV2>
void MVdot(const MV&  A,
              const MV1& B,
              MV2& C)
{
  TEUCHOS_FUNC_TIME_MONITOR("CUSP dot");
  cusp::detail::device::spmm_MVdot(A,B,C);
  cudaDeviceSynchronize();
}

/////////////////////////////////////////
//// Compute a*X + b * Y where X, Y are array2d
/////////////////////////////////////////////
template <typename ValueType,
         typename MV1,
         typename MV2>
void axpby(const ValueType&  A, const MV1& X, const ValueType&  B,
              const MV1& Y,
              MV2& Z)
{
  TEUCHOS_FUNC_TIME_MONITOR("CUSP axpby");
  cusp::detail::device::spmm_axpby(A,X,B,Y,Z);
  cudaDeviceSynchronize();
}


/////////////////
// Entry Point //
/////////////////
template <typename Matrix,
         typename MatrixOrVector1,
         typename MatrixOrVector2>
void MVmultiply(const Matrix&  A,
              const MatrixOrVector1& B,
              MatrixOrVector2& C)
{
  cusp::detail::device::MVmultiply(A, B, C,
                                   typename Matrix::format(),
                                   typename MatrixOrVector1::format(),
                                   typename MatrixOrVector2::format());
}

template <typename Matrix,
         typename MatrixOrVector1,
         typename MatrixOrVector2>
void OVmultiply(const Matrix&  A,
              const MatrixOrVector1& B,
              MatrixOrVector2& C)
{
  cusp::detail::device::OVmultiply(A, B, C,
                                   typename Matrix::format(),
                                   typename MatrixOrVector1::format(),
                                   typename MatrixOrVector2::format());
}

} // end namespace device
} // end namespace detail
} // end namespace cusp
