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

#include <cusp/array1d.h>
#include <cusp/array2d.h>
#include <cusp/linear_operator.h>

#include <cmath>

namespace cusp
{
namespace detail
{

template <typename IndexType, typename ValueType, typename MemorySpace, typename OrientationA>
int block_lu_factor(cusp::array2d<ValueType,MemorySpace,OrientationA>& A,
                    cusp::array1d<IndexType,MemorySpace>& pivot)
{
  const int n = A.num_rows;

  // For each row and column, k = 0, ..., n-1,
  for (int k = 0; k < n; k++)
  {
    // find the pivot row
    pivot[k] = k;
    ValueType max = std::fabs(A(k,k));

    for (int j = k + 1; j < n; j++)
    {
      if (max < std::fabs(A(j,k)))
      {
        max = std::fabs(A(j,k));
        pivot[k] = j;
      }
    }

    // and if the pivot row differs from the current row, then
    // interchange the two rows.
    if (pivot[k] != k)
      for (int j = 0; j < n; j++)
        std::swap(A(k,j), A(pivot[k],j));

    // and if the matrix is singular, return error
    if (A(k,k) == 0.0)
      return -1;

    // otherwise find the lower triangular matrix elements for column k.
    for (int i = k + 1; i < n; i++)
      A(i,k) /= A(k,k);

    // update remaining matrix
    for (int i = k + 1; i < n; i++)
      for (int j = k + 1; j < n; j++)
        A(i,j) -= A(i,k) * A(k,j);
  }
  return 0;
}


//LU solve for multiple right hand sides
template <typename IndexType, typename ValueType, typename MemorySpace,
          typename OrientationA, typename OrientationB>
int block_lu_solve(const cusp::array2d<ValueType,MemorySpace,OrientationA>& A,
                   const cusp::array1d<IndexType,MemorySpace>& pivot,
                   const cusp::array2d<ValueType,MemorySpace,OrientationB>& b,
                   cusp::array2d<ValueType,MemorySpace,OrientationB>& x,
                   cusp::array2d_format)
{
  const int n = A.num_rows;
  const int numRHS = b.num_cols;
  // copy rhs to x
  x = b;
  // Solve the linear equation Lx = b for x, where L is a lower triangular matrix

  for (int k = 0; k < n; k++)
  {
    if (pivot[k] != k){//swap row k of x with row pivot[k]
      for (int j = 0; j < numRHS; j++)
        std::swap(x(k,j),x(pivot[k],j));
    }

    for (int i = 0; i < k; i++){
      for (int j = 0; j< numRHS; j++)
        x(k,j) -= A(k,i) * x(i,j);
    }
  }

  // Solve the linear equation Ux = y, where y is the solution
  // obtained above of Lx = b and U is an upper triangular matrix.
  for (int k = n - 1; k >= 0; k--)
  {
    for (int j = 0; j< numRHS; j++){
      for (int i = k + 1; i < n; i++){
        x(k, j) -= A(k,i) * x(i, j);
      }
      if (A(k,k) == 0)
        return -1;
      x(k,j) /= A(k,k);
    }

  }
  return 0;
}




template <typename ValueType, typename MemorySpace>
class block_lu_solver : public cusp::linear_operator<ValueType,MemorySpace>
{
  cusp::array2d<ValueType,cusp::host_memory> lu;
  cusp::array1d<int,cusp::host_memory>       pivot;

public:
  block_lu_solver() : linear_operator<ValueType,MemorySpace>() {}

  template <typename MatrixType>
  block_lu_solver(const MatrixType& A) :
  linear_operator<ValueType,MemorySpace>(A.num_rows, A.num_cols, A.num_entries)
  {
    CUSP_PROFILE_SCOPED();

    lu = A;
    pivot.resize(A.num_rows);
    block_lu_factor(lu,pivot);
  }

  template <typename VectorType1, typename VectorType2>
  void operator()(const VectorType1& b, VectorType2& x) const
  {
    CUSP_PROFILE_SCOPED();
    block_lu_solve(lu, pivot, b, x, typename VectorType2::format());
  }
};

} // end namespace detail
} // end namespace cusp
