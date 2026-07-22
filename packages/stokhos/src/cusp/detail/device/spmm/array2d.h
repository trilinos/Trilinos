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
#include <cusp/array2d.h>
#include <cusp/array1d.h>
#include <cuda.h>
#include <iostream>

namespace cusp
{
namespace detail
{
namespace device
{


template <typename IndexType,
          typename ValueType>
__global__ void
row_axpby_kernel(const ValueType a, const ValueType b, const IndexType num_rows,
                 const IndexType num_cols,
                 const ValueType * x,
                 const ValueType * y,
                 ValueType * z)
{
  for (IndexType row = threadIdx.y + blockDim.y * blockIdx.x; row < num_rows;
       row += gridDim.x * blockDim.y)
    for (IndexType j = threadIdx.x; j < num_cols; j+=blockDim.x )
      z[j+num_cols*row]= a * x[j+num_cols*row] + b * y[j+num_cols*row];
}

template <typename IndexType,
          typename ValueType>
__global__ void
col_axpby_kernel(const ValueType a, const ValueType b, const IndexType num_rows,
                 const IndexType num_cols,
                 const ValueType * x,
                 const ValueType * y,
                 ValueType * z)
{
  for (IndexType j = threadIdx.y; j < num_cols; j+=blockDim.y )
    for (IndexType row = threadIdx.x + blockDim.x * blockIdx.x; row < num_rows;
         row += gridDim.x * blockDim.x)
      z[j*num_rows+row]= a * x[j*num_rows+row] + b * y[j*num_rows+row];
}

template <typename IndexType,
          typename ValueType>
__global__ void
row_spmm_dense_diag_kernel(const IndexType Anum_rows,
                           const IndexType Anum_cols,
                           const ValueType * Aval,
                           const ValueType * x,
                           ValueType * y)
{
  for (IndexType row = threadIdx.y + blockDim.y * blockIdx.x; row < Anum_rows;
       row += gridDim.x * blockDim.y)
    for (IndexType j = threadIdx.x; j < Anum_cols; j+=blockDim.x )
      y[j+Anum_cols*row]= x[j] * Aval[j+Anum_cols*row];
}

template <typename IndexType,
          typename ValueType>
__global__ void
col_spmm_dense_diag_kernel(const IndexType Anum_rows,
                           const IndexType Anum_cols,
                           const ValueType * Aval,
                           const ValueType * x,
                           ValueType * y)
{
  for (IndexType j = threadIdx.y; j < Anum_cols; j+=blockDim.y )
    for (IndexType row = threadIdx.x + blockDim.x * blockIdx.x; row < Anum_rows;
         row += gridDim.x * blockDim.x)
      y[j*Anum_rows+row]= x[j] * Aval[j*Anum_rows+row];
}

template <typename IndexType,
          typename ValueType>
__global__ void
row_spmm_MVdot_kernel(IndexType ROWS_PER_BLOCK, const IndexType Anum_rows,
                      const IndexType Anum_cols,
                      const ValueType * Aval,
                      const ValueType * x,
                      ValueType * temp)
{
  extern __shared__ int sh[];
  ValueType * const sdata = (ValueType*) sh;

  for (IndexType col = threadIdx.x; col < Anum_cols; col += blockDim.x){
    sdata[col + Anum_cols*threadIdx.y] = 0;
  }
  for (IndexType row = threadIdx.y + blockDim.y * blockIdx.x; row < Anum_rows;
       row += gridDim.x * blockDim.y) {
    for (IndexType col = threadIdx.x; col < Anum_cols; col += blockDim.x){
      sdata[col + Anum_cols*threadIdx.y] +=
        Aval[col + Anum_cols * row] * x[col + Anum_cols * row];
    }
  }
  //sum all local sums together to get
  __syncthreads();
  IndexType nwarp = blockDim.y / 2;
  while (nwarp > 0) {
    for (IndexType col = threadIdx.x; col < Anum_cols; col += blockDim.x){
      IndexType j = threadIdx.y;
      if (j < nwarp){
        IndexType j2 = j+nwarp;
        sdata[col+j*Anum_cols] +=  sdata[col+j2*Anum_cols];
      }
      __syncthreads();
    }
    nwarp /= 2;
  }
  __syncthreads();
  for (IndexType col = threadIdx.x; col < Anum_cols; col+=blockDim.x){
    temp[ Anum_cols * blockIdx.x + col ] = sdata[col];
  }

}

template <typename IndexType,
          typename ValueType>
__global__ void
col_spmm_MVdot_kernel(const IndexType Anum_rows,
                      const IndexType Anum_cols,
                      const ValueType * Aval,
                      const ValueType * x,
                      ValueType * temp)
{
  extern __shared__ int sh[];
  volatile ValueType * const sdata = (ValueType*) sh;

  IndexType tid = threadIdx.x + threadIdx.y*blockDim.x;

  for (IndexType col = threadIdx.y; col < Anum_cols; col += blockDim.y) {

    ValueType s = 0.0;
    for (IndexType row = threadIdx.x+blockDim.x*blockIdx.x; row < Anum_rows;
         row += gridDim.x*blockDim.x) {
      s += Aval[col*Anum_rows+row] * x[col*Anum_rows+row];
    }

    // sum all local sums together to get
    sdata[tid] = s;
    if (threadIdx.x+16 < blockDim.x) sdata[tid] += sdata[tid + 16];
    if (threadIdx.x+ 8 < blockDim.x) sdata[tid] += sdata[tid +  8];
    if (threadIdx.x+ 4 < blockDim.x) sdata[tid] += sdata[tid +  4];
    if (threadIdx.x+ 2 < blockDim.x) sdata[tid] += sdata[tid +  2];
    if (threadIdx.x+ 1 < blockDim.x) sdata[tid] += sdata[tid +  1];

    // store results
    if (threadIdx.x == 0)
      temp[col*gridDim.x+blockIdx.x] = sdata[tid];
  }

}

//kernel to compute final dot product by adding the sums across the blocks
template <typename IndexType,
          typename ValueType>
__global__ void
row_reduce_kernel(const IndexType num_rows,
                  const IndexType num_cols,
                  const ValueType * temp,
                  ValueType * y)
{
  extern __shared__ int sh[];
  ValueType * const sdata = (ValueType*) sh;
  for (IndexType col = threadIdx.x; col < num_cols; col += blockDim.x){
    sdata[col + num_cols*threadIdx.y] = 0;
  }
  for (IndexType row = threadIdx.y + blockDim.y * blockIdx.x; row < num_rows;
       row += gridDim.x * blockDim.y) {
    for (IndexType col = threadIdx.x; col < num_cols; col += blockDim.x){
      sdata[col + num_cols*threadIdx.y] +=  temp[col + num_cols * row];
    }
  }
  __syncthreads();
  IndexType nwarp = blockDim.y / 2;
  while (nwarp > 0) {
    for (IndexType col = threadIdx.x; col < num_cols; col += blockDim.x){
      IndexType j = threadIdx.y;
      if (j < nwarp){
        IndexType j2 = j+nwarp;
        sdata[col+j*num_cols] +=  sdata[col+j2*num_cols];
      }
      __syncthreads();
    }
    nwarp /= 2;
  }
  __syncthreads();
  for (IndexType col = threadIdx.x; col < num_cols; col+=blockDim.x){
    y[col + num_cols * blockIdx.x] = sdata[col];
  }
}

//kernel to compute final dot product by adding the sums across the blocks
template <typename IndexType,
          typename ValueType>
__global__ void
col_reduce_kernel(const IndexType num_rows,
                  const IndexType num_cols,
                  const ValueType * temp,
                  ValueType * y)
{
  extern __shared__ int sh[];
  volatile ValueType * const sdata = (ValueType*) sh;

  IndexType tid = threadIdx.x + threadIdx.y*blockDim.x;

  for (IndexType col = threadIdx.y; col < num_cols; col += blockDim.y) {

    ValueType s = 0.0;
    for (IndexType row = threadIdx.x; row < num_rows; row += blockDim.x) {
      s += temp[col*num_rows+row];
    }

    // sum all local sums together to get
    sdata[tid] = s;
    if (threadIdx.x+16 < blockDim.x) sdata[tid] += sdata[tid + 16];
    if (threadIdx.x+ 8 < blockDim.x) sdata[tid] += sdata[tid +  8];
    if (threadIdx.x+ 4 < blockDim.x) sdata[tid] += sdata[tid +  4];
    if (threadIdx.x+ 2 < blockDim.x) sdata[tid] += sdata[tid +  2];
    if (threadIdx.x+ 1 < blockDim.x) sdata[tid] += sdata[tid +  1];

    // store results
    if (threadIdx.x == 0)
      y[col] = sdata[tid];
  }
}


template <typename Vector1,
          typename Vector2,
          typename Vector3>
void __spmm_MVdot(const int numRHS, const Vector1& A,
                  const Vector2& x,
                  Vector3& y, cusp::row_major)
{
  CUSP_PROFILE_SCOPED();
  typedef typename Vector2::index_type   IndexType;
  typedef typename Vector2::value_type   ValueType;
  typedef typename Vector2::memory_space MemorySpace;
  const size_t BLOCK_SIZE = 32;

  const size_t ROWS_PER_BLOCK = 8;

  dim3 block(BLOCK_SIZE, ROWS_PER_BLOCK);
  const size_t shared =    block.y * numRHS * sizeof(ValueType);


  const size_t MAX_BLOCKS = 96;
  const size_t NUM_BLOCKS = std::min(MAX_BLOCKS, DIVIDE_INTO(A.num_rows, ROWS_PER_BLOCK));
  dim3 grid(NUM_BLOCKS, 1);


  cusp::array2d<ValueType, MemorySpace, cusp::row_major> temp(NUM_BLOCKS, A.num_cols, 0);

  row_spmm_MVdot_kernel<IndexType,ValueType> <<<grid, block, shared>>>
    (ROWS_PER_BLOCK, A.num_rows, A.num_cols,
     thrust::raw_pointer_cast(&A.values[0]),
     thrust::raw_pointer_cast(&x.values[0]),
     thrust::raw_pointer_cast(&temp.values[0]));

  // add rows of temp (sum from each block) for final answer
  IndexType num_rows = NUM_BLOCKS;
  IndexType num_blocks = 16;
  cusp::array2d<ValueType, MemorySpace, cusp::row_major> sum_temp(num_blocks, numRHS, 0);
  //Need size to be power of 2
  while (num_blocks > 0){
    dim3 sum_grid(num_blocks, 1);
    IndexType size = 1;
    while (size < num_rows / num_blocks){
      size <<=1;
    }
    dim3 sum_block( 32, size);
    size_t sum_shared = sum_block.y * numRHS * sizeof(ValueType);
    sum_temp.resize(num_blocks, numRHS);
    row_reduce_kernel<IndexType,ValueType> <<<sum_grid, sum_block, sum_shared>>>
      (temp.num_rows, temp.num_cols, thrust::raw_pointer_cast(&temp.values[0]),
       thrust::raw_pointer_cast(&sum_temp.values[0]));
    cusp::copy(sum_temp, temp);
    num_rows = num_blocks;
    num_blocks /= 2;
  }
  for (IndexType i = 0; i < numRHS; i++){
    y[i] = sum_temp(0, i);
  }
}

template <typename Vector1,
          typename Vector2,
          typename Vector3>
void __spmm_MVdot(const int numRHS, const Vector1& A,
                  const Vector2& x,
                  Vector3& y, cusp::column_major)
{
#if 0
  typedef typename Vector2::index_type   IndexType;
  for (IndexType col=0; col<A.num_cols; ++col) {
    y[col] = cusp::blas::dotc(A.column(col), x.column(col));
  }
#else
  CUSP_PROFILE_SCOPED();
  typedef typename Vector2::index_type   IndexType;
  typedef typename Vector2::value_type   ValueType;
  typedef typename Vector2::memory_space MemorySpace;

  const size_t BLOCK_SIZE = 32;
  const size_t COLS_PER_BLOCK = 8;

  dim3 block(BLOCK_SIZE, COLS_PER_BLOCK);
  const size_t shared = block.x * block.y * sizeof(ValueType);

  const size_t MAX_BLOCKS = 96;
  const size_t NUM_BLOCKS = std::min(MAX_BLOCKS, DIVIDE_INTO(A.num_rows, BLOCK_SIZE));
  dim3 grid(NUM_BLOCKS, 1);


  cusp::array2d<ValueType, MemorySpace, cusp::column_major> temp(NUM_BLOCKS, A.num_cols, 0);

  col_spmm_MVdot_kernel<IndexType,ValueType> <<<grid, block, shared>>>
    (A.num_rows, A.num_cols,
     thrust::raw_pointer_cast(&A.values[0]),
     thrust::raw_pointer_cast(&x.values[0]),
     thrust::raw_pointer_cast(&temp.values[0]));

  // add rows of temp (sum from each block) for final answer
  dim3 sum_grid(1, 1);
  col_reduce_kernel<IndexType,ValueType> <<<sum_grid, block, shared>>>
    (NUM_BLOCKS, A.num_cols,
     thrust::raw_pointer_cast(&temp.values[0]),
     thrust::raw_pointer_cast(&y[0]));
#endif
}

template <typename Vector1,
          typename Vector2,
          typename Vector3>
void __spmm_dense_diag(const Vector1& A,
                       const Vector2& x,
                       Vector3& y,
                       cusp::row_major)
{
  CUSP_PROFILE_SCOPED();
  typedef typename Vector3::index_type   IndexType;
  typedef typename Vector3::value_type   ValueType;
  typedef typename Vector3::memory_space MemorySpace;
  const size_t BLOCK_SIZE = 32;
  const size_t ROWS_PER_BLOCK = 8;
  const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(row_spmm_dense_diag_kernel<IndexType, ValueType>, BLOCK_SIZE, (size_t) 0);
  const size_t NUM_BLOCKS = std::min(MAX_BLOCKS, DIVIDE_INTO(A.num_rows, ROWS_PER_BLOCK));
  dim3 block(BLOCK_SIZE, ROWS_PER_BLOCK);
  dim3 grid(NUM_BLOCKS, 1);

  row_spmm_dense_diag_kernel<IndexType,ValueType> <<<grid, block>>>
    (A.num_rows, A.num_cols,
     thrust::raw_pointer_cast(&A.values[0]),
     thrust::raw_pointer_cast(&x[0]),
     thrust::raw_pointer_cast(&(y.values)[0]));
}

template <typename Vector1,
          typename Vector2,
          typename Vector3>
void __spmm_dense_diag(const Vector1& A,
                       const Vector2& x,
                       Vector3& y,
                       cusp::column_major)
{
#if 0
  typedef typename Vector1::index_type IndexType;
  typedef typename Vector1::value_type ValueType;
  for (IndexType col=0; col<A.num_cols; ++col) {
    cusp::blas::axpby(A.column(col), A.column(col), y.column(col),
                      x[col], ValueType(0.0));
  }
#else
  CUSP_PROFILE_SCOPED();
  typedef typename Vector3::index_type   IndexType;
  typedef typename Vector3::value_type   ValueType;
  typedef typename Vector3::memory_space MemorySpace;
  const size_t BLOCK_SIZE = 32;
  const size_t COLS_PER_BLOCK = 8;
  const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(col_spmm_dense_diag_kernel<IndexType, ValueType>, BLOCK_SIZE, (size_t) 0);
  const size_t NUM_BLOCKS = std::min(MAX_BLOCKS, DIVIDE_INTO(A.num_rows, BLOCK_SIZE));
  dim3 block(BLOCK_SIZE, COLS_PER_BLOCK, 1);
  dim3 grid(NUM_BLOCKS, 1);

  col_spmm_dense_diag_kernel<IndexType,ValueType> <<<grid, block>>>
    (A.num_rows, A.num_cols,
     thrust::raw_pointer_cast(&A.values[0]),
     thrust::raw_pointer_cast(&x[0]),
     thrust::raw_pointer_cast(&(y.values)[0]));
#endif
}

template <typename ValueType,
          typename Vector1,
          typename Vector2>
void __spmm_axpby(const ValueType& a,
                  const Vector1& x,
                  const ValueType& b,
                  const Vector1& y,
                  Vector2& z,
                  cusp::row_major)
{
  CUSP_PROFILE_SCOPED();
  typedef typename Vector2::index_type   IndexType;
  typedef typename Vector2::memory_space MemorySpace;
  const size_t BLOCK_SIZE = 32;
  const size_t ROWS_PER_BLOCK = 6;
  const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(row_axpby_kernel<IndexType, ValueType>, BLOCK_SIZE, (size_t) 0);
  const size_t NUM_BLOCKS = std::min(MAX_BLOCKS, DIVIDE_INTO(x.num_rows, ROWS_PER_BLOCK));
  dim3 block(BLOCK_SIZE, ROWS_PER_BLOCK);
  dim3 grid(NUM_BLOCKS, 1);

  row_axpby_kernel<IndexType,ValueType> <<<grid, block>>>
    (a,b, x.num_rows, x.num_cols,
     thrust::raw_pointer_cast(&x.values[0]),
     thrust::raw_pointer_cast(&(y.values)[0]),
     thrust::raw_pointer_cast(&(z.values)[0]));
}

template <typename ValueType,
          typename Vector1,
          typename Vector2>
void __spmm_axpby(const ValueType& a,
                  const Vector1& x,
                  const ValueType& b,
                  const Vector1& y,
                  Vector2& z,
                  cusp::column_major)
{
#if 0
  typedef typename Vector1::index_type IndexType;
  for (IndexType col=0; col<x.num_cols; ++col) {
    cusp::blas::axpby(x.column(col), y.column(col), z.column(col), a, b);
  }
#else
  CUSP_PROFILE_SCOPED();
  typedef typename Vector2::index_type   IndexType;
  typedef typename Vector2::memory_space MemorySpace;
  const size_t BLOCK_SIZE = 32;
  const size_t COLS_PER_BLOCK = 6;
  const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(col_axpby_kernel<IndexType, ValueType>, BLOCK_SIZE, (size_t) 0);
  const size_t NUM_BLOCKS = std::min(MAX_BLOCKS, DIVIDE_INTO(x.num_rows, BLOCK_SIZE));
  dim3 block(BLOCK_SIZE, COLS_PER_BLOCK, 1);
  dim3 grid(NUM_BLOCKS, 1);

  col_axpby_kernel<IndexType,ValueType> <<<grid, block>>>
    (a,b, x.num_rows, x.num_cols,
     thrust::raw_pointer_cast(&x.values[0]),
     thrust::raw_pointer_cast(&(y.values)[0]),
     thrust::raw_pointer_cast(&(z.values)[0]));
#endif
}

template <typename Vector1,
          typename Vector2,
          typename Vector3>
void spmm_MVdot(const Vector1& A,
                const Vector2& x,
                Vector3& y)
{
  //Determine if row-wise or col-wise then call appropriate multiply
  __spmm_MVdot(A.num_cols, A, x, y, typename Vector2::orientation());
}

template <typename Vector1,
          typename Vector2,
          typename Vector3>
void spmm_dense_diag(const Vector1& A,
                     const Vector2& x,
                     Vector3& y)
{
  //Determine if row-wise or col-wise then call appropriate multiply
  __spmm_dense_diag(A, x, y, typename Vector1::orientation());
}

template <typename ValueType,
          typename Vector1,
          typename Vector2>
void spmm_axpby(const ValueType& a, const Vector1& x, const ValueType& b,
                const Vector1& y,
                Vector2& z)
{
  __spmm_axpby(a, x, b, y, z, typename Vector1::orientation());
}


} // end namespace device
} // end namespace detail
} // end namespace cusp
