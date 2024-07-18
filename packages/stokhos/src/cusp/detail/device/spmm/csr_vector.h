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

#include <cusp/detail/device/arch.h>
#include <cusp/detail/device/common.h>
#include <cusp/detail/device/utils.h>
#include <cusp/detail/device/texture.h>
#include <thrust/device_ptr.h>
#include <cudaProfiler.h>
#include <cuda_profiler_api.h>
#include <stdio.h>
#include <cusp/detail/device/arch.h>

#include "Stokhos_config.h"
#if 0 && defined(HAVE_STOKHOS_CUSPARSE)
#define USE_CUSPARSE_ROW 0
#define USE_CUSPARSE_COL 1
#else
#define USE_CUSPARSE_ROW 0
#define USE_CUSPARSE_COL 0
#endif

#if USE_CUSPARSE_ROW || USE_CUSPARSE_COL
#include <sstream>
#include <stdexcept>
#include <cuda_runtime.h>
#include <cusparse.h>
#endif

namespace cusp
{
namespace detail
{
namespace device
{

#if USE_CUSPARSE_ROW || USE_CUSPARSE_COL

class CudaSparseSingleton {
public:

  cusparseStatus_t   status;
  cusparseHandle_t   handle;
  cusparseMatDescr_t descra;

  static CudaSparseSingleton & singleton();

private:

  CudaSparseSingleton()
  {
    status = cusparseCreate(&handle);
    if(status != CUSPARSE_STATUS_SUCCESS)
    {
      throw std::runtime_error( std::string("ERROR - CUSPARSE Library Initialization failed" ) );
    }

    status = cusparseCreateMatDescr(&descra);
    if(status != CUSPARSE_STATUS_SUCCESS)
    {
      throw std::runtime_error( std::string("ERROR - CUSPARSE Library Matrix descriptor failed" ) );
    }

    cusparseSetMatType(descra , CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descra , CUSPARSE_INDEX_BASE_ZERO);
  }

  CudaSparseSingleton( const CudaSparseSingleton & );
  CudaSparseSingleton & operator = ( const CudaSparseSingleton & );
};

CudaSparseSingleton & CudaSparseSingleton::singleton()
{
  static CudaSparseSingleton s ; return s ;
}

#endif

#if USE_CUSPARSE_ROW

void __spmm_csr_vector(
  const csr_matrix<int,double,device_memory>& A,
  const array2d<double, device_memory, column_major>& x,
  array2d<double, device_memory, column_major>& y,
  cusp::row_major)
{
  CudaSparseSingleton & s = CudaSparseSingleton::singleton();
  const double alpha = 1 , beta = 0 ;
  cusparseStatus_t status =
    cusparseDcsrmm2(s.handle,
                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                    CUSPARSE_OPERATION_TRANSPOSE,
                    A.num_rows, x.num_cols, A.num_cols, A.num_entries,
                    alpha,
                    s.descra,
                    thrust::raw_pointer_cast(&A.values[0]),
                    thrust::raw_pointer_cast(&A.row_offsets[0]),
                    thrust::raw_pointer_cast(&A.column_indices[0]),
                    thrust::raw_pointer_cast(&(x.values)[0]),
                    x.num_rows,
                    beta,
                    thrust::raw_pointer_cast(&(y.values)[0]),
                    y.num_rows);

  if ( CUSPARSE_STATUS_SUCCESS != status ) {
    throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
  }
}

#else

template <typename IndexType, typename ValueType, unsigned MAX_NNZ_PER_ROW>
__global__ void
spmm_csr_vector_kernel_row(const IndexType Anum_rows,
                           const IndexType xnum_rows,
                           const IndexType xnum_cols,
                           const IndexType * Ar,
                           const IndexType * Ac,
                           const ValueType * Aval,
                           const ValueType * x,
                           ValueType * y)
{

  // Allocate storage for shared A row
  extern __shared__ int sh[];
  volatile ValueType * const sh_Aval = (ValueType*) sh;
  volatile IndexType * const sh_Ac   = (IndexType*) (sh_Aval + MAX_NNZ_PER_ROW * blockDim.y);

  // Global row
  const IndexType row = threadIdx.y + blockDim.y * blockIdx.x;
  if (row < Anum_rows) {
    const IndexType row_start = Ar[row];
    const IndexType row_end   = Ar[row+1];

    // Initialize y for this row
    for (IndexType j=threadIdx.x; j<xnum_cols; j+=blockDim.x) {
      y[j+xnum_cols*row] = 0.0;
    }

    // Loop over cols of A MAX_NNZ_PER_ROW at a time
    for (IndexType block_col=row_start; block_col<row_end;
         block_col+=MAX_NNZ_PER_ROW) {
      const IndexType r = block_col + MAX_NNZ_PER_ROW < row_end ?
       MAX_NNZ_PER_ROW : row_end-block_col;

      // Read row of A into shared mem using all threads in block
      // Store in shared-memory in a transposed layout to avoid
      // bank conflicts between warps (the kernel is completely bandwidth
      // limited, so this doesn't really help performance).
      // And we don't need synchronization since blockDim.x <= warp_size
      for (IndexType i=threadIdx.x; i<r; i+=blockDim.x){
        sh_Aval[i*blockDim.y+threadIdx.y] = Aval[i+block_col];
        sh_Ac[  i*blockDim.y+threadIdx.y] = Ac  [i+block_col];
      }

      // Loop over cols of x
      for (IndexType j=threadIdx.x; j<xnum_cols; j+=blockDim.x){

        // Loop over cols of A
        ValueType sum = 0.0;
        for (IndexType jj=0; jj<r; jj++){
          IndexType J = sh_Ac[jj*blockDim.y+threadIdx.y];
          sum += sh_Aval[jj*blockDim.y+threadIdx.y] * x[j+xnum_cols*J];
        }
        y[j+xnum_cols*row] += sum;

      } // Loop over cols of x

    } // Loop over block cols of A
  }
}

template <typename Matrix, typename Vector2, typename Vector3>
void __spmm_csr_vector(const Matrix& A, const Vector2& x, Vector3& y,
                       cusp::row_major)
{
  typedef typename Matrix::index_type IndexType;
  typedef typename Matrix::value_type ValueType;
  typedef typename Matrix::memory_space MemorySpace;

  // Here we do a half-warp for columns of x to get a little better
  // granularity with the number of columns not being divisible by 32.
  // These numbers work out to 4 warps/block and 16 blocks/SM (based both on
  // shared memory and registers for Kepler).  This results in 100% occupancy.  
  const unsigned MAX_NNZ_PER_ROW = 32;
  const size_t COLS_PER_BLOCK = 16;
  const size_t ROWS_PER_BLOCK = 8;
  const size_t shared =
    ROWS_PER_BLOCK * MAX_NNZ_PER_ROW  * (sizeof(IndexType) + sizeof(ValueType));
  const size_t NUM_BLOCKS = (A.num_rows + ROWS_PER_BLOCK-1) / ROWS_PER_BLOCK;

  dim3 block(COLS_PER_BLOCK, ROWS_PER_BLOCK, 1);
  dim3 grid(NUM_BLOCKS, 1);

  spmm_csr_vector_kernel_row<IndexType, ValueType, MAX_NNZ_PER_ROW> <<<grid, block, shared>>>
    (A.num_rows, x.num_rows, x.num_cols,
     thrust::raw_pointer_cast(&A.row_offsets[0]),
     thrust::raw_pointer_cast(&A.column_indices[0]),
     thrust::raw_pointer_cast(&A.values[0]),
     thrust::raw_pointer_cast(&(x.values)[0]),
     thrust::raw_pointer_cast(&(y.values)[0]));
}

#endif

#if USE_CUSPARSE_COL

void __spmm_csr_vector(
  const csr_matrix<int,double,device_memory>& A,
  const array2d<double, device_memory, column_major>& x,
  array2d<double, device_memory, column_major>& y,
  cusp::column_major)
{
  CudaSparseSingleton & s = CudaSparseSingleton::singleton();
  const double alpha = 1 , beta = 0 ;
  cusparseStatus_t status =
    cusparseDcsrmm(s.handle,
                   CUSPARSE_OPERATION_NON_TRANSPOSE,
                   A.num_rows, x.num_cols, A.num_cols,
                   alpha,
                   s.descra,
                   thrust::raw_pointer_cast(&A.values[0]),
                   thrust::raw_pointer_cast(&A.row_offsets[0]),
                   thrust::raw_pointer_cast(&A.column_indices[0]),
                   thrust::raw_pointer_cast(&(x.values)[0]),
                   x.num_rows,
                   beta,
                   thrust::raw_pointer_cast(&(y.values)[0]),
                   y.num_rows);

  if ( CUSPARSE_STATUS_SUCCESS != status ) {
    throw std::runtime_error( std::string("ERROR - cusparseDcsrmv " ) );
  }
}

#else

#if 1

template <typename IndexType, typename ValueType, unsigned int VECTORS_PER_BLOCK, unsigned int THREADS_PER_VECTOR>
__launch_bounds__(VECTORS_PER_BLOCK * THREADS_PER_VECTOR,1)
__global__ void
spmm_csr_vector_kernel_col(const IndexType Anum_rows,
                           const IndexType xnum_rows,
                           const IndexType xnum_cols,
                           const IndexType * Ap,
                           const IndexType * Aj,
                           const ValueType * Ax,
                           const ValueType * x,
                           ValueType * y)
{
  __shared__ volatile ValueType sdata[VECTORS_PER_BLOCK * THREADS_PER_VECTOR + THREADS_PER_VECTOR / 2];  // padded to avoid reduction conditionals
  __shared__ volatile IndexType ptrs[VECTORS_PER_BLOCK][2];

  const IndexType THREADS_PER_BLOCK = VECTORS_PER_BLOCK * THREADS_PER_VECTOR;

  const IndexType thread_id   = THREADS_PER_BLOCK * blockIdx.x + threadIdx.x;    // global thread index
  const IndexType thread_lane = threadIdx.x & (THREADS_PER_VECTOR - 1);          // thread index within the vector
  const IndexType vector_id   = thread_id   /  THREADS_PER_VECTOR;               // global vector index
  const IndexType vector_lane = threadIdx.x /  THREADS_PER_VECTOR;               // vector index within the block
  const IndexType num_vectors = VECTORS_PER_BLOCK * gridDim.x;                   // total number of active vectors

  for(IndexType row = vector_id; row < Anum_rows; row += num_vectors)
  {
    // use two threads to fetch Ap[row] and Ap[row+1]
    // this is considerably faster than the straightforward version
    if(thread_lane < 2)
      ptrs[vector_lane][thread_lane] = Ap[row + thread_lane];

    const IndexType row_start = ptrs[vector_lane][0]; //same as: row_start = Ap[row];
    const IndexType row_end   = ptrs[vector_lane][1]; //same as: row_end   = Ap[row+1];

    //loop over cols of x
    for (IndexType j = 0; j < xnum_cols; j++){

      // initialize local sum
      ValueType sum = 0;

      if (THREADS_PER_VECTOR == 32 && row_end - row_start > 32)
      {
        // ensure aligned memory access to Aj and Ax

        IndexType jj = row_start - (row_start & (THREADS_PER_VECTOR - 1)) + thread_lane;

        // accumulate local sums
        if(jj >= row_start && jj < row_end)
          sum += Ax[jj] * x[Aj[jj]+xnum_rows*j];

        // accumulate local sums
        for(jj += THREADS_PER_VECTOR; jj < row_end; jj += THREADS_PER_VECTOR)
          sum += Ax[jj] * x[Aj[jj]+xnum_rows*j];
      }
      else
      {
        // accumulate local sums
        for(IndexType jj = row_start + thread_lane; jj < row_end; jj += THREADS_PER_VECTOR)
          sum += Ax[jj] * x[Aj[jj]+xnum_rows*j];
      }

      // store local sum in shared memory
      sdata[threadIdx.x] = sum;

      // reduce local sums to row sum
      if (THREADS_PER_VECTOR > 16) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x + 16];
      if (THREADS_PER_VECTOR >  8) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  8];
      if (THREADS_PER_VECTOR >  4) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  4];
      if (THREADS_PER_VECTOR >  2) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  2];
      if (THREADS_PER_VECTOR >  1) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  1];

      // first thread writes the result
      if (thread_lane == 0)
        y[j*Anum_rows+row] = sdata[threadIdx.x];

    }
  }
}

#else

template <typename IndexType, typename ValueType, unsigned int VECTORS_PER_BLOCK, unsigned int THREADS_PER_VECTOR>
__launch_bounds__(VECTORS_PER_BLOCK * THREADS_PER_VECTOR,1)
__global__ void
spmm_csr_vector_kernel_col(const IndexType Anum_rows,
                           const IndexType xnum_rows,
                           const IndexType xnum_cols,
                           const IndexType * Ar,
                           const IndexType * Ac,
                           const ValueType * Aval,
                           const ValueType * x,
                           ValueType * y)
{
  __shared__ volatile ValueType sdata[VECTORS_PER_BLOCK * THREADS_PER_VECTOR + THREADS_PER_VECTOR / 2];  // padded to avoid reduction conditionals
  __shared__ volatile IndexType ptrs[VECTORS_PER_BLOCK][2];

  extern __shared__ int sha[];
  ValueType * const sh_Aval = (ValueType*) sha;
  IndexType * const sh_Ac   = (IndexType*) (sh_Aval + 32 * VECTORS_PER_BLOCK);

  const IndexType THREADS_PER_BLOCK = VECTORS_PER_BLOCK * THREADS_PER_VECTOR;

  const IndexType thread_id   = THREADS_PER_BLOCK * blockIdx.x + threadIdx.x;    // global thread index
  const IndexType thread_lane = threadIdx.x & (THREADS_PER_VECTOR - 1);          // thread index within the vector
  const IndexType vector_id   = thread_id   /  THREADS_PER_VECTOR;               // global vector index i.e. global warp id
  const IndexType vector_lane = threadIdx.x /  THREADS_PER_VECTOR;               // vector index within the block
  const IndexType num_vectors = VECTORS_PER_BLOCK * gridDim.x;                   // total number of active vectors

  for(IndexType row = vector_id; row < Anum_rows; row += num_vectors) {
    // use two threads to fetch Ar[row] and Ar[row+1]
    // this is considerably faster than the straightforward version

    if(thread_lane < 2)
      ptrs[vector_lane][thread_lane] = Ar[row + thread_lane];

    const IndexType row_start = ptrs[vector_lane][0];                   //same as: row_start = Ar[row];
    const IndexType row_end   = ptrs[vector_lane][1];                   //same as: row_end   = Ar[row+1];
    //nnz in row
    const IndexType r = row_end - row_start;

    //fetch Aval and Ac for current row of A and store in shared mem
    for (IndexType i = thread_lane; i < r; i+= THREADS_PER_VECTOR){
      sh_Aval[vector_lane*32+i] = Aval[i+row_start];
      sh_Ac[vector_lane*32+i] = Ac[i+row_start];
    }
    //loop over cols of x
    for (IndexType j = 0; j < xnum_cols; j++){


      // initialize local sum

      ValueType sum = 0;

      if (THREADS_PER_VECTOR == 32 && row_end - row_start > 32)
      {
        // ensure aligned memory access to Ac and Aval
        IndexType jj = row_start - (row_start & (THREADS_PER_VECTOR - 1)) + thread_lane;
        // accumulate local sums
        if(jj >= row_start && jj < row_end)
          sum += Aval[jj] * x[j*xnum_rows+Ac[jj]];
        // accumulate local sums
        for(jj += THREADS_PER_VECTOR; jj < row_end; jj += THREADS_PER_VECTOR)
          sum += Aval[jj] * x[j*xnum_rows+Ac[jj]];
      }
      else
      {
        //accumulate local sums
        for (IndexType jj = thread_lane; jj < r; jj+=THREADS_PER_VECTOR)
          //sum += sh_Aval[jj+vector_lane*32] * x[j * xnum_rows + sh_Ac[vector_lane*32+jj]];
          sum += Aval[jj] * x[j*xnum_rows+Ac[jj]];

      }
      // store local sum in shared memory
      sdata[threadIdx.x] = sum;
      // reduce local sums to row sum
      if (THREADS_PER_VECTOR > 16) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x + 16];
      if (THREADS_PER_VECTOR >  8) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  8];
      if (THREADS_PER_VECTOR >  4) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  4];
      if (THREADS_PER_VECTOR >  2) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  2];
      if (THREADS_PER_VECTOR >  1) sdata[threadIdx.x] = sum = sum + sdata[threadIdx.x +  1];
      // first thread writes the result
      if (thread_lane == 0)
        y[j*Anum_rows+row] = sdata[threadIdx.x];
    }
  }
}

#endif

template <bool UseCache, unsigned int THREADS_PER_VECTOR,
          typename Matrix, typename Vector2, typename Vector3>
void __spmm_csr_vector_col(const Matrix& A, const Vector2& x, Vector3& y)
{
  typedef typename Matrix::index_type IndexType;
  typedef typename Matrix::value_type ValueType;

  const size_t THREADS_PER_BLOCK  = 256;
  const size_t VECTORS_PER_BLOCK  = THREADS_PER_BLOCK / THREADS_PER_VECTOR;
  //const size_t SHARED = VECTORS_PER_BLOCK * 32 * (sizeof(IndexType)+sizeof(ValueType));
  const size_t SHARED = 0;

  const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(spmm_csr_vector_kernel_col<IndexType, ValueType, VECTORS_PER_BLOCK, THREADS_PER_VECTOR>, THREADS_PER_BLOCK, SHARED);
  const size_t NUM_BLOCKS = std::min<size_t>(MAX_BLOCKS, DIVIDE_INTO(A.num_rows, VECTORS_PER_BLOCK));

  spmm_csr_vector_kernel_col<IndexType, ValueType, VECTORS_PER_BLOCK, THREADS_PER_VECTOR> <<<NUM_BLOCKS, THREADS_PER_BLOCK, SHARED>>>
    (A.num_rows, x.num_rows, x.num_cols,
     thrust::raw_pointer_cast(&A.row_offsets[0]),
     thrust::raw_pointer_cast(&A.column_indices[0]),
     thrust::raw_pointer_cast(&A.values[0]),
     thrust::raw_pointer_cast(&(x.values)[0]),
     thrust::raw_pointer_cast(&(y.values)[0]));
}

template <typename Matrix, typename Vector2, typename Vector3>
void __spmm_csr_vector(const Matrix& A, const Vector2& x, Vector3& y,
                       cusp::column_major)
{
#if 0
  typedef typename Vector2::index_type IndexType;
  for (IndexType col=0; col<x.num_cols; ++col) {
    multiply(A, x.column(col), y.column(col));
  }
#else
  typedef typename Matrix::index_type IndexType;

  const IndexType nnz_per_row = A.num_entries / A.num_rows;

  if (nnz_per_row <=  2) { __spmm_csr_vector_col<false, 2>(A, x, y); return; }
  if (nnz_per_row <=  4) { __spmm_csr_vector_col<false, 4>(A, x, y); return; }
  if (nnz_per_row <=  8) { __spmm_csr_vector_col<false, 8>(A, x, y); return; }
  if (nnz_per_row <= 16) { __spmm_csr_vector_col<false,16>(A, x, y); return; }

  __spmm_csr_vector_col<false,32>(A, x, y);
#endif
}

#endif

template <typename Matrix, typename Vector2, typename Vector3>
void spmm_csr_vector(const Matrix& A, const Vector2& x, Vector3& y)
{
  y.resize(A.num_rows, x.num_cols);
  __spmm_csr_vector(A, x, y, typename Vector2::orientation());
}

} // end namespace device
} // end namespace detail
} // end namespace cusp
