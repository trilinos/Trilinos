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


namespace cusp
{
namespace detail
{
namespace device
{

template <typename IndexType,
          typename ValueType>
__global__ void
_spmm_csr_vector_kernel(const int nnz, const IndexType Anum_rows,
                        const IndexType xnum_rows,
                        const IndexType xnum_cols,
                       const IndexType * Ar,
                       const IndexType * Ac,
                       const ValueType * Aval,
                       const ValueType * x,
                             ValueType * y)
{
   
    //allocate storage for shared A row 
    extern __shared__ int sh[];
    volatile ValueType * const sh_Aval = (ValueType*) sh;			
    volatile IndexType * const sh_Ac   = (IndexType*) (sh_Aval + nnz * blockDim.z);

    for (IndexType row = threadIdx.z + blockDim.z * blockIdx.x; row < Anum_rows; row += gridDim.x * blockDim.z){
	const IndexType row_start = Ar[row];
        const IndexType row_end   = Ar[row+1];
	const IndexType r = row_end - row_start;

	__syncthreads();
	//read row of A into shared mem using all threads in block
	for (IndexType i = threadIdx.x + blockDim.x * threadIdx.y; i < r; i+=blockDim.x * blockDim.y){
		sh_Aval[i + nnz*threadIdx.z] = Aval[i+row_start];
		sh_Ac[i+ nnz*threadIdx.z] = Ac[i+row_start];
	}
	__syncthreads();
	//loop over cols of x
	for (IndexType j = threadIdx.x + blockDim.x * threadIdx.y; j < xnum_cols; j+=blockDim.x * blockDim.y){
		ValueType sum = 0.0;
      		for (IndexType jj = 0; jj < r; jj++){
			IndexType J = sh_Ac[jj+nnz*threadIdx.z];
			sum +=  sh_Aval[jj+ nnz*threadIdx.z]* x[j+xnum_cols*J];
		}
		y[j+xnum_cols*row]= sum; 
	}
   }
	

}

 

template <typename IndexType, typename ValueType, unsigned int VECTORS_PER_BLOCK, unsigned int THREADS_PER_VECTOR>
__launch_bounds__(VECTORS_PER_BLOCK * THREADS_PER_VECTOR,1)
__global__ void
column_spmm_csr_vector_kernel(const IndexType Anum_rows,
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

 for(IndexType row = vector_id; row < Anum_rows; row += num_vectors)
    {
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
			sum += sh_Aval[jj+vector_lane*32] * x[j * xnum_rows + sh_Ac[vector_lane*32+jj]];

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


template <typename Matrix, typename Vector2,
          typename Vector3>
void __spmm_csr_vector(const int nnz, const Matrix&    A,  
                       const Vector2& x, 
                             Vector3& y,
                        cusp::row_major)
{
	typedef typename Matrix::index_type IndexType;
	typedef typename Matrix::value_type ValueType;
	typedef typename Matrix::memory_space MemorySpace;
	
	
	const size_t BLOCK_SIZE =32;
	const size_t ROWS_PER_BLOCK = 7;
	dim3 block(BLOCK_SIZE, ROWS_PER_BLOCK,2);
	const size_t shared =  block.z * nnz * (sizeof(IndexType)+sizeof(ValueType));


	const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(_spmm_csr_vector_kernel<IndexType, ValueType>, BLOCK_SIZE, shared);
	const size_t NUM_BLOCKS = std::min( MAX_BLOCKS, DIVIDE_INTO(A.num_rows, ROWS_PER_BLOCK));
	
	dim3 grid(NUM_BLOCKS, 1);

        cudaProfilerStart();
	
	_spmm_csr_vector_kernel<IndexType, ValueType> <<<grid, block, shared>>> 
        (nnz, A.num_rows, x.num_rows, x.num_cols,
         thrust::raw_pointer_cast(&A.row_offsets[0]),
         thrust::raw_pointer_cast(&A.column_indices[0]),
         thrust::raw_pointer_cast(&A.values[0]),
         thrust::raw_pointer_cast(&(x.values)[0]),
         thrust::raw_pointer_cast(&(y.values)[0]));
	
	cudaProfilerStop();
}

template <unsigned int THREADS_PER_VECTOR, typename Matrix, typename Vector2,
          typename Vector3>
void __spmm_csr_vector(const Matrix&    A, const size_t nnz,
                       const Vector2& x,
                             Vector3& y,
                        cusp::column_major)
{
    typedef typename Matrix::index_type IndexType;
    typedef typename Matrix::value_type ValueType;

    const size_t THREADS_PER_BLOCK  = 128;
    const size_t VECTORS_PER_BLOCK  = THREADS_PER_BLOCK / THREADS_PER_VECTOR;

    const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(column_spmm_csr_vector_kernel<IndexType, ValueType, VECTORS_PER_BLOCK, THREADS_PER_VECTOR>, THREADS_PER_BLOCK, (size_t) 0);
    const size_t NUM_BLOCKS = std::min<size_t>(MAX_BLOCKS, DIVIDE_INTO(A.num_rows, VECTORS_PER_BLOCK));
    const size_t SHARED = VECTORS_PER_BLOCK * 32 * (sizeof(IndexType)+sizeof(ValueType));    


    column_spmm_csr_vector_kernel<IndexType, ValueType, VECTORS_PER_BLOCK, THREADS_PER_VECTOR> <<<NUM_BLOCKS, THREADS_PER_BLOCK, SHARED>>>
        (A.num_rows, x.num_rows, x.num_cols,
         thrust::raw_pointer_cast(&A.row_offsets[0]),
         thrust::raw_pointer_cast(&A.column_indices[0]),
         thrust::raw_pointer_cast(&A.values[0]),
         thrust::raw_pointer_cast(&(x.values)[0]),
         thrust::raw_pointer_cast(&(y.values)[0]));



}

template <typename IndexType>
__global__ void 
get_nnz_kernel(const IndexType Anum_rows, const IndexType * Ar, IndexType * nnz)
{
	IndexType r = 0; 
	IndexType max = 1; 
	for (IndexType i = 0; i < Anum_rows; i++){
		r = Ar[i+1] - Ar[i];
		if (r > max)
			max = r;
	}
	
	nnz[0] =  max;
}

template <typename Matrix,
          typename Vector2,
          typename Vector3>
void spmm_csr_vector(const Matrix&    A, 
			const Vector2& x,
			 Vector3& y)
{
    typedef typename Matrix::index_type IndexType;
    typedef typename Matrix::memory_space MemorySpace;
    y.resize(A.num_rows, x.num_cols);
    dim3 block(1);
    dim3 grid(1);

    cusp::array1d<IndexType, MemorySpace> nnz(1);
    get_nnz_kernel<IndexType> <<<block, grid>>>
        (A.num_rows, thrust::raw_pointer_cast(&A.row_offsets[0]), thrust::raw_pointer_cast(&nnz[0]));
    const IndexType nnz_per_row = nnz[0];
     __spmm_csr_vector(nnz_per_row, A, x, y, typename Vector2::orientation()); 
}


} // end namespace device
} // end namespace detail
} // end namespace cusp

