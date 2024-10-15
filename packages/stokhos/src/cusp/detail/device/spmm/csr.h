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

#include <cusp/array1d.h>

#include <cusp/detail/format_utils.h>

//#include <iostream> 
//#include <cusp/print.h>
namespace cusp
{
namespace detail
{
namespace device
{

//Row-wise AX where X is a 2d array with one thread per row (based on scalar model)
template <typename IndexType,
          typename ValueType>
__global__ void
row_spmm_csr_scalar_kernel(const IndexType Anum_rows, 
			const IndexType xnum_row,
			const IndexType xnum_cols,
                       const IndexType * Ar,
                       const IndexType * Ac,
                       const ValueType * Aval,
                       const ValueType * x,
                             ValueType * y)
{
    const IndexType thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    const IndexType grid_size = gridDim.x * blockDim.x;


    for(IndexType row = thread_id; row < Anum_rows; row += grid_size)
    {
        const IndexType row_start = Ar[row];
        const IndexType row_end   = Ar[row+1];
	const IndexType r = row_end - row_start;	
        
	for (IndexType j = 0; j < xnum_cols; j++){

                ValueType sum = 0.0;
                for (IndexType jj = row_start; jj < row_end; jj++)
                    sum += Aval[jj] * x[j+xnum_cols*Ac[jj]];
		y[j+xnum_cols*row]=sum;
	}
        
    }
}



//Col-wise with one thread per row
template <typename IndexType,
          typename ValueType>
__global__ void
column_spmm_csr_scalar_kernel(const IndexType Anum_rows, 
			const IndexType xnum_rows,
                        const IndexType xnum_cols,
                       const IndexType * Ar,
                       const IndexType * Ac,
                       const ValueType * Aval,
                       const ValueType * x,
                             ValueType * y)
{
    const IndexType thread_id = blockDim.x * blockIdx.x + threadIdx.x;
    const IndexType grid_size = gridDim.x * blockDim.x;
    for(IndexType row = thread_id; row < Anum_rows; row += grid_size){
	const IndexType row_start = Ar[row];
        const IndexType row_end   = Ar[row+1];

	for (IndexType j = 0; j < xnum_cols; j++){
    	
                ValueType sum = 0;
                for (IndexType jj = row_start; jj < row_end; jj++)
                    sum += Aval[jj] * x[Ac[jj]+xnum_rows*j];//x(Ac[jj], j)
		y[j*Anum_rows+row]=sum;
    	}
    }
}


////////////////////////////////////////////////////////////////////////
//// CSR SpMM kernels based on a scalar model (one thread per row)
/////////////////////////////////////////////////////////////////////////
////
//// spmm_csr_scalar_device
////   Straightforward translation of standard CSR SpMV to CUDA
////   where each thread computes y[i] = A[i,:] * x 
////   (the dot product of the i-th row of A with the x vector)


template <typename Matrix1,
          typename Vector2,
          typename Vector3>
void __spmm_csr_scalar(const Matrix1& A,
              const Vector2& x,
                    Vector3& y, cusp::row_major)
{
	CUSP_PROFILE_SCOPED();
	typedef typename Vector3::index_type   IndexType;
	typedef typename Vector3::value_type   ValueType;
	typedef typename Vector3::memory_space MemorySpace;
	const size_t BLOCK_SIZE = 256;
	const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(row_spmm_csr_scalar_kernel<IndexType, ValueType>, BLOCK_SIZE, (size_t) 0);
	const size_t NUM_BLOCKS = std::min(MAX_BLOCKS, DIVIDE_INTO(A.num_rows, BLOCK_SIZE));



	row_spmm_csr_scalar_kernel<IndexType,ValueType> <<<NUM_BLOCKS, BLOCK_SIZE >>>
        (A.num_rows, x.num_rows, x.num_cols, 
         thrust::raw_pointer_cast(&A.row_offsets[0]),
         thrust::raw_pointer_cast(&A.column_indices[0]),
         thrust::raw_pointer_cast(&A.values[0]),
         thrust::raw_pointer_cast(&(x.values)[0]),
	 thrust::raw_pointer_cast(&(y.values)[0]));
	
}

template <typename Matrix1,
          typename Vector2,
          typename Vector3>
void __spmm_csr_scalar(const Matrix1& A,
              const Vector2& x,
                    Vector3& y, cusp::column_major)
{
        CUSP_PROFILE_SCOPED();
        typedef typename Vector3::index_type   IndexType;
        typedef typename Vector3::value_type   ValueType;
        typedef typename Vector3::memory_space MemorySpace;
        const size_t BLOCK_SIZE = 256;
        const size_t MAX_BLOCKS = cusp::detail::device::arch::max_active_blocks(column_spmm_csr_scalar_kernel<IndexType, ValueType>, BLOCK_SIZE, (size_t) 0);
        const size_t NUM_BLOCKS = std::min(MAX_BLOCKS, DIVIDE_INTO(A.num_rows, BLOCK_SIZE));
        column_spmm_csr_scalar_kernel<IndexType,ValueType> <<<NUM_BLOCKS, BLOCK_SIZE>>>
        (A.num_rows, x.num_rows, x.num_cols,
         thrust::raw_pointer_cast(&A.row_offsets[0]),
         thrust::raw_pointer_cast(&A.column_indices[0]),
         thrust::raw_pointer_cast(&A.values[0]),
         thrust::raw_pointer_cast(&(x.values)[0]),
         thrust::raw_pointer_cast(&(y.values)[0]));

}


template <typename Matrix1,
          typename Vector2,
          typename Vector3>
void spmm_csr_scalar(const Matrix1& A,
              const Vector2& x,
                    Vector3& y)
{
//Determine if row-wise or col-wise then call appropriate multiply
 __spmm_csr_scalar(A, x, y, typename Vector2::orientation());
}
} // end namespace device
} // end namespace detail
} // end namespace cusp
