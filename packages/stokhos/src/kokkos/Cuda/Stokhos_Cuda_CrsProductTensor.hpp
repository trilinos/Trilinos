// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_CRSPRODUCTTENSOR_HPP
#define STOKHOS_CUDA_CRSPRODUCTTENSOR_HPP

#include <iostream>

#include "KokkosArray_Cuda.hpp"
#include "Cuda/KokkosArray_Cuda_Parallel.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_CrsProductTensor.hpp"

#include "cuda_profiler_api.h"

namespace Stokhos {

//----------------------------------------------------------------------------

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< CrsProductTensor< TensorScalar, KokkosArray::Cuda >,
                  MatrixScalar, KokkosArray::Cuda >,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  DefaultSparseMatOps >
{
public:
  
  typedef KokkosArray::Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef CrsProductTensor< TensorScalar , device_type >       tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef KokkosArray::View< VectorScalar** , KokkosArray::LayoutLeft , KokkosArray::Cuda >           vector_type ;

  class ProductTensorLoop {
  public:

    const matrix_type m_A ;
    const vector_type m_x ;
    const vector_type m_y ;
    const size_type BlockSize ;

    ProductTensorLoop( const matrix_type & A ,
                       const vector_type & x ,
                       const vector_type & y ,
		       const size_type block_size )
      : m_A( A ), m_x( x ), m_y( y ), BlockSize(block_size) {}

    __device__
    void operator()(void) const
    {
      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();
      // const size_type rem = dim % 32;
      // const size_type dim_align = rem > 0 ? dim + rem : dim ;
      const size_type dim_align = dim;

      VectorScalar * const sh_x = 
	kokkos_impl_cuda_shared_memory<VectorScalar>();
      VectorScalar * const sh_A = sh_x + BlockSize*dim_align ;
      VectorScalar * const sh_y = sh_A + BlockSize*dim_align ;
      volatile VectorScalar * const sh_t = sh_y + dim_align ;

      const size_type nid = 
	KokkosArray::Impl::CudaTraits::WarpSize * blockDim.y ;
      const size_type tid = 
	threadIdx.x + KokkosArray::Impl::CudaTraits::WarpSize * threadIdx.y ;

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
      size_type numBlock = (iBlockEntryEnd-iBlockEntryBeg) / BlockSize;
      const size_type remBlock = (iBlockEntryEnd-iBlockEntryBeg) % BlockSize;
      if (remBlock > 0) ++numBlock;

      // Zero y
      for ( size_type i = tid ; i < dim ; i += nid ) {
	sh_y[i] = 0.0;
      }

      // Loop over columns in the discrete (finite element) system.

      size_type iBlockEntry = iBlockEntryBeg ;
      for ( size_type block = 0 ; block < numBlock ; ++block ) {
	const size_type block_size = 
	  (block == numBlock-1 && remBlock > 0) ? remBlock : BlockSize;

	// Coalesced read blocks of X and A into shared memory
	for ( size_type iBlock = 0; iBlock < block_size ; ++iBlock ) {

	  const size_type iBlockColumn = 
	    m_A.graph.entries( iBlockEntry + iBlock );

	  const VectorScalar * const x = & m_x(        0 , iBlockColumn );
	  const MatrixScalar * const A = & m_A.values( 0 , iBlockEntry + iBlock );

	  // Wait for X and A to be used in the previous iteration before reading new values.
        __syncthreads();

	  // Coalesced read by the whole block from global memory:
	  for ( size_type i = tid ; i < dim ; i += nid ) {
	    sh_x[iBlock+i*block_size] = x[i] ; // m_x(        i , iBlockColumn );
	    sh_A[iBlock+i*block_size] = A[i] ; // m_A.values( i , iBlockEntry );
	  }

	}

	__syncthreads(); // wait for X and A to be read

	// This cuda block is responsible for computing all values of 'y'
	//
	for ( size_type iyInner = threadIdx.y ; iyInner < dim ; 
	      iyInner += blockDim.y) {

	  VectorScalar y = 0 ;

	  // Product tensor entries which this warp will iterate:
	  //
	  const size_type iBeg = m_A.block.entry_begin( iyInner ) ;
	  const size_type iEnd = m_A.block.entry_end(   iyInner ) ;

	  // Loop through sparse tensor contributions with coalesced reads.

	  for ( size_type i = iBeg + threadIdx.x ; i < iEnd ; 
		i += blockDim.x ) {

	    // Read 'CudaTraits::WarpSize' entries from the tensor
	    const int j = m_A.block.coord( i , 0 ); // coalesced read
	    const int k = m_A.block.coord( i , 1 ); // coalesced read
	    const MatrixScalar v = m_A.block.value(i);    // coalesced read

	    for ( size_type iBlock = 0; iBlock < block_size ; ++iBlock ) {
	      const size_type jj = iBlock+j*block_size;
	      const size_type kk = iBlock+k*block_size;
	      y += v * ( sh_A[jj] * sh_x[kk] + sh_A[kk] * sh_x[jj] ) ;
	    }

	  }

	  //__syncthreads(); // wait for X and A to be used.

	  // Reduction of 'y' within 'CudaTraits::WarpSize'

	  sh_t[ tid ] = y ;

	  if ( threadIdx.x + 16 < KokkosArray::Impl::CudaTraits::WarpSize ) sh_t[tid] += sh_t[tid+16];
	  if ( threadIdx.x +  8 < KokkosArray::Impl::CudaTraits::WarpSize ) sh_t[tid] += sh_t[tid+ 8];
	  if ( threadIdx.x +  4 < KokkosArray::Impl::CudaTraits::WarpSize ) sh_t[tid] += sh_t[tid+ 4];
	  if ( threadIdx.x +  2 < KokkosArray::Impl::CudaTraits::WarpSize ) sh_t[tid] += sh_t[tid+ 2];
	  if ( threadIdx.x +  1 < KokkosArray::Impl::CudaTraits::WarpSize ) sh_t[tid] += sh_t[tid+ 1];

	  if ( threadIdx.x == 0 )
	    sh_y[iyInner] += sh_t[tid];
	  
	}

	iBlockEntry += block_size;
      }

      __syncthreads();
	
      // Store result back in global memory
      for ( size_type i = tid ; i < dim ; i += nid ) {
	m_y( i , blockIdx.x ) = sh_y[ i ] ;
      }
    }
  };

  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_type row_count = A.graph.row_map.dimension(0) - 1 ;
    const size_type tensor_dimension = A.block.dimension();
    // const size_type rem = tensor_dimension % 32;
    // const size_type tensor_align = 
    //   rem > 0 ? tensor_dimension + rem : tensor_dimension ;
    const size_type tensor_align = tensor_dimension ;

    const size_type nWarp = 16 ; 
    const dim3 dBlock( KokkosArray::Impl::CudaTraits::WarpSize , nWarp , 1 );
    const dim3 dGrid( row_count , 1 , 1 );

    // Use at most half of shared memory to get 2 blocks per SMP
    const size_type shcap = KokkosArray::Impl::CudaTraits::SharedMemoryCapacity / 2;
    int block_size = ((shcap / sizeof(VectorScalar) - dBlock.x*dBlock.y) / tensor_align - 1) / 2;
    block_size = std::min( block_size, 9 );
    if (block_size % 2 == 0)
      --block_size;
    //const int block_size = 5;
    const size_type shmem = 
      sizeof(VectorScalar) * ((2*block_size+1) * tensor_align + dBlock.x*dBlock.y);

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl 
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
	      << "  block_size(" << block_size << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tensor_entry_max(" << A.block.entry_maximum() << ")" << std::endl
              ;
#endif
    //cudaProfilerStart();
    KokkosArray:: Impl::cuda_parallel_launch_local_memory<<< dGrid , dBlock , shmem >>>
      ( ProductTensorLoop( A , x , y, block_size ) );
    //cudaProfilerStop();
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_PRODUCTTENSOR_HPP */

