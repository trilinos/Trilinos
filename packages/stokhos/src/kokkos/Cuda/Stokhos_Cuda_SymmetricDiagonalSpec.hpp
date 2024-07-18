// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_CUDA_SYMMETRIC_DIAGONAL_SPEC_HPP
#define KOKKOS_CUDA_SYMMETRIC_DIAGONAL_SPEC_HPP

#include <stdexcept>

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_SymmetricDiagonalSpec.hpp"

//----------------------------------------------------------------------------

namespace Stokhos {

template<>
class BlockMultiply< SymmetricDiagonalSpec< Kokkos::Cuda > >
{
public:
  typedef Kokkos::Cuda execution_space ;
  typedef execution_space::size_type size_type ;
  typedef SymmetricDiagonalSpec< Kokkos::Cuda >  block_type ;

  __host__
  static dim3 thread_block( const block_type & block )
  {
    const int d = block.dimension();
    const int warp_size = Kokkos::Impl::CudaTraits::WarpSize;
    auto const maxWarpCount = std::min<unsigned>(
        execution_space().cuda_device_prop().maxThreadsPerBlock / warp_size,
        warp_size);
    const int y = ( maxWarpCount * warp_size ) / d ;

    if ( 0 == y ) {
      throw std::runtime_error( std::string("Stokhos::Multiply< SymmetricDiagonalSpec<Cuda> > ERROR: block too large") );
    }

    // dimension X #diagonals to concurrently process
    return dim3( d , std::min( y , ( 1 + d ) / 2 ) , 1 );
  }

  template< typename VectorValue >
  __host__
  static size_type shmem_size( const block_type & block )
  {
    const dim3 d = thread_block( block );

    return sizeof(VectorValue) * d.x * d.y ;
  }

  __host__
  static size_type matrix_size( const block_type & block )
    { return block.matrix_size(); }

  // Required: blockDim.x == block.dimension()
  // Required: blockDim.y > 1
  // Computing: Y[ threadIdx.x ]
  template< typename MatrixValue , typename VectorValue >
  __device__
  static VectorValue apply( const block_type  & block ,
                            const MatrixValue * const a ,
                            const VectorValue * const x )
  {
    const int  dimension = block.dimension();
    const int  dim_half  = ( dimension + 1 ) >> 1 ;

    VectorValue * const shX = kokkos_impl_cuda_shared_memory<VectorValue>();

    int ia = -1 ;
    int ix = -1 ;

    VectorValue y = 0 ;

    if ( 0 == threadIdx.y ) {
      // Load vector and multiply the main diagonal (first diagonal)

      shX[ threadIdx.x ] = x[ threadIdx.x ]; // Load 'x'

      y = shX[ threadIdx.x ] * a[ threadIdx.x ];
    }

    __syncthreads(); // Wait for load of 'x'

    if ( 0 == threadIdx.y && ! ( dimension & 01 ) ) {

      // If even number of diagonals then the last diagonal is half-length
      // and only used once.

      ix = threadIdx.x ;
      ia = threadIdx.x + dim_half * dimension ;

      if ( threadIdx.x < dim_half ) {
        ix += dim_half ;
      }
      else {
        ix -= dim_half ;
        ia -= dim_half ;
      }
      y += shX[ ix ] * a[ ia ];
    }

    //------------------------------------

    const int A_stride = blockDim.y * dimension ;

    int d = 1 + threadIdx.y ;

    const MatrixValue * A = a + d * dimension ;

    for ( ; d < dim_half ; d += blockDim.y , A += A_stride ) {

      ix = threadIdx.x + d ; if ( dimension <= ix ) ix -= dimension ;
      ia = threadIdx.x - d ; if ( ia < 0 ) ia += dimension ;

      // A 'threadIdx.y' group accesses the matrix diagonal
      // A[ 0 .. dimension - 1 ] twice in the following statement.
      // Spatial and temporal locality are excellent for L1 cache reuse
      // for the second access.

      y += shX[ ix ] * A[ threadIdx.x ] +
           shX[ ia ] * A[ ia ];
    }

    if ( 0 < threadIdx.y ) {
      shX[ threadIdx.x + threadIdx.y * dimension ] = y ;
    }

    __syncthreads(); // Wait for load of contributions to 'y'

    for ( ix = 1 ; ix < blockDim.y ; ++ix ) {
      y += shX[ threadIdx.x + ix * dimension ];
    }

    return y ;
  }
};

} /* namespace Stokhos */

//----------------------------------------------------------------------------

#endif /* #ifndef STOKHOS_CUDA_SYMMETRIC_DIAGONAL_SPEC_HPP */
