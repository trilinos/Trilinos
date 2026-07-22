// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_BLOCKCRSMATRIX_HPP
#define STOKHOS_CUDA_BLOCKCRSMATRIX_HPP

#include <utility>
#include <sstream>
#include <stdexcept>

#include "Kokkos_Core.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"

namespace Stokhos {

template< class BlockSpec , typename MatrixValue , typename VectorValue >
class Multiply<
  BlockCrsMatrix< BlockSpec , MatrixValue , Kokkos::Cuda > ,
  Kokkos::View< VectorValue** , Kokkos::LayoutLeft , Kokkos::Cuda > ,
  Kokkos::View< VectorValue** , Kokkos::LayoutLeft , Kokkos::Cuda > >
{
public:

  typedef Kokkos::Cuda                         execution_space ;
  typedef execution_space::size_type                    size_type ;
  typedef Kokkos::View< VectorValue** ,Kokkos::LayoutLeft , Kokkos::Cuda > block_vector_type ;
  typedef BlockCrsMatrix< BlockSpec , MatrixValue , execution_space >  matrix_type ;

  const matrix_type  m_A ;
  const block_vector_type  m_x ;
  const block_vector_type  m_y ;

  Multiply( const matrix_type & A ,
            const block_vector_type & x ,
            const block_vector_type & y )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  {}

  //--------------------------------------------------------------------------
  //  A( storage_size( m_A.block.size() ) , m_A.graph.row_map.size() );
  //  x( m_A.block.dimension() , m_A.graph.row_map.first_count() );
  //  y( m_A.block.dimension() , m_A.graph.row_map.first_count() );
  //

  __device__
  void operator()(void) const
  {
    const size_type blockCount = m_A.graph.row_map.extent(0) - 1 ;

    for ( size_type iBlock = blockIdx.x ;
                    iBlock < blockCount ; iBlock += gridDim.x ) {
      VectorValue y = 0 ;

      const size_type iEntryEnd = m_A.graph.row_map[iBlock+1];
            size_type iEntry    = m_A.graph.row_map[iBlock];

      for ( ; iEntry < iEntryEnd ; ++iEntry ) {
        const VectorValue * const x = & m_x( 0 , m_A.graph.entries(iEntry) );
        const MatrixValue * const a = & m_A.values( 0 , iEntry );

        y += BlockMultiply< BlockSpec >::apply( m_A.block , a , x );
      }

      if ( threadIdx.x + blockDim.x * threadIdx.y < m_A.block.dimension() ) {
        m_y(threadIdx.x,iBlock) = y ;
      }
    }
  }

  static void apply( const matrix_type & A ,
                     const block_vector_type & x ,
                     const block_vector_type & y )
  {
    const int warp_size = Kokkos::Impl::CudaTraits::WarpSize;
    auto const maxWarpCount = std::min<unsigned>(
        execution_space().cuda_device_prop().maxThreadsPerBlock / warp_size,
        warp_size);

    const size_type thread_max =
      maxWarpCount * warp_size ;

    const size_type row_count = A.graph.row_map.extent(0) - 1 ;

    const size_type maxGridSizeX = execution_space().cuda_device_prop().maxGridSize[0];
    const dim3 grid(
      std::min( row_count , maxGridSizeX ) , 1 , 1 );
    const dim3 block = BlockMultiply<BlockSpec>::thread_block( A.block );

    const size_type shmem =
      BlockMultiply<BlockSpec>::template shmem_size<block_vector_type>( A.block );

    if ( thread_max < block.x * block.y ) {
      std::ostringstream msg ;
      msg << "Kokkos::Impl::Multiply< BlockCrsMatrix< Block , Value , Cuda > , ... >"
          << " ERROR: block dimension = " << block.x * block.y
          << " > " << thread_max << "== maximum Cuda threads per block" ;
      throw std::runtime_error(msg.str());
    }

    Kokkos::Impl::cuda_parallel_launch_local_memory<<< grid , block , shmem >>>( Multiply(A,x,y) );
  }
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_BLOCKCRSMATRIX_HPP */
