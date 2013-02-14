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

#ifndef STOKHOS_CUDA_BLOCKCRSMATRIX_HPP
#define STOKHOS_CUDA_BLOCKCRSMATRIX_HPP

#include <utility>
#include <sstream>
#include <stdexcept>

#include "KokkosArray_Cuda.hpp"
#include "Cuda/KokkosArray_Cuda_Parallel.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"

namespace Stokhos {

template< class BlockSpec , typename MatrixValue , typename VectorValue >
class Multiply<
  BlockCrsMatrix< BlockSpec , MatrixValue , KokkosArray::Cuda > ,
  KokkosArray::View< VectorValue** ,
                     KokkosArray::LayoutLeft ,
                     KokkosArray::Cuda > ,
  KokkosArray::View< VectorValue** ,
                     KokkosArray::LayoutLeft ,
                     KokkosArray::Cuda > >
{
public:

  typedef KokkosArray::Cuda                         device_type ;
  typedef device_type::size_type                    size_type ;
  typedef KokkosArray::View< VectorValue** ,KokkosArray::LayoutLeft , KokkosArray::Cuda > block_vector_type ;
  typedef BlockCrsMatrix< BlockSpec , MatrixValue , device_type >  matrix_type ;

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
    const size_type blockCount = m_A.graph.row_map.dimension(0) - 1 ;

    for ( size_type iBlock = blockIdx.x ;
                    iBlock < blockCount ; iBlock += gridDim.x ) {
      VectorValue y = 0 ;

      const size_type iEntryEnd = m_A.graph.row_map[iBlock+1];
            size_type iEntry    = m_A.graph.row_map[iBlock];

      for ( ; iEntry < iEntryEnd ; ++iEntry ) {
        const VectorValue * const x = & m_x( 0 , m_A.graph.entries(iEntry) );
        const MatrixValue * const a = & m_A.values( 0 , iEntry );

        y += Multiply< BlockSpec >::apply( m_A.block , a , x );
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
    const size_type thread_max =
      KokkosArray::Impl::cuda_internal_maximum_warp_count() * KokkosArray::Impl::CudaTraits::WarpSize ;

    const size_type row_count = A.graph.row_map.dimension(0) - 1 ;

    const dim3 grid(
      std::min( row_count , KokkosArray::Impl::cuda_internal_maximum_grid_count() ) , 1 , 1 );
    const dim3 block = Multiply<BlockSpec>::thread_block( A.block );

    const size_type shmem =
      Multiply<BlockSpec>::template shmem_size<block_vector_type>( A.block );

    if ( thread_max < block.x * block.y ) {
      std::ostringstream msg ;
      msg << "KokkosArray::Impl::Multiply< BlockCrsMatrix< Block , Value , Cuda > , ... >"
          << " ERROR: block dimension = " << block.x * block.y
          << " > " << thread_max << "== maximum Cuda threads per block" ;
      throw std::runtime_error(msg.str());
    }

    KokkosArray::Impl::cuda_parallel_launch_local_memory<<< grid , block , shmem >>>( Multiply(A,x,y) );
  }
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_BLOCKCRSMATRIX_HPP */

