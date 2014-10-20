// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
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

  typedef Kokkos::Cuda                         device_type ;
  typedef device_type::size_type                    size_type ;
  typedef Kokkos::View< VectorValue** ,Kokkos::LayoutLeft , Kokkos::Cuda > block_vector_type ;
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
    const size_type blockCount = m_A.graph.row_map.dimension_0() - 1 ;

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
    const size_type thread_max =
      Kokkos::Impl::cuda_internal_maximum_warp_count() * Kokkos::Impl::CudaTraits::WarpSize ;

    const size_type row_count = A.graph.row_map.dimension_0() - 1 ;

    const dim3 grid(
      std::min( row_count , Kokkos::Impl::cuda_internal_maximum_grid_count() ) , 1 , 1 );
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
