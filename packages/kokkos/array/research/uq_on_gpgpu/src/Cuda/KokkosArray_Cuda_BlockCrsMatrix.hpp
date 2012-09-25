/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_CUDA_BLOCKCRSMATRIX_HPP
#define KOKKOSARRAY_CUDA_BLOCKCRSMATRIX_HPP

#include <utility>
#include <sstream>
#include <stdexcept>

namespace KokkosArray {
namespace Impl {

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
  typedef KokkosArray::Cuda                              device_type ;
  typedef device_type::size_type                    size_type ;
  typedef View< VectorValue** , LayoutLeft , Cuda > block_vector_type ;
  typedef BlockCrsMatrix< BlockSpec , MatrixValue , device_type >  matrix_type ;
  typedef Impl::Multiply< BlockSpec , void , void > block_matrix_type ;

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
      cuda_internal_maximum_warp_count() * Impl::CudaTraits::WarpSize ;

    const size_type row_count =A.graph.row_map.dimension(0) - 1 ;

    const dim3 grid(
      std::min( row_count , cuda_internal_maximum_grid_count() ) , 1 , 1 );
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

    Impl::cuda_parallel_launch_local_memory<<< grid , block , shmem >>>( Multiply(A,x,y) );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_CUDA_BLOCKCRSMATRIX_HPP */

