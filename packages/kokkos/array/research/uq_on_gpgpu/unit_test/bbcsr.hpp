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

#ifndef KOKKOSARRAY_BBCRS_HPP
#define KOKKOSARRAY_BBCRS_HPP

#include <stdexcept>
#include <sstream>

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief  Big Block CRS Graph
 *
 *  CRS Matrix storage with big NxN blocks instead of scalars.
 *  
 *  block_system_size :
 *    Number of row blocks, equal to number of column blocks
 *  
 *  block_size :
 *    Each dense block is block_size X block_size
 *    Thus N = block_system_size * block_size 
 *    and the system is N X N
 *  
 *  block_column_offset[ 0 .. system_size ] :
 *    Offsets into the block_column_index array for each row-block
 *
 *  block_column_index[ 0 .. block_column_offset[ system_size ] - 1 ] :
 *    Index into the matrix coefficent array for the beginning
 *    of the block.
 */
template< class Device >
class BigBlockCRSGraph {
public:
  typedef Device                     device_type ;
  typedef typename Device::size_type index_type ;

  typedef View< index_type[] , device_type > vector_type ;

  /** \brief  Index to the begining of the block */
  index_type offset_block( index_type block_row , index_type block_column ) const ;

  /** \brief  Index to a particular entry */
  index_type offset_entry( index_type row , index_type column ) const
  {
    const int block_offset = offset_block( row / block_size , column / block_size );
    // Column major ordering
    return block_offset + row % block_size + block_size * ( column % block_size );
  }

  index_type  block_system_size ;
  index_type  block_size ;
  index_type  block_stride ;
  index_type  block_maximum_columns ; ///< Maximum number of block-columns in a block-row
  vector_type block_column_offset ;
  vector_type block_column_index ;
};

//----------------------------------------------------------------------------

template< class Device , typename MatrixScalar , typename VectorScalar >
class BigBlockCRSMultiply ;

template< class Device , typename MatrixScalar , typename VectorScalar >
void multiply( const BigBlockCRSGraph<Device> & A_graph ,
               const View<MatrixScalar[],Device> & A_coeff ,
               const View<VectorScalar[],Device> & input ,
               const View<VectorScalar[],Device> & output )
{
  BigBlockCRSMultiply< Device , MatrixScalar , VectorScalar >
    ( A_graph , A_coeff , input , output );
}

//----------------------------------------------------------------------------

template< typename MatrixScalar , typename VectorScalar >
class BigBlockCRSMultiply< Cuda , MatrixScalar , VectorScalar > {
public:
  typedef Cuda::size_type                        index_type ;
  typedef BigBlockCRSGraph< Cuda >               graph_type ;
  typedef View< MatrixScalar[] , Cuda > matrix_type ;
  typedef View< VectorScalar[] , Cuda > vector_type ;

  const graph_type  m_graph ;
  const matrix_type m_matrix ;
  const vector_type m_input ;
  const vector_type m_output ;
  const index_type  m_shared_index_offset ;

  BigBlockCRSMultiply( const graph_type  & graph ,
                       const matrix_type & matrix ,
                       const vector_type & input ,
                       const vector_type & output )
  : m_graph( graph )
  , m_matrix( matrix )
  , m_input(  input )
  , m_output( output )
  , m_shared_index_offset( ( graph.block_size * sizeof(VectorScalar) ) / sizeof(index_type) )
  {
    // Require: 2 < m_shared_index_offset
    // Require: m_graph.block_size <= blockDim.x
    // Shared memory use = sizeof(VectorScalar) * m_graph.block_size +
    //                     sizeof(index_type) *   m_graph.block_maximum_columns ;

    if ( m_shared_index_offset < 2 ) {
      throw std::runtime_error( std::string("block is too small") );
    }

    const index_type warp_count = m_graph.block_size / Impl::CudaTraits::WarpSize +
                                ( m_graph.block_size % Impl::CudaTraits::WarpSize ? 1 : 0 );

    if ( Cuda::maximum_warp_count() < warp_count ) {
      throw std::runtime_error( std::string("block is too big") );
    }

/*
    if ( Cuda::maximum_grid_count() < m_graph.block_system_size ) {
      std::ostringstream msg ;
      msg << "Too many CUDA blocks: max = "
          << Cuda::maximum_grid_count()
          << " , requested = " << m_graph.block_system_size ;

      throw std::runtime_error( msg.str() );
    }
*/

    const index_type thread_block_count = m_graph.block_system_size ;
    const index_type thread_per_block   = warp_count * Impl::CudaTraits::WarpSize ;
    const index_type shmem_size = sizeof(VectorScalar) * m_graph.block_size +
                                  sizeof(index_type) *   m_graph.block_maximum_columns ;

    Impl::cuda_parallel_launch_local_memory<<< thread_block_count , thread_per_block , shmem_size >>>( *this );
  }

  __device__
  void operator()(void) const
  {
    extern __shared__ Cuda::size_type shared_data[] ;

    VectorScalar * const shared_input = (VectorScalar *) shared_data ;

    // Required: 2 < m_graph.block_size <= blockDim.x

    const index_type this_thread_has_work = threadIdx.x < m_graph.block_size ;

    // Read the span of indices for this block:
    if ( threadIdx.x < 2 ) {
      shared_data[threadIdx.x] =
        m_graph.block_column_offset( blockIdx.x + threadIdx.x );
    }

    __syncthreads();

    const index_type block_begin = shared_data[0] ;
    const index_type block_count = shared_data[1] - block_begin ;

    // Read the starting indices for all of the row-blocks' coefficients
    for ( index_type i = threadIdx.x ; i < block_count ; i += blockDim.x ) {
      shared_data[ m_shared_index_offset + i ] =
        m_graph.block_column_index( block_begin + i ) * m_graph.block_size ;
    }

    __syncthreads();

    // Loop over matrix entries
    index_type matrix_offset = m_graph.block_stride *
                               m_graph.block_size * block_begin ;

    VectorScalar y = 0 ;

    for ( index_type i = 0 ; i < block_count ; ++i ) {

      // Load the input vector for this block
      if ( this_thread_has_work ) {
        // Pre-loaded from global memory into shared memory:
        //   shared_data[ m_shared_index_offset + i ] == m_graph.block_column_index( block_begin + i )
        const index_type vector_offset = shared_data[ m_shared_index_offset + i ];
        shared_input[ threadIdx.x ] = m_input( vector_offset + threadIdx.x );
      }

      __syncthreads();

      if ( this_thread_has_work ) {
        // This thread performs its row of the dense block multiplication
        for ( index_type k = 0 ; k < m_graph.block_size ;
              ++k , matrix_offset += m_graph.block_stride ) {

          y += m_matrix( matrix_offset + threadIdx.x ) * shared_input[k] ;
        }
      }

      __syncthreads();
    }

    if ( this_thread_has_work ) {
      m_output( blockIdx.x * m_graph.block_size + threadIdx.x ) = y ;
    }
  }
};

//----------------------------------------------------------------------------

}

#endif /* #ifndef KOKKOSARRAY_BBCRS_HPP */

