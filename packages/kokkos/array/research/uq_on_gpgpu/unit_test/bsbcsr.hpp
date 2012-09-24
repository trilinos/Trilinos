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

#ifndef KOKKOSARRAY_BSBCSR_HPP
#define KOKKOSARRAY_BSBCSR_HPP

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief  Dense symmetric matrix diagonal storage.
 *
 *  A dense symmetric NxN matrix is stored in
 *  a dense NxD matrix where D = number of diagonals stored.
 *
 *  If N is even then only the first half of the last diagonal is used.
 */
class DenseSymmetricMatrixDiagonalStorageMap {
public:

  const int m_length ;     // Actual length
  const int m_stride ;     // Stride for columns 
  const int m_diag_count ; // Number of diagonals to store
  const int m_size ;       // Size of storage

  // Offset in the storage
  __host__ __device__
  unsigned offset( const int row , const int column ) const
  {
    const int diag  = column - row ;

    int offset = 0 ;

    if ( ( 0 <= diag && diag < m_diag_count ) || ( diag <= - m_diag_count ) ) {
      offset = row + m_stride * ( ( m_length + diag ) % m_length );
    }
    else {
      offset = column + m_stride * ( ( m_length - diag ) % m_length );
    }

    return offset ;
  }

  // Storage maps to full matrix as:
  //   ( row , diagonal ) -> ( row , column( row , diagonal ) )

  __host__ __device__
  unsigned column( const int row , const int diagonal ) const
  {
    const int col = row + diagonal ;
    return m_length < col ? col - m_length : col ;
  }

  template< typename MatrixScalar , typename VectorScalar >
  void multiply( const MatrixScalar * const AD ,
                 const VectorScalar * const x ,
                       VectorScalar * const y )
  {
    const int even      = m_length & 01 ? 0 : 1 ;
    const int diag_full = m_diag_count - even ;

    // Main diagonal:
    for ( int row = 0 ; row < m_length ; ++row ) {
      y[row] = AD[row] * x[row] ;
    }
    // Full diagonals:
    const MatrixScalar * ad = AD ;
    for ( int diag = 1 ; diag < diag_full ; ++diag ) {
      int kx  = diag ;
      int kxr = diag + m_length - 2 * diag ;
      ad += m_stride ;

      for ( int row = 0 ; row < m_length ; ++row ) {
        y[row] += ad[ row ] * x[kx] ;
        y[row] += ad[ kxr ] * x[kxr] ;
        if ( m_length == ++kx )  kx  = 0 ;
        if ( m_length == ++kxr ) kxr = 0 ;
      }
    }
    // Final half diagonal:
    if ( even ) {
      const int half = m_length / 2 ;
      int kx  = diag_full ;

      ad += m_stride ;

      for ( int row = 0 ; row < half ; ++row ) {
        y[row] += ad[ row ] * x[kx] ;
        if ( m_length == ++kx ) kx = 0 ;
      }

      for ( int row = half ; row < m_length ; ++row ) {
        y[row] += ad[ row - half ] * x[kx] ;
        if ( m_length == ++kx ) kx = 0 ;
      }
    }
  }

  __host__ __device__
  DenseSymmetricMatrixDiagonalStorageMap( unsigned length , unsigned stride )
  : m_length( length )
  , m_stride( stride )
  , m_diag_count( 1 + length / 2 )
  , m_size( stride * m_diag_count )
  {}

};

//----------------------------------------------------------------------------
/** \brief  Big Symmetric Block CSR Graph
 *
 *  CSR Matrix storage with big NxN blocks instead of scalars.
 *  
 *  block_system_size :
 *    Number of row blocks, equal to number of column blocks
 *  
 *  block_length :
 *    Each dense block is block_length X diag_count
 *  
 *  block_column_offset[ 0 .. system_size ] :
 *    Offsets into the block_column_index array for each row-block
 *
 *  block_column_index[ 0 .. block_column_offset[ system_size ] - 1 ] :
 *    Index into the matrix coefficent array for the beginning
 *    of the block.
 */
template< class Device >
class BigSymmetricBlockCSRGraph {
public:
  typedef Device                     device_type ;
  typedef typename Device::size_type index_type ;

  typedef View< index_type[] , device_type > vector_type ;

  index_type  block_system_size ;
  index_type  block_length ;
  index_type  diag_stride ;
  index_type  diag_count ;
  vector_type block_column_offset ;
  vector_type block_column_index ;
};

//----------------------------------------------------------------------------

template< class Device , typename MatrixScalar , typename VectorScalar >
class BigSymmetricBlockCSRMultiply ;

template< class Device , typename MatrixScalar , typename VectorScalar >
void multiply( const BigSymmetricBlockCSRGraph<Device> & A_graph ,
               const View<MatrixScalar[],Device> & A_coeff ,
               const View<VectorScalar[],Device> & input ,
               const View<VectorScalar[],Device> & output )
{
  BigSymmetricBlockCSRMultiply< Device , MatrixScalar , VectorScalar >
    ( A_graph , A_coeff , input , output );
}

//----------------------------------------------------------------------------

template< typename MatrixScalar , typename VectorScalar >
class BigSymmetricBlockCSRMultiply< Cuda , MatrixScalar , VectorScalar > {
public:
  typedef Cuda::size_type                        index_type ;
  typedef BigSymmetricBlockCSRGraph< Cuda >      graph_type ;
  typedef View< MatrixScalar[] , Cuda > matrix_type ;
  typedef View< VectorScalar[] , Cuda > vector_type ;

  const graph_type  m_graph ;
  const matrix_type m_matrix ;
  const vector_type m_input ;
  const vector_type m_output ;
  const index_type  m_block_stride ;
  const index_type  m_diag_full ; // Number of full diagionals
  const index_type  m_shared_offset ;

  BigSymmetricBlockCSRMultiply( const graph_type  & graph ,
                                const matrix_type & matrix ,
                                const vector_type & input ,
                                const vector_type & output )
  : m_graph( graph )
  , m_matrix( matrix )
  , m_input(  input )
  , m_output( output )
  , m_block_stride( graph.diag_stride * graph.diag_count )
  , m_diag_full( graph.diag_count - ( graph.block_length & 01 ? 0 : 1 ) )
  , m_shared_offset( graph.block_length +
      ( graph.block_length % Impl::CudaTraits::WarpSize ?
      ( graph.block_length - 
        ( graph.block_length % Impl::CudaTraits::WarpSize ) ) : 0 ) )
  {
    const index_type thread_max = Cuda::maximum_warp_count() * Impl::CudaTraits::WarpSize ;

    if ( thread_max < graph.block_length ) {
      throw std::runtime_error( std::string("block is too big") );
    }

/*
    if ( Cuda::maximum_grid_count() < m_graph.block_system_size ) {
      throw std::runtime_error( std::string("too many blocks") );
    }
*/

    const index_type shmem_size = 2 * m_shared_offset * sizeof(VectorScalar);

    Impl::cuda_parallel_launch_local_memory<<< graph.block_system_size , graph.block_length , shmem_size >>>( *this );
  }

  //--------------------------------------------------------------------------
  // Shared memory usage:
  //   [ vector ]
  //   [ matrix column ]

  __device__
  void multiply( VectorScalar & y ,
                 const index_type vector_offset ,
                 const index_type matrix_offset ) const
  {
    extern __shared__ VectorScalar shared_input[] ;

    // Row of matrix == threadIdx.x

    // Load the input vector into share memory
    // using a coellesced memory read
    shared_input[ threadIdx.x ] = m_input( vector_offset + threadIdx.x );

    __syncthreads();

    int diag_offset = matrix_offset ;
    int kx = threadIdx.x ;

    // Multiply the main diagonal with the input vector

    y += m_matrix( threadIdx.x + diag_offset ) * shared_input[ kx ] ;

    // Multiply all full diagonals with the input vector

    for ( int diag = 1 ; diag < m_diag_full ; ++diag ) {

      diag_offset += m_graph.diag_stride ;

      // Load the diagonal into shared memory via coellesced memory read.
      // Full diagonals are accessed twice.

      shared_input[ threadIdx.x + m_shared_offset ] =
          m_matrix( diag_offset + threadIdx.x );

      __syncthreads();

      if ( m_graph.block_length == ++kx ) kx = 0 ;

      int kxr = kx + m_graph.block_length - ( diag << 1 ); // * 2
      if ( m_graph.block_length <= kxr ) kxr -= m_graph.block_length ;

      y += shared_input[threadIdx.x + m_shared_offset] * shared_input[kx];
      y += shared_input[kxr         + m_shared_offset] * shared_input[kxr];
    }

    if ( m_diag_full < m_graph.diag_count ) {
      const int half = m_graph.block_length >> 1 ;

      if ( m_graph.block_length == ++kx ) kx = 0 ;

      diag_offset += m_graph.diag_stride ;

      int i = threadIdx.x ; 

      // Load the half diagonal into shared memory via coellesced memory read
      // Each entry is accessed by two threads.

      if ( threadIdx.x < half ) {
        shared_input[ threadIdx.x + m_graph.block_length ] =
          m_matrix( diag_offset + threadIdx.x );
      }
      else {
        i -= half ;
      }
      __syncthreads();

      y += shared_input[ i + m_graph.block_length ] * shared_input[kx] ;
    }
  }

  __device__
  void execute_on_device() const
  {
    // Require: m_graph.block_length == blockDim.x

    // Global memory read broadcast to threads in the block:
    const index_type block_begin = m_graph.block_column_offset(blockIdx.x);
    const index_type block_end   = m_graph.block_column_offset(blockIdx.x + 1);

    index_type matrix_offset = block_begin * m_block_stride ;

    VectorScalar y = 0 ;

    for ( index_type i_block = block_begin ; i_block < block_end ;
          ++i_block , matrix_offset += m_block_stride ) {

      const index_type vector_offset =
        m_graph.block_column_index( i_block ) * m_graph.block_length ;

      multiply( y , vector_offset , matrix_offset );
    }

    m_output( threadIdx.x + blockIdx.x * m_graph.block_length ) = y ;
  }
};

//----------------------------------------------------------------------------

}

#endif /* #ifndef KOKKOSARRAY_BSBCSR_HPP */

