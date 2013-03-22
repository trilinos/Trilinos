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

#ifndef KOKKOSARRAY_HOST_BLOCKCRSMATRIX_HPP
#define KOKKOSARRAY_HOST_BLOCKCRSMATRIX_HPP

namespace KokkosArray {
namespace Impl {

template< class BlockSpec , typename MatrixValue , typename VectorValue >
class Multiply<
  BlockCrsMatrix< BlockSpec , MatrixValue , Host > ,
  KokkosArray::View< VectorValue** , LayoutLeft , Host > ,
  KokkosArray::View< VectorValue** , LayoutLeft , Host > >
{
public:
  typedef Host                                      device_type ;
  typedef device_type::size_type                    size_type ;
  typedef View< VectorValue** , LayoutLeft , Host > block_vector_type ;
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

  inline
  void operator()( const size_type iBlockRow ) const
  {
    // Prefer that y[ m_A.block.dimension() ] be scratch space
    // on the local thread, but cannot dynamically allocate
    VectorValue * const y = & m_y(0,iBlockRow);

    const size_type iEntryBegin = m_A.graph.row_map[ iBlockRow ];
    const size_type iEntryEnd   = m_A.graph.row_map[ iBlockRow + 1 ];

    // Leading dimension guaranteed contiguous for LayoutLeft
    for ( size_type j = 0 ; j < m_A.block.dimension() ; ++j ) { y[j] = 0 ; }

    for ( size_type iEntry = iEntryBegin ; iEntry < iEntryEnd ; ++iEntry ) {
      const VectorValue * const x = & m_x( 0 , m_A.graph.entries(iEntry) );
      const MatrixValue * const a = & m_A.values( 0 , iEntry );

      Multiply< BlockSpec , void , void >::apply( m_A.block , a , x , y );
    }
  }

  static void apply( const matrix_type & A ,
                     const block_vector_type & x ,
                     const block_vector_type & y )
  {
    const size_t row_count = A.graph.row_map.dimension_0() - 1 ;
    parallel_for( row_count , Multiply(A,x,y) );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_HOST_BLOCKCRSMATRIX_HPP */

