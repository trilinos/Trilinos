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

#ifndef STOKHOS_HOST_BLOCKCRSMATRIX_HPP
#define STOKHOS_HOST_BLOCKCRSMATRIX_HPP

#include "KokkosArray_Host.hpp"
#include "KokkosArray_ParallelFor.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"

namespace Stokhos {

template< class BlockSpec , typename MatrixValue , typename VectorValue >
class Multiply<
  BlockCrsMatrix< BlockSpec , MatrixValue , KokkosArray::Host > ,
  KokkosArray::View< VectorValue** , KokkosArray::LayoutLeft , KokkosArray::Host > ,
  KokkosArray::View< VectorValue** , KokkosArray::LayoutLeft , KokkosArray::Host > ,
  DefaultSparseMatOps >
{
public:

  typedef KokkosArray::Host                                      device_type ;
  typedef device_type::size_type                    size_type ;
  typedef KokkosArray::View< VectorValue** , KokkosArray::LayoutLeft , KokkosArray::Host > block_vector_type ;
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

      Multiply< BlockSpec , void , void , DefaultSparseMatOps >::apply( m_A.block , a , x , y );
    }
  }

  static void apply( const matrix_type & A ,
                     const block_vector_type & x ,
                     const block_vector_type & y )
  {
    const size_t row_count = A.graph.row_map.dimension(0) - 1 ;
    KokkosArray::parallel_for( row_count , Multiply(A,x,y) );
  }
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_HOST_BLOCKCRSMATRIX_HPP */

