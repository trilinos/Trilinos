// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_BLOCKCRSMATRIX_HPP
#define STOKHOS_BLOCKCRSMATRIX_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_StaticCrsGraph.hpp"

#include "Stokhos_Multiply.hpp"

namespace Stokhos {

/** \brief  CRS matrix of dense blocks.
 *
 *  Matrix coefficients are stored by block and then by Crs entry.
 *    m_values( block.size() , m_graph.entry_count() )
 *
 *  Vectors are conformally stored as
 *    View( block.dimension() , m_graph.row_map.length() )
 */
template <typename BlockSpec, typename ValueType, class Device>
class BlockCrsMatrix {
public:

  typedef Device execution_space;
  typedef typename execution_space::size_type size_type;
  typedef ValueType value_type;
  typedef BlockSpec block_spec;
  typedef Kokkos::StaticCrsGraph< size_type , execution_space > graph_type;
  typedef Kokkos::View< value_type**, Kokkos::LayoutLeft, execution_space > block_vector_type ;

  block_vector_type  values;
  graph_type         graph;
  block_spec         block;
};

template <typename BlockSpec,
          typename MatrixValue,
          typename VectorValue,
          typename Device>
class Multiply< BlockCrsMatrix< BlockSpec, MatrixValue, Device >,
                Kokkos::View< VectorValue**, Kokkos::LayoutLeft, Device >,
                Kokkos::View< VectorValue**, Kokkos::LayoutLeft, Device > >
{
public:

  typedef Device execution_space ;
  typedef typename BlockSpec::size_type size_type ;
  typedef Kokkos::View< VectorValue**, Kokkos::LayoutLeft, Device > block_vector_type ;
  typedef BlockCrsMatrix< BlockSpec, MatrixValue, Device >  matrix_type ;

  const matrix_type m_A;
  const block_vector_type m_x;
  block_vector_type m_y;

  Multiply( const matrix_type& A,
            const block_vector_type& x,
            block_vector_type& y )
  : m_A( A )
  , m_x( x )
  , m_y( y )
  {}

  //--------------------------------------------------------------------------
  //  A( storage_size( m_A.block.size() ) , m_A.graph.row_map.size() );
  //  x( m_A.block.dimension() , m_A.graph.row_map.first_count() );
  //  y( m_A.block.dimension() , m_A.graph.row_map.first_count() );
  //

  KOKKOS_INLINE_FUNCTION
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

      BlockMultiply< BlockSpec >::apply( m_A.block , a , x , y );
    }

  }

  static void apply( const matrix_type& A,
                     const block_vector_type& x,
                     block_vector_type& y )
  {
    const size_t row_count = A.graph.row_map.extent(0) - 1;
    Kokkos::parallel_for( row_count , Multiply(A,x,y) );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_BLOCKCRSMATRIX_HPP */
