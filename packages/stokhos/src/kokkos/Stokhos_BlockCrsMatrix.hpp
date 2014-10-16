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

#ifndef STOKHOS_BLOCKCRSMATRIX_HPP
#define STOKHOS_BLOCKCRSMATRIX_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_CrsArray.hpp"

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

  typedef Device device_type;
  typedef typename device_type::size_type size_type;
  typedef ValueType value_type;
  typedef BlockSpec block_spec;
  typedef Kokkos::CrsArray< size_type , device_type > graph_type;
  typedef Kokkos::View< value_type**, Kokkos::LayoutLeft, device_type > block_vector_type ;

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

  typedef Device device_type ;
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
    const size_t row_count = A.graph.row_map.dimension_0() - 1;
    Kokkos::parallel_for( row_count , Multiply(A,x,y) );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_BLOCKCRSMATRIX_HPP */
