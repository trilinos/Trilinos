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

#include <stdexcept>
#include <sstream>
#include <vector>
#include <KokkosArray_SymmetricDiagonalSpec.hpp>
#include <KokkosArray_BlockCrsMatrix.hpp>

namespace unit_test {

template< typename IntType >
inline
IntType map_coord( const IntType & N ,
                   const IntType & i ,
                   const IntType & j ,
                   const IntType & k )
{
  return k + N * ( j + N * i );
}

// Some arbitrary values for testing.
inline
long generate_matrix_value( const long inner_row ,
                            const long inner_column ,
                            const long outer_row ,
                            const long outer_column )
{
  return std::min( inner_column , inner_row ) +
         2 * std::max( inner_column , inner_row ) +
         10 * outer_row + 20 * outer_column ;
}

template< class Device >
void generate_matrix(
  const size_t M ,
  const size_t N ,
  KokkosArray::BlockCrsMatrix<KokkosArray::SymmetricDiagonalSpec<Device>,long,Device> & matrix )
{
  typedef KokkosArray::BlockCrsMatrix<KokkosArray::SymmetricDiagonalSpec<Device>,long,Device> matrix_type ;
  typedef KokkosArray::View< long** , KokkosArray::LayoutLeft , Device >
     block_vector_type ;

  typedef typename matrix_type::graph_type  graph_type ;
  typedef typename matrix_type::block_spec  block_type ;

  const size_t outer_length = N * N * N ;

  std::vector< std::vector<size_t> > graph( outer_length );

  size_t total = 0 ;

  for ( int i = 0 ; i < (int) N ; ++i ) {
  for ( int j = 0 ; j < (int) N ; ++j ) {
  for ( int k = 0 ; k < (int) N ; ++k ) {
    const size_t row = map_coord((int)N,i,j,k);
    graph[row].reserve(27);

    for ( int ii = -1 ; ii < 2 ; ++ii ) {
    for ( int jj = -1 ; jj < 2 ; ++jj ) {
    for ( int kk = -1 ; kk < 2 ; ++kk ) {
      if ( 0 <= i + ii && i + ii < (int) N &&
           0 <= j + jj && j + jj < (int) N &&
           0 <= k + kk && k + kk < (int) N ) {
        size_t col = map_coord((int)N,i+ii,j+jj,k+kk);

        graph[row].push_back(col); 
      }
    }}}
    total += graph[row].size();
  }}}

  matrix.block  = block_type( M );

  const size_t block_size = KokkosArray::Impl::Multiply< block_type >::matrix_size( matrix.block );

  matrix.graph  = KokkosArray::create_crsarray<graph_type>( std::string("test crs graph") , graph );
  matrix.values = block_vector_type( "matrix" , block_size , total );

  typename graph_type::HostMirror h_graph =
    KokkosArray::create_mirror( matrix.graph );

  typename block_vector_type::HostMirror h_values =
     KokkosArray::create_mirror( matrix.values );

  for ( size_t outer_row = 0 ; outer_row < outer_length ; ++outer_row ) {
    const size_t outer_entry_begin = h_graph.row_map[outer_row];
    const size_t outer_entry_end   = h_graph.row_map[outer_row+1];

    for ( size_t outer_entry = outer_entry_begin ;
                 outer_entry < outer_entry_end ; ++outer_entry ) {
      const size_t outer_column = h_graph.entries(outer_entry);

      for ( size_t inner_row = 0 ; inner_row < M ; ++inner_row ) {
        for ( size_t inner_column = 0 ; inner_column <= inner_row ; ++inner_column ) {

          const size_t inner_entry = matrix.block.matrix_offset(inner_row,inner_column);

          h_values(inner_entry,outer_entry) =
            generate_matrix_value( inner_row , inner_column , outer_row , outer_column );
        }
      }
    }
  }

  KokkosArray::deep_copy( matrix.values , h_values );
}


template< class Device >
void test_block_crs_matrix( const size_t M , const size_t N )
{
  const size_t outer_length = N * N * N ;

  typedef long value_type ; // to avoid comparison round-off differences
  typedef KokkosArray::View<value_type**,KokkosArray::LayoutLeft,Device> block_vector_type ;

  typedef KokkosArray::SymmetricDiagonalSpec< Device > block_spec ;
  typedef KokkosArray::BlockCrsMatrix< block_spec , value_type , Device > matrix_type ;

  typedef typename matrix_type::graph_type  graph_type ;
  typedef typename graph_type::HostMirror   host_graph_type ;

  matrix_type matrix ;

  generate_matrix( M , N , matrix );

  block_vector_type x = block_vector_type( "x" , M , outer_length );
  block_vector_type y = block_vector_type( "y" , M , outer_length );

  typename block_vector_type::HostMirror hx = KokkosArray::create_mirror( x );
  typename block_vector_type::HostMirror hy = KokkosArray::create_mirror( y );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < M ; ++j ) {
      hx(j,i) = 1 + j + 10 * i ;
    }
  }

  KokkosArray::deep_copy( x , hx );

  KokkosArray::multiply( matrix , x , y );

  KokkosArray::deep_copy( hy , y );

  // Check answer:

  host_graph_type h_graph  = KokkosArray::create_mirror( matrix.graph );

  for ( size_t outer_row = 0 ; outer_row < outer_length ; ++outer_row ) {
    const size_t outer_entry_begin = h_graph.row_map[outer_row];
    const size_t outer_entry_end   = h_graph.row_map[outer_row+1];

    for ( size_t inner_row = 0 ; inner_row < M ; ++inner_row ) {

      value_type value = 0 ;

      for ( size_t outer_entry = outer_entry_begin ;
                   outer_entry < outer_entry_end ; ++outer_entry ) {

        const size_t outer_column = h_graph.entries( outer_entry );

        for ( size_t inner_column = 0 ; inner_column < M ; ++inner_column ) {

          value += hx( inner_column , outer_column ) *
                   generate_matrix_value( inner_row , inner_column ,
                                          outer_row , outer_column );
        }
      }

      if ( value != hy(inner_row,outer_row) ) {
        std::ostringstream msg ;
        msg << "correctness test failed "
            << value << " != " << hy(inner_row,outer_row)
            << " = hy( " << inner_row << "," << outer_row << ")" ;
        throw std::runtime_error(msg.str());
      }
    }
  }
}


}

