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



#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <KokkosArray_ProductTensor.hpp>
#include <KokkosArray_CrsProductTensorLegendre.hpp>
#include <KokkosArray_BlockCrsMatrix.hpp>

#include <TestGenerateGraph.hpp>

namespace unit_test_tensor {

// Some arbitrary values for testing.
template< typename IntType >
inline
IntType generate_matrix_value( const IntType inner_row ,
                              const IntType outer_row ,
                              const IntType outer_column )
{
  return 1 + inner_row + 10 * outer_row + 20 * outer_column ;
}

template< typename Scalar >
inline
void generate_tensor( const size_t M ,
                      std::map< KokkosArray::ProductTensorIndex<3> , Scalar > & tensor )
{
  typedef KokkosArray::ProductTensorIndex<3> index_type ;

  for ( size_t i = 0 ; i < M ; ++i ) {
    if ( 0 < i ) {
      tensor[ index_type(i,0,0) ] = i + 1 ;
      tensor[ index_type(0,i,0) ] = i + 2 ;
      tensor[ index_type(0,0,i) ] = i + 3 ;
    }
    tensor[ index_type(i,i,i) ] = i + 4 ;
  }
}

template< typename ScalarType , class Device >
void generate_matrix(
  const size_t M ,
  const size_t N ,
  KokkosArray::BlockCrsMatrix<KokkosArray::SparseProductTensor<3,ScalarType,Device>,ScalarType,Device> & matrix )
{
  typedef KokkosArray::BlockCrsMatrix<KokkosArray::SparseProductTensor<3,ScalarType,Device>,ScalarType,Device> matrix_type ;

  typedef typename matrix_type::block_vector_type    block_vector_type ;
  typedef typename matrix_type::graph_type           graph_type ;
  typedef typename matrix_type::block_spec           tensor_type ;
  typedef KokkosArray::ProductTensorIndex<3>         index_type ;

  typedef KokkosArray::Impl::Multiply< tensor_type > tensor_multiply ;

  std::vector< std::vector<size_t> > graph ;

  std::map< index_type , ScalarType > tensor_input ;

  generate_tensor( M , tensor_input );

  matrix.block = KokkosArray::create_product_tensor<tensor_type>( tensor_input );

  const size_t total      = unit_test::generate_fem_graph( N , graph );
  const size_t block_size = tensor_multiply::matrix_size( matrix.block );

  matrix.graph  = KokkosArray::create_crsarray<graph_type>( std::string("test crs graph") , graph );
  matrix.values = block_vector_type( "matrix_values" , block_size , total );

  typename graph_type::HostMirror h_graph =
    KokkosArray::create_mirror( matrix.graph );

  typename block_vector_type::HostMirror h_values =
    KokkosArray::create_mirror( matrix.values );

  for ( size_t outer_row = 0 ; outer_row < N*N*N ; ++outer_row ) {
    const size_t outer_entry_begin = h_graph.row_map[outer_row];
    const size_t outer_entry_end   = h_graph.row_map[outer_row+1];

    for ( size_t outer_entry = outer_entry_begin ;
                 outer_entry < outer_entry_end ; ++outer_entry ) {
      const size_t outer_column = h_graph.entries(outer_entry);

      for ( size_t inner_entry = 0 ; inner_entry < M ; ++inner_entry ) {
        h_values(inner_entry,outer_entry) =
          generate_matrix_value( inner_entry , outer_row , outer_column );
      }
    }
  }

  KokkosArray::deep_copy( matrix.values , h_values );
}


template< class Device , typename IntType >
void test_tensor_crs_matrix( const size_t M , const size_t N , const bool print = false )
{
  const size_t length = N * N * N ;

  typedef IntType value_type ; // to avoid comparison round-off differences
  typedef KokkosArray::View<value_type**,KokkosArray::LayoutLeft,Device>
    block_vector_type ;

  typedef KokkosArray::SparseProductTensor< 3 , IntType , Device > block_spec ;
  typedef KokkosArray::BlockCrsMatrix< block_spec , value_type , Device > matrix_type ;

  typedef typename matrix_type::graph_type    graph_type ;
  typedef KokkosArray::ProductTensorIndex<3>  index_type ;

  matrix_type matrix ;

  generate_matrix( M , N , matrix );

  block_vector_type x = block_vector_type( "x" , M , length );
  block_vector_type y = block_vector_type( "y" , M , length );

  typename block_vector_type::HostMirror hx = KokkosArray::create_mirror( x );
  typename block_vector_type::HostMirror hy = KokkosArray::create_mirror( y );

  for ( size_t i = 0 ; i < length ; ++i ) {
    for ( size_t j = 0 ; j < M ; ++j ) {
      hx(j,i) = 1 + j + 10 * i ;
    }
  }

  KokkosArray::deep_copy( x , hx );

  KokkosArray::multiply( matrix , x , y );

  KokkosArray::deep_copy( hy , y );

  // Check answer:

  typename graph_type::HostMirror h_graph =
    KokkosArray::create_mirror( matrix.graph );

  std::map< index_type , IntType > tensor_input ;

  generate_tensor( M , tensor_input );

  for ( size_t outer_row = 0 ; outer_row < length ; ++outer_row ) {
    const size_t outer_entry_begin = h_graph.row_map[outer_row];
    const size_t outer_entry_end   = h_graph.row_map[outer_row+1];

    for ( size_t inner_row = 0 ; inner_row < M ; ++inner_row ) {

      value_type value = 0 ;

      for ( size_t outer_entry = outer_entry_begin ;
                   outer_entry < outer_entry_end ; ++outer_entry ) {

        const size_t outer_column = h_graph.entries( outer_entry );

        for ( typename std::map< index_type , IntType >::iterator
              iter =  tensor_input.begin() ;
              iter != tensor_input.end() ; ++iter ) {

          const index_type index = (*iter).first ;
          const IntType    coeff = (*iter).second ;

          size_t i , j ;

          if ( inner_row == index.coord(2) ) {
            i = index.coord(0);
            j = index.coord(1);
          }
          else if ( inner_row == index.coord(1) ) {
            i = index.coord(0);
            j = index.coord(2);
          }
          else if ( inner_row == index.coord(0) ) {
            i = index.coord(1);
            j = index.coord(2);
          }
          else {
            continue ;
          }

          const IntType ai = generate_matrix_value( i , outer_row , outer_column );
          const IntType aj = generate_matrix_value( j , outer_row , outer_column );

          value += coeff * ( i == j ? ( ai * hx( i , outer_column ) )
                                    : ( ai * hx( j , outer_column ) +
                                        aj * hx( i , outer_column ) ) );
        }
      }

      if ( value != hy(inner_row,outer_row) ) {
        std::ostringstream msg ;
        msg << "correctness test failed "
            << value << " != " << hy(inner_row,outer_row)
            << " = hy( " << inner_row << "," << outer_row << ")" ;
        throw std::runtime_error(msg.str());
      }
      else if ( print ) {
        std::cout << " = hy( " << inner_row << "," << outer_row << ") = "
                  << hy( inner_row , outer_row ) << std::endl ;
      }
    }
  }
}

}

