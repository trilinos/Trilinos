/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/



#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <Kokkos_ProductTensor.hpp>
#include <Kokkos_BlockCrsMatrix.hpp>

namespace unit_test_tensor {

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
template< typename IntType >
inline
IntType generate_matrix_value( const IntType inner_row ,
                              const IntType outer_row ,
                              const IntType outer_column )
{
  return 1 + inner_row + 10 * outer_row + 20 * outer_column ;
}

template< typename Scalar , class Device >
inline
void generate_tensor( const size_t M , std::map< Kokkos::ProductTensorIndex<3,Device> , Scalar > & tensor )
{
  typedef Kokkos::ProductTensorIndex<3,Device> index_type ;

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
  Kokkos::BlockCrsMatrix<Kokkos::SparseProductTensor<3,ScalarType,Device>,ScalarType,Device> & matrix )
{
  typedef Kokkos::MultiVector< ScalarType , Device > values_type ;
  typedef Kokkos::CrsMap< Device >                   graph_type ;
  typedef Kokkos::SparseProductTensor<3,ScalarType,Device> tensor_type ;
  typedef Kokkos::ProductTensorIndex<3,Device>             index_type ;

  typedef typename values_type::HostMirror host_values_type ;
  typedef typename graph_type ::HostMirror host_graph_type ;

  std::vector< std::vector<size_t> > graph( N * N * N );

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

  std::map< index_type , ScalarType > tensor_input ;

  generate_tensor( M , tensor_input );

  matrix.block  = Kokkos::create_sparse_product_tensor<3,ScalarType,Device>( tensor_input );
  matrix.graph  = Kokkos::create_labeled_crsmap<Device>( std::string("test crs graph") , graph );
  matrix.values = Kokkos::create_multivector<ScalarType,Device>( matrix.block.size() , total );

  host_graph_type  h_graph  = Kokkos::create_mirror( matrix.graph );
  host_values_type h_values = Kokkos::create_mirror( matrix.values );

  for ( size_t outer_row = 0 ; outer_row < N*N*N ; ++outer_row ) {
    const size_t outer_entry_begin = h_graph.row_range_begin( outer_row );
    const size_t outer_entry_end   = h_graph.row_range_end( outer_row );

    for ( size_t outer_entry = outer_entry_begin ;
                 outer_entry < outer_entry_end ; ++outer_entry ) {
      const size_t outer_column = h_graph.column(outer_entry);

      for ( size_t inner_entry = 0 ; inner_entry < M ; ++inner_entry ) {
        h_values(inner_entry,outer_entry) =
          generate_matrix_value( inner_entry , outer_row , outer_column );
      }
    }
  }

  Kokkos::deep_copy( matrix.values , h_values );
}


template< class Device , typename IntType >
void test_tensor_crs_matrix( const size_t M , const size_t N , const bool print = false )
{
  const size_t length = N * N * N ;

  typedef IntType value_type ; // to avoid comparison round-off differences

  typedef Kokkos::CrsMap< Device >        graph_type ;
  typedef typename graph_type::HostMirror host_graph_type ;
  typedef Kokkos::SparseProductTensor< 3 , IntType , Device > block_spec ;
  typedef Kokkos::ProductTensorIndex<3,Device>             index_type ;

  Kokkos::BlockCrsMatrix< block_spec , value_type , Device > matrix ;

  generate_matrix( M , N , matrix );

  Kokkos::MultiVector<value_type,Device> x = Kokkos::create_multivector<value_type,Device>( M , length );
  Kokkos::MultiVector<value_type,Device> y = Kokkos::create_multivector<value_type,Device>( M , length );

  typename Kokkos::MultiVector<value_type,Device>::HostMirror hx = Kokkos::create_mirror( x );
  typename Kokkos::MultiVector<value_type,Device>::HostMirror hy = Kokkos::create_mirror( y );

  for ( size_t i = 0 ; i < length ; ++i ) {
    for ( size_t j = 0 ; j < M ; ++j ) {
      hx(j,i) = 1 + j + 10 * i ;
    }
  }

  Kokkos::deep_copy( x , hx );

  Kokkos::multiply( matrix , x , y );

  Kokkos::deep_copy( hy , y );

  // Check answer:

  host_graph_type h_graph  = Kokkos::create_mirror( matrix.graph );

  std::map< index_type , IntType > tensor_input ;

  generate_tensor( M , tensor_input );

  for ( size_t outer_row = 0 ; outer_row < length ; ++outer_row ) {
    const size_t outer_entry_begin = h_graph.row_range_begin( outer_row );
    const size_t outer_entry_end   = h_graph.row_range_end( outer_row );

    for ( size_t inner_row = 0 ; inner_row < M ; ++inner_row ) {

      value_type value = 0 ;

      for ( size_t outer_entry = outer_entry_begin ;
                   outer_entry < outer_entry_end ; ++outer_entry ) {

        const size_t outer_column = h_graph.column( outer_entry );

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

