
#include <utility>
#include <cmath>
#include <iostream>

#include <impl/KokkosArray_Timer.hpp>

namespace unit_test {

//----------------------------------------------------------------------------

template< typename VectorScalar , typename MatrixScalar , typename TensorScalar , class Device >
std::pair<size_t,double>
test_product_tensor_legendre(
  const std::vector<int> & arg_var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef KokkosArray::View< VectorScalar** ,
                             KokkosArray::LayoutLeft ,
                             Device > vector_type ;

  typedef KokkosArray::CrsProductTensorLegendre< TensorScalar , Device >  tensor_type ;

  typedef KokkosArray::BlockCrsMatrix< tensor_type , MatrixScalar , Device > matrix_type ;

  typedef typename matrix_type::graph_type graph_type ;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  std::vector<unsigned> var_degree( arg_var_degree.size() , 0u );

  unsigned maximum_degree = 2 ; // Minimum degree to capture pair-wise covariance
  for ( unsigned i = 0 ; i < arg_var_degree.size() ; ++i ) {
    var_degree[i] = arg_var_degree[i];
    maximum_degree = std::max( maximum_degree , var_degree[i] );
  }

  matrix_type matrix ;

  matrix.block = tensor_type( var_degree , maximum_degree );

  matrix.graph = KokkosArray::create_crsarray<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();
  const size_t inner_matrix_size = matrix.block.dimension();

  matrix.values = vector_type( "matrix" , inner_matrix_size , graph_length );

  vector_type x = vector_type( "x" , inner_length , outer_length );
  vector_type y = vector_type( "y" , inner_length , outer_length );

  typename vector_type::HostMirror hM = KokkosArray::create_mirror( matrix.values );
  
  for ( size_t i = 0 ; i < graph_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1 + i ;
    }
  }
  
  KokkosArray::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:
  
  typename vector_type::HostMirror hx = KokkosArray::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1 + j + 10 * i ;
    }
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  const KokkosArray::Impl::Multiply< matrix_type , vector_type , vector_type > op( matrix , x , y );

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    op.run();
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );

  //------------------------------

  if ( print_flag ) {
    typename vector_type::HostMirror hy = KokkosArray::create_mirror( y );

    KokkosArray::deep_copy( hy , y );

    std::cout << std::endl << "test_product_tensor_matrix" << std::endl ;
    for ( size_t i = 0 ; i < outer_length ; ++i ) {
      std::cout << "hy(:," << i << ") =" ;
      for ( size_t j = 0 ; j < inner_length ; ++j ) {
        std::cout << " " << hy(j,i);
      }
      std::cout << std::endl ;
    }
  }

  return std::pair<size_t,double>( outer_length * inner_length , seconds_per_iter );
}

//----------------------------------------------------------------------------

} // namespace unit_test



