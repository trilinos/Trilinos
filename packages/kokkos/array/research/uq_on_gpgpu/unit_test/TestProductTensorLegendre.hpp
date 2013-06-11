
#include <utility>
#include <cmath>
#include <iostream>

#include <TestGenerateTensor.hpp>
#include <TestGenerateGraph.hpp>
#include <impl/KokkosArray_Timer.hpp>

namespace unit_test {

//----------------------------------------------------------------------------

template< class TensorType , typename VectorScalar , typename MatrixScalar >
std::vector<double>
test_product_tensor_legendre(
  const std::vector<int> & arg_var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool check )
{
  typedef TensorType                         tensor_type ;
  typedef typename tensor_type::device_type  device_type ;

  typedef KokkosArray::View< VectorScalar** ,
                             KokkosArray::LayoutLeft ,
                             device_type > vector_type ;

  typedef KokkosArray::BlockCrsMatrix< tensor_type , MatrixScalar , device_type > matrix_type ;

  typedef typename matrix_type::graph_type graph_type ;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;
  const size_t fem_graph_length = unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  const std::vector<unsigned> var_degree( arg_var_degree.begin() , arg_var_degree.end() );

  const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
    tensor( var_degree );

  const size_t stoch_length = tensor.bases_count();

  std::vector< std::vector< size_t > > stoch_graph( stoch_length );

  for ( size_t i = 0 ; i < stoch_length ; ++i ) {
    for ( size_t j = 0 ; j < stoch_length ; ++j ) {
      if ( KokkosArray::matrix_nonzero(tensor,i,j) ) {
        stoch_graph[i].push_back(j);
      }
    }
  }

  //------------------------------
  // Generate input multivector:
  
  vector_type x = vector_type( "x" , stoch_length , fem_length );
  vector_type y = vector_type( "y" , stoch_length , fem_length );

  typename vector_type::HostMirror hx        = KokkosArray::create_mirror( x );
  typename vector_type::HostMirror hy_result = KokkosArray::create_mirror( y );

  for ( size_t iColFEM = 0 ;   iColFEM < fem_length ;   ++iColFEM ) {
  for ( size_t iColStoch = 0 ; iColStoch < stoch_length ; ++iColStoch ) {
    hx(iColStoch,iColFEM) =
      generate_vector_coefficient( fem_length , stoch_length ,
                                   iColFEM , iColStoch );
  }}

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.block = tensor_type( var_degree );

  matrix.graph = KokkosArray::create_crsarray<graph_type>( std::string("test crs graph") , fem_graph );

  if ( stoch_length != matrix.block.dimension() ) {
    throw std::runtime_error("test_crs_product_tensor_legendre matrix sizing error");
  }

  matrix.values = vector_type( "matrix" , stoch_length , fem_graph_length );

  typename vector_type::HostMirror hM = KokkosArray::create_mirror( matrix.values );

  for ( size_t iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < fem_length ; ++iRowFEM ) {
    for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
      const size_t iColFEM = fem_graph[iRowFEM][iRowEntryFEM] ;

      for ( size_t k = 0 ; k < stoch_length ; ++k ) {
        hM(k,iEntryFEM) = generate_matrix_coefficient( fem_length , stoch_length , iRowFEM , iColFEM , k );
      }
    }
  }

  KokkosArray::deep_copy( matrix.values , hM );

  //------------------------------

  if (check) {
    for ( size_t iRowStoch = 0 ; iRowStoch < stoch_length ; ++iRowStoch ) {
      for ( size_t iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < fem_length ; ++iRowFEM ) {

	double y = 0 ;

	for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < fem_graph[ iRowFEM ].size() ; ++iRowEntryFEM , ++iEntryFEM ) {

	  const size_t iColFEM = fem_graph[iRowFEM][iRowEntryFEM] ;

	  for ( size_t iRowEntryStoch = 0 ; iRowEntryStoch < stoch_graph[iRowStoch].size() ; ++iRowEntryStoch ) {

	    const size_t iColStoch = stoch_graph[iRowStoch][iRowEntryStoch];

	    double value = 0 ;
	    for ( unsigned k = 0 ; k < stoch_length ; ++k ) {

	      const double A_fem_k = generate_matrix_coefficient( fem_length , stoch_length , iRowFEM , iColFEM , k );

	      if ( 1.0e-6 < std::abs( hM(k,iEntryFEM) - A_fem_k ) ) {
		std::cout << "test_crs_product_tensor_legendre error: Matrix entry"
			  << "  A(" << k << ",(" << iRowFEM << "," << iColFEM << ")) = " << hM(k,iEntryFEM) 
			  << " , error = " << hM(k,iEntryFEM) - A_fem_k
			  << std::endl ;
	      }
	      
	      value += tensor(iRowStoch,iColStoch,k) * A_fem_k ;
	    }
	    
	    y += value * hx( iColStoch , iColFEM );
	  }
	}

	hy_result( iRowStoch , iRowFEM ) = y ;
      }
    }
  }

  //------------------------------

  const KokkosArray::Impl::Multiply< matrix_type , vector_type , vector_type > op( matrix , x , y );

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    op.run();
  }
  device_type::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = matrix.block.multiply_add_flops();
  const double flops = 1.0e-9*fem_graph_length*flops_per_block / seconds_per_iter;

  //------------------------------
  // Verify result

  if (check)
  {
    const double tol = KokkosArray::Impl::is_same<double,VectorScalar>::value ? 1.0e-13 : 1.0e-5 ;
    const size_t error_max = 10 ;

    KokkosArray::deep_copy( hx , y );

    size_t error_count = 0 ;

    for ( size_t iRowFEM = 0 ; iRowFEM < fem_length ; ++iRowFEM ) {
      for ( size_t iRowStoch = 0 ; iRowStoch < stoch_length ; ++iRowStoch ) {
        const double mag   = std::abs( hy_result(iRowStoch,iRowFEM) );
        const double error = std::abs( hx(iRowStoch,iRowFEM) - hy_result(iRowStoch,iRowFEM) );
        if ( tol < error && tol < error / mag ) {
          if ( error_count < error_max ) {
            std::cout << "test_product_tensor_legendre error:"
                      << " y(" << iRowStoch << "," << iRowFEM << ") = " << hx(iRowStoch,iRowFEM)
                      << " , error = " << ( hx(iRowStoch,iRowFEM) - hy_result(iRowStoch,iRowFEM) )
                      << std::endl ;
          }
          ++error_count ;
        }
      }
    }
    if ( error_count ) {
      std::cout << "test_crs_product_tensor_legendre error_count = " << error_count << std::endl ;
    }
  }

  //------------------------------

  std::vector<double> perf(3) ;

  perf[0] = fem_length * stoch_length ;
  perf[1] = seconds_per_iter ;
  perf[2] = flops ;

  return perf ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace unit_test



