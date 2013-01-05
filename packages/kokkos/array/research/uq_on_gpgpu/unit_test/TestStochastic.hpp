
#include <utility>
#include <cmath>
#include <iostream>

#include <impl/KokkosArray_Timer.hpp>
#include <TestProductTensorLegendre.hpp>
#include <TestGenerateTensor.hpp>

#ifdef HAVE_KOKKOSARRAY_STOKHOS
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#endif

namespace unit_test {

template <typename Scalar> struct ScalarTolerances {};

template <> struct ScalarTolerances<float> {
  typedef float scalar_type;
  static scalar_type sparse_cijk_tol() { return 1e-5; }
};

template <> struct ScalarTolerances<double> {
  typedef double scalar_type;
  static scalar_type sparse_cijk_tol() { return 1e-12; }
};

template< typename ScalarType , class Device ,
          template< unsigned , typename , class > class TensorType >
std::vector<double>
test_product_tensor_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool check ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type** ,
                             KokkosArray::LayoutLeft ,
                             Device > block_vector_type ;

  typedef KokkosArray::NormalizedLegendrePolynomialBases<4> polynomial ;

  typedef KokkosArray::StochasticProductTensor< value_type , polynomial , Device , TensorType > tensor_type ;

  typedef KokkosArray::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  matrix_type matrix ;

  matrix.block = KokkosArray::create_product_tensor< tensor_type >( var_degree );
  matrix.graph = KokkosArray::create_crsarray<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();
  const size_t inner_matrix_size = matrix.block.dimension();

  matrix.values = block_vector_type( "matrix" , inner_matrix_size , graph_length );

  block_vector_type x = block_vector_type( "x" , inner_length , outer_length );
  block_vector_type y = block_vector_type( "y" , inner_length , outer_length );

  typename block_vector_type::HostMirror hM = KokkosArray::create_mirror( matrix.values );

  // Loop over spatial rows
  size_t outer_nz = 0;
  for ( size_t iOuterRow = 0 ; iOuterRow < outer_length ; ++iOuterRow ) {

    //Loop over non-zeros in a row
    const size_t iOuterNZ = graph[iOuterRow].size();
    for ( size_t iOuterEntry = 0 ; iOuterEntry < iOuterNZ ; ++iOuterEntry ) {
      //const size_t iOuterCol = graph[iOuterRow][iOuterEntry];
  
      // Build PC expansion for this nonzero
      for ( size_t j = 0 ; j < inner_length ; ++j ) {
	hM(j,outer_nz) = 1 + j + 0.10 * outer_nz ;
      }
      outer_nz++;
    }
  }
  
  // for ( size_t i = 0 ; i < graph_length ; ++i ) {
  //   for ( size_t j = 0 ; j < inner_length ; ++j ) {
  //     hM(j,i) = 1 + j + 10 * i ;
  //   }
  // }
  
  KokkosArray::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:
  
  typename block_vector_type::HostMirror hx = KokkosArray::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1 + j + 10 * i ;
    }
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );

  //------------------------------

  if ( print_flag ) {
    typename block_vector_type::HostMirror hy = KokkosArray::create_mirror( y );

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

  std::vector<double> perf(2) ;

  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter ;

  return perf ;
}

//----------------------------------------------------------------------------

template< typename ScalarType , class Device >
std::vector<double>
test_product_tensor_diagonal_matrix(
  const std::vector<int> & arg_var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool check )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type**,
                             KokkosArray::LayoutLeft ,
                             Device > block_vector_type ;

  //------------------------------

  typedef KokkosArray::BlockCrsMatrix< KokkosArray::SymmetricDiagonalSpec< Device > ,
                                  value_type , Device > matrix_type ;

  typedef typename matrix_type::graph_type  graph_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;
  const size_t fem_graph_length = unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate product tensor from variables' degrees

  const std::vector<unsigned> var_degree( arg_var_degree.begin() , arg_var_degree.end() );

  const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
    tensor( var_degree );

  const size_t stoch_length = tensor.bases_count();

  //------------------------------

  block_vector_type x = block_vector_type( "x" , stoch_length , fem_length );
  block_vector_type y = block_vector_type( "y" , stoch_length , fem_length );

  typename block_vector_type::HostMirror hx        = KokkosArray::create_mirror( x );
  typename block_vector_type::HostMirror hy_result = KokkosArray::create_mirror( y );

  for ( size_t iColFEM = 0 ;   iColFEM < fem_length ;   ++iColFEM ) {
  for ( size_t iColStoch = 0 ; iColStoch < stoch_length ; ++iColStoch ) {
    hx(iColStoch,iColFEM) =
      generate_vector_coefficient( fem_length , stoch_length ,
                                   iColFEM , iColStoch );
  }}

  KokkosArray::deep_copy( x , hx );

  //------------------------------
  // Generate CRS matrix of blocks with symmetric diagonal storage

  matrix_type matrix ;

  matrix.block  = KokkosArray::SymmetricDiagonalSpec< Device >( stoch_length );
  matrix.graph  = KokkosArray::create_crsarray<graph_type>( std::string("test product tensor graph") , fem_graph );
  matrix.values = block_vector_type( "matrix" , matrix.block.matrix_size() , fem_graph_length );

  {
    typename block_vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t iRowStoch = 0 ; iRowStoch < stoch_length ; ++iRowStoch ) {
      for ( size_t iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < fem_length ; ++iRowFEM ) {

        double y = 0 ;

        for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
          const size_t iColFEM = fem_graph[iRowFEM][iRowEntryFEM];

          for ( size_t iColStoch = 0 ; iColStoch < stoch_length ; ++iColStoch ) {

            const size_t offset = matrix.block.matrix_offset( iRowStoch , iColStoch );

            double value = 0 ;

            for ( size_t k = 0 ; k < stoch_length ; ++k ) {
              value += tensor( iRowStoch , iColStoch , k ) *
                       generate_matrix_coefficient( fem_length , stoch_length ,
                                                    iRowFEM , iColFEM , k );
            }

            y += value * hx(iColStoch,iColFEM);

            hM( offset , iEntryFEM ) = value ;
          }
        }

        hy_result( iRowStoch , iRowFEM ) = y ;
      }
    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1e-9*fem_graph_length*stoch_length*stoch_length / seconds_per_iter;

 //------------------------------
  // Verify result

  if (check)
  {
    const double tol = 1.0e-13 ;

    KokkosArray::deep_copy( hx , y );

    for ( size_t iRowFEM = 0 ; iRowFEM < fem_length ; ++iRowFEM ) {
      for ( size_t iRowStoch = 0 ; iRowStoch < stoch_length ; ++iRowStoch ) {
        const double mag   = std::abs( hy_result(iRowStoch,iRowFEM) );
        const double error = std::abs( hx(iRowStoch,iRowFEM) - hy_result(iRowStoch,iRowFEM) );
        if ( tol < error && tol < error / mag ) {
          std::cout << "test_product_tensor_diagonal_matrix error:"
                    << " y(" << iRowStoch << "," << iRowFEM << ") = " << hx(iRowStoch,iRowFEM)
                    << " , error = " << ( hx(iRowStoch,iRowFEM) - hy_result(iRowStoch,iRowFEM) )
                    << std::endl ;
        }
      }
    }
  }

  //------------------------------

  std::vector<double> perf(3);
  perf[0] = fem_length * stoch_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops;
  return perf;
}

//----------------------------------------------------------------------------
// Flatten to a plain CRS matrix
//
//  Outer DOF == fem
//  Inner DOF == stochastic

template< typename ScalarType , class Device >
std::vector<double>
test_product_flat_commuted_matrix(
  const std::vector<int> & arg_var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool check )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type[] , Device > vector_type ;

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;

  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Tensor evaluation:

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
  // Generate flattened graph with FEM outer and stochastic inner

  const size_t flat_length = fem_length * stoch_length ;

  std::vector< std::vector<size_t> > flat_graph( flat_length );

  for ( size_t iOuterRow = 0 ; iOuterRow < fem_length ; ++iOuterRow ) {

    const size_t iOuterRowNZ = fem_graph[iOuterRow].size();

    for ( size_t iInnerRow = 0 ; iInnerRow < stoch_length ; ++iInnerRow ) {

      const size_t iInnerRowNZ = stoch_graph[ iInnerRow ].size(); ;
      const size_t iFlatRowNZ  = iOuterRowNZ * iInnerRowNZ ;
      const size_t iFlatRow    = iInnerRow + iOuterRow * stoch_length ;

      flat_graph[iFlatRow].resize( iFlatRowNZ );

      size_t iFlatEntry = 0 ;

      for ( size_t iOuterEntry = 0 ; iOuterEntry < iOuterRowNZ ; ++iOuterEntry ) {

        const size_t iOuterCol = fem_graph[iOuterRow][iOuterEntry];

        for ( size_t iInnerEntry = 0 ; iInnerEntry < iInnerRowNZ ; ++iInnerEntry ) {

          const size_t iInnerCol   = stoch_graph[iInnerRow][iInnerEntry] ;
          const size_t iFlatColumn = iInnerCol + iOuterCol * stoch_length ;

          flat_graph[iFlatRow][iFlatEntry] = iFlatColumn ;

          ++iFlatEntry ;
        }
      }
    }
  }

  //------------------------------

  vector_type x = vector_type( "x" , flat_length );
  vector_type y = vector_type( "y" , flat_length );

  typename vector_type::HostMirror hx        = KokkosArray::create_mirror( x );
  typename vector_type::HostMirror hy_result = KokkosArray::create_mirror( y );

  for ( size_t iCol = 0 ; iCol < flat_length ; ++iCol ) {
    const size_t iColFEM   = iCol / stoch_length ;
    const size_t iColStoch = iCol % stoch_length ;

    hx(iCol) = generate_vector_coefficient( fem_length , stoch_length ,
                                            iColFEM , iColStoch );
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.dimension(0);

  matrix.values = vector_type( "matrix" , flat_graph_length );
  {
    typename vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t iRow = 0 , iEntry = 0 ; iRow < flat_length ; ++iRow ) {
      const size_t iRowFEM   = iRow / stoch_length ;
      const size_t iRowStoch = iRow % stoch_length ;

      double y = 0 ;

      for ( size_t iRowEntry = 0 ; iRowEntry < flat_graph[ iRow ].size() ; ++iRowEntry , ++iEntry ) {
        const size_t iCol = flat_graph[ iRow ][ iRowEntry ];
        const size_t iColFEM   = iCol / stoch_length ;
        const size_t iColStoch = iCol % stoch_length ;

        double value = 0 ;
        for ( unsigned k = 0 ; k < stoch_length ; ++k ) {
          const double A_fem_k = generate_matrix_coefficient( fem_length , stoch_length ,
                                                              iRowFEM, iColFEM, k );
          value += tensor(iRowStoch,iColStoch,k) * A_fem_k ;
        }
        hM( iEntry ) = value ;

        y += value * hx( iCol );
      }

      hy_result( iRow ) = y ;
    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  //KokkosArray::write_matrix_market(matrix, "flat_commuted.mm");

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1e-9*flat_graph_length / seconds_per_iter;

  //------------------------------
  // Verify result

  if (check)
  {
    const double tol = 1.0e-13 ;

    KokkosArray::deep_copy( hx , y );

    for ( size_t i = 0 ; i < flat_length ; ++i ) {
      const double mag   = std::abs( hy_result(i) );
      const double error = std::abs( hx(i) - hy_result(i) );
      if ( tol < error && tol < error / mag ) {
        std::cout << "test_product_flat_commuted_matrix error:"
                  << " y[" << i << "] = " << hx(i)
                  << " , error = " << ( hx(i) - hy_result(i) )
                  << std::endl ;
      }
    }
  }

  std::vector<double> perf(4);
  perf[0] = flat_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops;
  perf[3] = flat_graph_length ;
  return perf;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Flatten to a plain CRS matrix
//
//  Outer DOF == stochastic
//  Inner DOF == fem

template< typename ScalarType , class Device >
std::vector<double>
test_product_flat_original_matrix(
  const std::vector<int> & arg_var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool check )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type[] , Device > vector_type ;

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;

  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Tensor evaluation:

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
  // Generate flattened graph with stochastic outer and FEM inner

  const size_t flat_length = fem_length * stoch_length ;

  std::vector< std::vector<size_t> > flat_graph( flat_length );

  for ( size_t iOuterRow = 0 ; iOuterRow < stoch_length ; ++iOuterRow ) {

    const size_t iOuterRowNZ = stoch_graph[iOuterRow].size();

    for ( size_t iInnerRow = 0 ; iInnerRow < fem_length ; ++iInnerRow ) {

      const size_t iInnerRowNZ = fem_graph[iInnerRow].size();
      const size_t iFlatRowNZ  = iOuterRowNZ * iInnerRowNZ ;
      const size_t iFlatRow    = iInnerRow + iOuterRow * fem_length ;

      flat_graph[iFlatRow].resize( iFlatRowNZ );

      size_t iFlatEntry = 0 ;

      for ( size_t iOuterEntry = 0 ; iOuterEntry < iOuterRowNZ ; ++iOuterEntry ) {

        const size_t iOuterCol = stoch_graph[ iOuterRow ][ iOuterEntry ];

        for ( size_t iInnerEntry = 0 ; iInnerEntry < iInnerRowNZ ; ++iInnerEntry ) {

          const size_t iInnerCol   = fem_graph[ iInnerRow][iInnerEntry];
          const size_t iFlatColumn = iInnerCol + iOuterCol * fem_length ;

          flat_graph[iFlatRow][iFlatEntry] = iFlatColumn ;
          ++iFlatEntry ;
        }
      }
    }
  }

  //------------------------------

  vector_type x = vector_type( "x" , flat_length );
  vector_type y = vector_type( "y" , flat_length );

  typename vector_type::HostMirror hx        = KokkosArray::create_mirror( x );
  typename vector_type::HostMirror hy_result = KokkosArray::create_mirror( y );

  for ( size_t iCol = 0 ; iCol < flat_length ; ++iCol ) {
    const size_t iColStoch = iCol / fem_length ;
    const size_t iColFEM   = iCol % fem_length ;

    hx(iCol) = generate_vector_coefficient( fem_length , stoch_length ,
                                            iColFEM , iColStoch );
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.dimension(0);

  matrix.values = vector_type( "matrix" , flat_graph_length );
  {
    typename vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t iRow = 0 , iEntry = 0 ; iRow < flat_length ; ++iRow ) {
      const size_t iRowStoch = iRow / fem_length ;
      const size_t iRowFEM   = iRow % fem_length ;

      double y = 0 ;

      for ( size_t iRowEntry = 0 ; iRowEntry < flat_graph[ iRow ].size() ; ++iRowEntry , ++iEntry ) {
        const size_t iCol = flat_graph[ iRow ][ iRowEntry ];
        const size_t iColStoch = iCol / fem_length ;
        const size_t iColFEM   = iCol % fem_length ;

        double value = 0 ;
        for ( unsigned k = 0 ; k < stoch_length ; ++k ) {
          const double A_fem_k =
            generate_matrix_coefficient( fem_length , stoch_length ,
                                         iRowFEM , iColFEM , k );
          value += tensor(iRowStoch,iColStoch,k) * A_fem_k ;
        }
        hM( iEntry ) = value ;

        y += value * hx( iCol );
      }

      hy_result( iRow ) = y ;
    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  //KokkosArray::write_matrix_market(matrix, "flat_original.mm");

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1e-9*flat_graph_length / seconds_per_iter;

  //------------------------------
  // Verify result

  if (check)
  {
    const double tol = 1.0e-13 ;

    KokkosArray::deep_copy( hx , y );

    for ( size_t i = 0 ; i < flat_length ; ++i ) {
      const double mag   = std::abs( hy_result(i) );
      const double error = std::abs( hx(i) - hy_result(i) );
      if ( tol < error && tol < error / mag ) {
        std::cout << "test_product_flat_original_matrix error:"
                  << " y[" << i << "] = " << hx(i)
                  << " , error = " << ( hx(i) - hy_result(i) )
                  << std::endl ;
      }
    }
  }

  std::vector<double> perf(4);
  perf[0] = flat_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops;
  perf[3] = flat_graph_length ;
  return perf;
}

//----------------------------------------------------------------------------
// Outer original matrix-free block algorithm
#ifdef HAVE_KOKKOSARRAY_STOKHOS
template< typename ScalarType , class Device >
std::vector<double>
test_original_matrix_free_block(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag ,
  const bool test_block ,
  const bool check )
{
  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::LegendreBasis<int,value_type> basis_type;
  typedef Stokhos::CompletePolynomialBasis<int,value_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL); 
  for (size_t i=0; i<num_KL; i++)
    bases[i] = Teuchos::rcp(new basis_type(var_degree[i],true));
  RCP<const product_basis_type> basis = 
    rcp(new product_basis_type(bases, 1e-12));
  const size_t outer_length = basis->size();
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor(outer_length);

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  
  typedef KokkosArray::View<value_type**, KokkosArray::LayoutLeft, Device> multi_vec_type ;
  typedef KokkosArray::View<value_type[], Device> vec_type ;

  std::vector<matrix_type> matrix( outer_length ) ;
  multi_vec_type x( "x" , inner_length , outer_length ) ;
  multi_vec_type y( "y" , inner_length , outer_length ) ;
  multi_vec_type tmp( "tmp" , inner_length ,  outer_length ) ;

  for (size_t block=0; block<outer_length; ++block) {
    matrix[block].graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );
    const size_t graph_length = matrix[block].graph.entries.dimension(0);
    matrix[block].values = vec_type( "matrix" , graph_length );

    typename vec_type::HostMirror hM =
      KokkosArray::create_mirror( matrix[block].values );
    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }
    KokkosArray::deep_copy( matrix[block].values , hM );
  }
  
  typename multi_vec_type::HostMirror hx =
    KokkosArray::create_mirror( x );
  for (size_t block=0; block<outer_length; ++block) {
    for ( size_t i = 0 ; i < inner_length ; ++i ) {
      hx( i , block) = 1 + i ;
    }
  }
  KokkosArray::deep_copy( x , hx );

  KokkosArray::Impl::Timer clock ;
  int n_apply = 0;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    // Original matrix-free multiply algorithm using a block apply
    n_apply = 0;
    typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
    typename Cijk_type::k_iterator k_end = Cijk->k_end();
    for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int nj = Cijk->num_j(k_it);
      if (nj > 0) {
	int k = index(k_it);
	typename Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
	typename Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
	std::vector<int> j_indices(nj);
	int jdx = 0;
	for (typename Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; 
	     ++j_it) {
	  j_indices[jdx++] = index(j_it);
	}
        KokkosArray::multiply( matrix[k] , x , tmp, j_indices , test_block );
        n_apply += nj;
	jdx = 0;
	for (typename Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; 
	     ++j_it) {
	  const vec_type tmp_view( tmp , jdx++ );
	  typename Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
	  typename Cijk_type::kji_iterator i_end =  Cijk->i_end(j_it);
	  for (typename Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; 
	       ++i_it) {
	    int i = index(i_it);
	    value_type c = value(i_it);
	    const vec_type y_view( y , i );
	    KokkosArray::update( value_type(1.0) , y_view , c , tmp_view );
	  }
	}
      }
    }

  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1.0e-9*n_apply*matrix[0].graph.entries.dimension(0);

  //------------------------------

  if ( print_flag ) {
    std::cout << std::endl << "test_product_flat_matrix"
              << std::endl ;
    typename multi_vec_type::HostMirror hy = KokkosArray::create_mirror( y );
    KokkosArray::deep_copy( hy , y );
    std::cout << "hy = " << std::endl ;
    for ( size_t i = 0 ; i < inner_length ; ++i ) {
      for ( size_t j = 0 ; j < outer_length ; ++j ) {
	std::cout << " " << hy(i,j);
      }
      std::cout << std::endl ;
    }
  }

  std::vector<double> perf(4);
  perf[0] = outer_length * inner_length;
  perf[1] = seconds_per_iter ;
  perf[2] = flops/seconds_per_iter;
  perf[3] = flops;

  return perf;
}

template< typename ScalarType , class Device >
std::vector<double>
test_original_matrix_free_vec(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag ,
  const bool test_block ,
  const bool check )
{
  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::LegendreBasis<int,value_type> basis_type;
  typedef Stokhos::CompletePolynomialBasis<int,value_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL); 
  for (size_t i=0; i<num_KL; i++)
    bases[i] = Teuchos::rcp(new basis_type(var_degree[i],true));
  RCP<const product_basis_type> basis = 
    rcp(new product_basis_type(
	  bases, ScalarTolerances<value_type>::sparse_cijk_tol()));
  const size_t outer_length = basis->size();
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor(outer_length);

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  
  typedef KokkosArray::View<value_type[],Device> vec_type ;

  std::vector<matrix_type> matrix( outer_length ) ;
  std::vector<vec_type> x( outer_length ) ;
  std::vector<vec_type> y( outer_length ) ;
  std::vector<vec_type> tmp( outer_length ) ;

  for (size_t block=0; block<outer_length; ++block) {
    matrix[block].graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );

    const size_t graph_length = matrix[block].graph.entries.dimension(0);

    matrix[block].values = vec_type( "matrix" , graph_length );

    x[block]   = vec_type( "x" , inner_length );
    y[block]   = vec_type( "y" , inner_length );
    tmp[block] = vec_type( "tmp" , inner_length );

    typename vec_type::HostMirror hM =
      KokkosArray::create_mirror( matrix[block].values );

    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }

    KokkosArray::deep_copy( matrix[block].values , hM );

    typename vec_type::HostMirror hx =
      KokkosArray::create_mirror( x[block] );

    for ( size_t i = 0 ; i < inner_length ; ++i ) {
      hx(i) = 1 + i ;
    }

    KokkosArray::deep_copy( x[block] , hx );
  }
  

  KokkosArray::Impl::Timer clock ;
  int n_apply = 0;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    // Original matrix-free multiply algorithm using a block apply
    n_apply = 0;
    typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
    typename Cijk_type::k_iterator k_end = Cijk->k_end();
    for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int nj = Cijk->num_j(k_it);
      if (nj > 0) {
	int k = index(k_it);
	typename Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
	typename Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
	std::vector<vec_type> xx(nj), yy(nj);
	int jdx = 0;
	for (typename Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; 
	     ++j_it) {
	  int j = index(j_it);
	  xx[jdx] = x[j];
	  yy[jdx] = tmp[j];
	  jdx++;
	}
        KokkosArray::multiply( matrix[k] , xx , yy, test_block );
        n_apply += nj;
	jdx = 0;
	for (typename Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; 
	     ++j_it) {
	  typename Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
	  typename Cijk_type::kji_iterator i_end =  Cijk->i_end(j_it);
	  for (typename Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; 
	       ++i_it) {
	    int i = index(i_it);
	    value_type c = value(i_it);
	    KokkosArray::update( value_type(1.0) , y[i] , c , yy[jdx] );
	  }
	  jdx++;
	}
      }
    }

  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1.0e-9*n_apply*matrix[0].graph.entries.dimension(0);

  //------------------------------

  if ( print_flag ) {
    std::cout << std::endl << "test_product_flat_matrix"
              << std::endl ;
    for ( size_t i = 0 ; i < outer_length ; ++i ) {
      typename vec_type::HostMirror hy = KokkosArray::create_mirror( y[i] );

      KokkosArray::deep_copy( hy , y[i] );

      std::cout << "hy(:," << i << ") =" ;
      for ( size_t j = 0 ; j < inner_length ; ++j ) {
        std::cout << " " << hy(j);
      }
      std::cout << std::endl ;
    }
  }

  std::vector<double> perf(4);
  perf[0] = outer_length * inner_length;
  perf[1] = seconds_per_iter ;
  perf[2] = flops/seconds_per_iter;
  perf[3] = flops;

  return perf;
}
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class Scalar, class Device >
void performance_test_driver_all( const int pdeg ,
				  const int minvar ,
				  const int maxvar ,
				  const int nGrid ,
				  const int nIter ,
				  const bool print ,
				  const bool test_block ,
				  const bool check )
{
  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_nonzeros = unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  std::cout.precision(8);

  //------------------------------

  std::cout << std::endl << "\"FEM NNZ = " << fem_nonzeros << "\"" << std::endl;

  std::cout << std::endl
	    << "\"#nGrid\" , "
            << "\"#Variable\" , \"PolyDegree\" , \"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
	    << "\"Original-Flat MXV-Time\" , "
	    << "\"Original-Flat MXV-Speedup\" , "
            << "\"Original-Flat MXV-GFLOPS\" , "
	    << "\"Commuted-Flat MXV-Speedup\" , "
            << "\"Commuted-Flat MXV-GFLOPS\" , "
    //<< "\"Commuted-Flat NNZ\" , "
#ifdef HAVE_KOKKOSARRAY_STOKHOS
            << "\"Original-Matrix-Free-Block MXV-Speedup\" , "
            << "\"Original-Matrix-Free-Block MXV-GFLOPS\" , "
#endif
	    << "\"Block-Diagonal MXV-Speedup\" , "
            << "\"Block-Diagonal MXV-GFLOPS\" , "
    //<< "\"Block-Coord-Tensor MXV-Speedup\" , "
	    << "\"Block-Crs-Tensor MXV-Speedup\" , "
	    << "\"Block-Crs-Tensor MXV-Time\" , "
#ifdef HAVE_KOKKOSARRAY_STOKHOS
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , "
#endif
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {

    std::vector<int> var_degree( nvar , pdeg );

    //------------------------------
    // Tensor evaluation:

    const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
      tensor( std::vector<unsigned>( nvar , unsigned(pdeg) ) );

    const size_t stoch_length = tensor.bases_count();
    size_t stoch_nonzero = 0 ;

    for ( size_t i = 0 ; i < stoch_length ; ++i ) {
    for ( size_t j = i ; j < stoch_length ; ++j ) {
    for ( size_t k = j ; k < stoch_length ; ++k ) {
      if ( tensor.is_non_zero(i,j,k) ) ++stoch_nonzero ;
    }}}

    //------------------------------

    const std::vector<double> perf_matrix =
      test_product_tensor_diagonal_matrix<Scalar,Device>( var_degree , nGrid , nIter , check );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_legendre<Scalar,Scalar,Scalar,Device>( var_degree , nGrid , nIter , check );

    const std::vector<double> perf_flat_commuted =
      test_product_flat_commuted_matrix<Scalar,Device>( var_degree , nGrid , nIter , check );

    const std::vector<double> perf_flat_original =
      test_product_flat_original_matrix<Scalar,Device>( var_degree , nGrid , nIter , check );

    if ( perf_flat_commuted[0] != perf_flat_original[0] ||
         perf_flat_commuted[3] != perf_flat_original[3] ) {
      std::cout << "ERROR: Original and commuted matrix sizes do not match" 
                << std::endl
                << "  original size = " << perf_flat_original[0]
                << " , nonzero = " << perf_flat_original[3]
                << std::endl
                << "  commuted size = " << perf_flat_commuted[0]
                << " , nonzero = " << perf_flat_commuted[3]
                << std::endl ;
    }

#ifdef HAVE_KOKKOSARRAY_STOKHOS
    const std::vector<double> perf_original_mat_free_block =
      test_original_matrix_free_vec<Scalar,Device>( var_degree , nGrid , nIter , print , test_block , check );
#endif

    std::cout << nGrid << " , " << nvar << " , " << pdeg << " , "
	      << tensor.bases_count() << " , "
	      << stoch_nonzero << " , "
	      << perf_flat_original[0] << " , "
	      << perf_flat_original[1] << " , "
	      << perf_flat_original[1] / perf_flat_original[1] << " , "
              << perf_flat_original[2] << " , "
	      << perf_flat_original[1] / perf_flat_commuted[1] << " , "
              << perf_flat_commuted[2] << " , "
      //<< perf_flat_commuted[3] << " , "
#ifdef HAVE_KOKKOSARRAY_STOKHOS
	      << perf_flat_original[1] / perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
#endif
	      << perf_flat_original[1] / perf_matrix[1] << " , "
              << perf_matrix[2] << " , "
      //<< perf_flat_original.second / perf_tensor.second << " , "
	      << perf_flat_original[1] / perf_crs_tensor[1] << " , "
	      << perf_crs_tensor[1] << " , "
#ifdef HAVE_KOKKOSARRAY_STOKHOS
	      << perf_original_mat_free_block[3] / perf_crs_tensor[1] << " , "
#endif
	      << std::endl ;
  }

  //------------------------------
}

#ifdef HAVE_KOKKOSARRAY_STOKHOS
template< class Scalar, class Device >
void performance_test_driver_poly( const int pdeg ,
				   const int minvar ,
				   const int maxvar ,
				   const int nGrid ,
				   const int nIter ,
				   const bool print ,
				   const bool test_block ,
				   const bool check )
{
  typedef KokkosArray::NormalizedLegendrePolynomialBases<8> polynomial ;
  typedef KokkosArray::StochasticProductTensor< Scalar , polynomial , Device , KokkosArray::CrsProductTensor > tensor_type ;

  std::cout.precision(8);

  //------------------------------

  std::vector< std::vector<size_t> > fem_graph ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );
  std::cout << std::endl << "\"FEM NNZ = " << graph_length << "\"" << std::endl;

  std::cout << std::endl
	    << "\"#nGrid\" , "
            << "\"#Variable\" , \"PolyDegree\" , \"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
	    << "\"Original-Matrix-Free-Block-MXV-Time\" , "
	    << "\"Original-Matrix-Free-Block-MXV-Speedup\" , "
            << "\"Original-Matrix-Free-Block-MXV-GFLOPS\" , "
	    << "\"Block-Legendre-Tensor MXV-Time\" , "
	    << "\"Block-Legendre-Tensor MXV-Speedup\" , "
            << "\"Block-Legendre-Tensor MXV-GFLOPS\" , "
	    << "\"Block-Crs-Tensor MXV-Time\" , "
	    << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , "
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {
    std::vector<int> var_degree( nvar , pdeg );

    const tensor_type tensor = KokkosArray::create_product_tensor< tensor_type >( var_degree );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_matrix<Scalar,Device,KokkosArray::CrsProductTensor>( var_degree , nGrid , nIter , check );

    const std::vector<double> perf_legendre_tensor =
      test_product_tensor_legendre<Scalar,Scalar,Scalar,Device>( var_degree , nGrid , nIter , check );

    const std::vector<double> perf_original_mat_free_block =
      test_original_matrix_free_vec<Scalar,Device>( var_degree , nGrid , nIter , print , test_block , check );

    std::cout << nGrid << " , "
	      << nvar << " , " << pdeg << " , "
	      << tensor.dimension() << " , "
	      << tensor.tensor().entry_count() << " , "
	      << perf_original_mat_free_block[0] << " , "
	      << perf_original_mat_free_block[1] << " , "
	      << perf_original_mat_free_block[1] / perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
	      << perf_legendre_tensor[1] << " , "
	      << perf_original_mat_free_block[1] / perf_legendre_tensor[1] << " , "
              << perf_original_mat_free_block[3] / perf_legendre_tensor[1] << " , "
	      << perf_crs_tensor[1] << " , "
	      << perf_original_mat_free_block[1] / perf_crs_tensor[1] << " , "
              << perf_original_mat_free_block[3] / perf_crs_tensor[1] << " , "

	      << std::endl ;
  }

  //------------------------------
}
#endif

template< class Scalar, class Device >
struct performance_test_driver {
  static void run(bool test_flat, bool test_orig, bool test_block,
		  bool check) {}
};

//----------------------------------------------------------------------------

}


