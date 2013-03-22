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
#include <utility>
#include <cmath>
#include <iostream>

#include "KokkosArray_Host.hpp"
#include "impl/KokkosArray_Timer.hpp"

#include "Stokhos_CrsMatrix.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_StochasticProductTensor.hpp"
#include "Stokhos_CrsProductTensor.hpp"
#include "Stokhos_FlatSparse3Tensor.hpp"
#include "Stokhos_FlatSparse3Tensor_kji.hpp"

#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_JacobiBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_TotalOrderBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

namespace unit_test {

template< typename IntType >
inline
IntType map_fem_graph_coord( const IntType & N ,
                             const IntType & i ,
                             const IntType & j ,
                             const IntType & k )
{ 
  return k + N * ( j + N * i );
} 

inline
size_t generate_fem_graph( size_t N , 
			   std::vector< std::vector<size_t> > & graph )
{
  graph.resize( N * N * N , std::vector<size_t>() );

  size_t total = 0 ;

  for ( int i = 0 ; i < (int) N ; ++i ) {
  for ( int j = 0 ; j < (int) N ; ++j ) {
  for ( int k = 0 ; k < (int) N ; ++k ) {

    const size_t row = map_fem_graph_coord((int)N,i,j,k);

    graph[row].reserve(27);

    for ( int ii = -1 ; ii < 2 ; ++ii ) {
    for ( int jj = -1 ; jj < 2 ; ++jj ) {
    for ( int kk = -1 ; kk < 2 ; ++kk ) {
      if ( 0 <= i + ii && i + ii < (int) N &&
           0 <= j + jj && j + jj < (int) N &&
           0 <= k + kk && k + kk < (int) N ) {
        size_t col = map_fem_graph_coord((int)N,i+ii,j+jj,k+kk);

        graph[row].push_back(col);
      }
    }}}
    total += graph[row].size();
  }}}

  return total ;
}

template <typename Scalar> struct ScalarTolerances {};

template <> struct ScalarTolerances<float> {
  typedef float scalar_type;
  static scalar_type sparse_cijk_tol() { return 1e-5; }
};

template <> struct ScalarTolerances<double> {
  typedef double scalar_type;
  static scalar_type sparse_cijk_tol() { return 1e-12; }
};

template< typename ScalarType , typename TensorType, class Device >
std::vector<double>
test_product_tensor_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type** ,
                             KokkosArray::LayoutLeft ,
                             Device > block_vector_type ;

   typedef Stokhos::StochasticProductTensor< value_type , TensorType , Device > tensor_type ;

    typedef Stokhos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::JacobiBasis<int,value_type> basis_type;
  typedef Stokhos::LexographicLess<Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,value_type,order_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL); 
  for (size_t i=0; i<num_KL; i++) {
    if (symmetric)
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,1.0,true));
    else
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,2.0,true));
  }
  RCP<const product_basis_type> basis = 
    rcp(new product_basis_type(
	  bases, ScalarTolerances<value_type>::sparse_cijk_tol()));
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor();

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  matrix_type matrix ;

  matrix.block = 
    Stokhos::create_stochastic_product_tensor< TensorType >( *basis,
							     *Cijk );
  matrix.graph = KokkosArray::create_crsarray<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();

  matrix.values = block_vector_type( "matrix" , inner_length , graph_length );

  block_vector_type x = block_vector_type( "x" , inner_length , outer_length );
  block_vector_type y = block_vector_type( "y" , inner_length , outer_length );

  typename block_vector_type::HostMirror hM = 
    KokkosArray::create_mirror( matrix.values );

  for ( size_t i=0 ; i < graph_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1.0;
    }
  }
  
  KokkosArray::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:
  
  typename block_vector_type::HostMirror hx = KokkosArray::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1.0 ;
    }
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = matrix.block.tensor().num_flops();
  const double flops = 1.0e-9*graph_length*flops_per_block;

  // std::cout << "tensor: flops = " << flops 
  // 	    << " time = " << seconds_per_iter << std::endl;

  std::vector<double> perf(5) ;

  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter ;
  perf[2] = flops / seconds_per_iter;
  perf[3] = matrix.block.tensor().entry_count();
  perf[4] = inner_length ;

  return perf ;
}

template< typename ScalarType , class Device , class SparseMatOps >
std::vector<double>
test_original_matrix_free_vec(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool test_block ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::JacobiBasis<int,value_type> basis_type;
  typedef Stokhos::LexographicLess<Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,value_type,order_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL); 
  for (size_t i=0; i<num_KL; i++) {
    if (symmetric)
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,1.0,true));
    else
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,2.0,true));
  }
  RCP<const product_basis_type> basis = 
    rcp(new product_basis_type(
	  bases, ScalarTolerances<value_type>::sparse_cijk_tol()));
  const size_t outer_length = basis->size();
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor();

  //------------------------------

  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = 
    unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  
  typedef KokkosArray::View<value_type[],Device> vec_type ;

  std::vector<matrix_type> matrix( outer_length ) ;
  std::vector<vec_type> x( outer_length ) ;
  std::vector<vec_type> y( outer_length ) ;
  std::vector<vec_type> tmp( outer_length ) ;

  for (size_t block=0; block<outer_length; ++block) {
    matrix[block].graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );

    matrix[block].values = vec_type( "matrix" , graph_length );

    x[block]   = vec_type( "x" , inner_length );
    y[block]   = vec_type( "y" , inner_length );
    tmp[block] = vec_type( "tmp" , inner_length );

    typename vec_type::HostMirror hM =
      KokkosArray::create_mirror( matrix[block].values );

    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      hM(i) = 1.0 ;
    }

    KokkosArray::deep_copy( matrix[block].values , hM );

    typename vec_type::HostMirror hx =
      KokkosArray::create_mirror( x[block] );

    for ( size_t i = 0 ; i < inner_length ; ++i ) {
      hx(i) = 1.0 ;
    }

    KokkosArray::deep_copy( x[block] , hx );
  }
  
  SparseMatOps smo;
  KokkosArray::Impl::Timer clock ;
  int n_apply = 0;
  int n_add = 0;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    // Original matrix-free multiply algorithm using a block apply
    n_apply = 0;
    n_add = 0;
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
        Stokhos::multiply( matrix[k] , xx , yy, test_block ,smo );
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
	    Stokhos::update( value_type(1.0) , y[i] , c , yy[jdx] );
	    ++n_add;
	  }
	  jdx++;
	}
      }
    }

  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 1.0e-9*(2.0*static_cast<double>(n_apply)*graph_length+
			       static_cast<double>(n_add)*inner_length);

  // std::cout << "mat-free: flops = " << flops 
  // 	    << " time = " << seconds_per_iter << std::endl;

  std::vector<double> perf(4);
  perf[0] = outer_length * inner_length;
  perf[1] = seconds_per_iter ;
  perf[2] = flops/seconds_per_iter;
  perf[3] = flops;

  return perf;
}

template< class Scalar, class Device , class SparseMatOps >
void performance_test_driver_poly( const int pdeg ,
				   const int minvar ,
				   const int maxvar ,
				   const int nGrid ,
				   const int nIter ,
				   const bool test_block ,
				   const bool symmetric )
{
  bool do_flat_sparse = 
    KokkosArray::Impl::is_same<Device,KokkosArray::Host>::value ;

  std::cout.precision(8);

  //------------------------------

  std::vector< std::vector<size_t> > fem_graph ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );
  std::cout << std::endl << "\"FEM NNZ = " << graph_length << "\"" << std::endl;

  std::cout << std::endl
	    << "\"#nGrid\" , "
            << "\"#Variable\" , "
	    << "\"PolyDegree\" , "
	    << "\"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
	    << "\"Original-Matrix-Free-Block-MXV-Time\" , "
	    << "\"Original-Matrix-Free-Block-MXV-Speedup\" , "
            << "\"Original-Matrix-Free-Block-MXV-GFLOPS\" , "
	    << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , ";
  if (do_flat_sparse)
    std::cout << "\"Block-Flat-Sparse-3-Tensor MXV-Speedup\" , "
	      << "\"Block-Flat-Sparse-3-Tensor MXV-GFLOPS\" , ";
  std::cout << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {
    std::vector<int> var_degree( nvar , pdeg );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_matrix<Scalar,Stokhos::CrsProductTensor<Scalar,Device>,Device>( 
    	var_degree , nGrid , nIter , symmetric );

    std::vector<double> perf_flat_sparse_3_tensor;
    if (do_flat_sparse) {
      perf_flat_sparse_3_tensor = 
	test_product_tensor_matrix<Scalar,Stokhos::FlatSparse3Tensor_kji<Scalar,Device>,Device>( var_degree , nGrid , nIter , symmetric );
    }

    const std::vector<double> perf_original_mat_free_block =
      test_original_matrix_free_vec<Scalar,Device,SparseMatOps>( 
    	var_degree , nGrid , nIter , test_block , symmetric );

    std::cout << nGrid << " , "
	      << nvar << " , " 
	      << pdeg << " , "
	      << perf_crs_tensor[4] << " , "
	      << perf_crs_tensor[3] << " , "
	      << perf_original_mat_free_block[0] << " , "
	      << perf_original_mat_free_block[1] << " , "
	      << perf_original_mat_free_block[1] / perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
	      << perf_original_mat_free_block[1] / perf_crs_tensor[1] << " , "
              << perf_crs_tensor[2] << " , ";
    if (do_flat_sparse)
      std::cout << perf_original_mat_free_block[1] / perf_flat_sparse_3_tensor[1] << " , "
		<< perf_flat_sparse_3_tensor[2] << " , ";

    std::cout << std::endl ;
  }

  //------------------------------
}

template< class Scalar, class Device >
struct performance_test_driver {
  static void run(bool test_flat, bool test_orig, bool test_block, bool symmetric) {}
};

//----------------------------------------------------------------------------

}


