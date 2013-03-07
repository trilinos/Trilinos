// $Id$ 
// $Source$ 
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

#ifndef STOKHOS_KOKKOS_ARRAY_KERNELS_UNIT_TEST_HPP
#define STOKHOS_KOKKOS_ARRAY_KERNELS_UNIT_TEST_HPP

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_TestForException.hpp"

#include "Stokhos_Epetra.hpp"
#include "EpetraExt_BlockUtility.h"
#include "Stokhos_UnitTestHelpers.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "KokkosArray_config.h"
#include <KokkosArray_Host.hpp>
#include <KokkosArray_ProductTensor.hpp>
#include <KokkosArray_LegendrePolynomial.hpp>
#include <KokkosArray_SymmetricDiagonalSpec.hpp>
#include <KokkosArray_StochasticProductTensor.hpp>
#include <KokkosArray_CrsMatrix.hpp>
#include <KokkosArray_CrsProductTensorLegendre.hpp>
#include <KokkosArray_BlockCrsMatrix.hpp>

namespace KokkosArrayKernelsUnitTest {

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;
  using Teuchos::ParameterList;
  
  template< typename IntType >
  inline
  IntType map_fem_graph_coord( const IntType & N ,
			       const IntType & i ,
			       const IntType & j ,
			       const IntType & k )
  { 
    return k + N * ( j + N * i );
  } 
  
  template < typename ordinal >
  inline
  ordinal generate_fem_graph( ordinal N , 
			      std::vector< std::vector<ordinal> > & graph )
  {
    graph.resize( N * N * N , std::vector<ordinal>() );
    
    ordinal total = 0 ;
    
    for ( int i = 0 ; i < (int) N ; ++i ) {
      for ( int j = 0 ; j < (int) N ; ++j ) {
	for ( int k = 0 ; k < (int) N ; ++k ) {
	  
	  const ordinal row = map_fem_graph_coord((int)N,i,j,k);
	  
	  graph[row].reserve(27);
	  
	  for ( int ii = -1 ; ii < 2 ; ++ii ) {
	    for ( int jj = -1 ; jj < 2 ; ++jj ) {
	      for ( int kk = -1 ; kk < 2 ; ++kk ) {
		if ( 0 <= i + ii && i + ii < (int) N &&
		     0 <= j + jj && j + jj < (int) N &&
		     0 <= k + kk && k + kk < (int) N ) {
		  ordinal col = map_fem_graph_coord((int)N,i+ii,j+jj,k+kk);
		  
		  graph[row].push_back(col);
		}
	      }}}
	  total += graph[row].size();
	}}}
    
    return total ;
  }

  template <typename scalar, typename ordinal>
  inline
  scalar generate_matrix_coefficient( const ordinal nFEM ,
				      const ordinal nStoch ,
				      const ordinal iRowFEM ,
				      const ordinal iColFEM ,
				      const ordinal iStoch )
  {
    const scalar A_fem = ( 10.0 + scalar(iRowFEM) / scalar(nFEM) ) +
      (  5.0 + scalar(iColFEM) / scalar(nFEM) );
    
    const scalar A_stoch = ( 1.0 + scalar(iStoch) / scalar(nStoch) );

    return A_fem + A_stoch ;
  }

  template <typename scalar, typename ordinal>
  inline
  scalar generate_vector_coefficient( const ordinal nFEM ,
				      const ordinal nStoch ,
				      const ordinal iColFEM ,
				      const ordinal iStoch )
  {
    const scalar X_fem = 100.0 + scalar(iColFEM) / scalar(nFEM);
    const scalar X_stoch =  1.0 + scalar(iStoch) / scalar(nStoch);
    return X_fem + X_stoch ;
  }

  struct UnitTestSetup {
    typedef double value_type ;
    typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
    typedef Stokhos::LegendreBasis<int,value_type> basis_type;
    typedef Stokhos::CompletePolynomialBasis<int,value_type> product_basis_type;
    typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

    int p, d, nGrid, fem_length, stoch_length, fem_graph_length;
    double rel_tol, abs_tol;
    std::vector< std::vector<int> > fem_graph ;
    RCP< product_basis_type> basis;
    RCP<Cijk_type> Cijk;
    RCP<Stokhos::EpetraVectorOrthogPoly> sg_x, sg_y;
    RCP<Stokhos::ProductEpetraVector> sg_y_commuted;
    Teuchos::Array<int> perm, inv_perm, perm2, inv_perm2;

    // Can't be a constructor because MPI will not be initialized
    void setup() {

      p = 5;
      d = 2;
      nGrid = 5;
      rel_tol = 1e-12;
      abs_tol = 1e-12;

      // Create a communicator for Epetra objects
      RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
      globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
      globalComm = rcp(new Epetra_SerialComm);
#endif
      //int MyPID = globalComm->MyPID();

      //------------------------------
      // Generate FEM graph:
      
      fem_length = nGrid * nGrid * nGrid ;
      fem_graph_length = generate_fem_graph( nGrid , fem_graph );

      // Create Stochastic Galerkin basis and expansion
      Array< RCP<const abstract_basis_type> > bases(d); 
      for (int i=0; i<d; i++)
	bases[i] = rcp(new basis_type(p,true));
      basis = rcp(new product_basis_type(bases, 1e-12));
      stoch_length = basis->size();
      Cijk = basis->computeTripleProductTensor();
      
      // Create stochastic parallel distribution
      ParameterList parallelParams;
      RCP<Stokhos::ParallelData> sg_parallel_data =
	rcp(new Stokhos::ParallelData(basis, Cijk, globalComm, parallelParams));
      RCP<const EpetraExt::MultiComm> sg_comm = 
	sg_parallel_data->getMultiComm();
      RCP<const Epetra_Comm> app_comm = 
	sg_parallel_data->getSpatialComm();
      RCP<const Stokhos::EpetraSparse3Tensor> epetraCijk =
	sg_parallel_data->getEpetraCijk();
      RCP<const Epetra_BlockMap> stoch_row_map = 
	epetraCijk->getStochasticRowMap();

      //------------------------------

      // Generate Epetra objects
      RCP<const Epetra_Map> x_map = 
	rcp(new Epetra_Map(fem_length, 0, *app_comm));
      RCP<const Epetra_Map> sg_x_map = 
	rcp(EpetraExt::BlockUtility::GenerateBlockMap(
	      *x_map, *stoch_row_map, *sg_comm));

      RCP<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy, *x_map, 27));
      int *my_GIDs = x_map->MyGlobalElements();
      int num_my_GIDs = x_map->NumMyElements();
      for (int i=0; i<num_my_GIDs; ++i) {
	int row = my_GIDs[i];
	int num_indices = fem_graph[row].size();
	int *indices = &(fem_graph[row][0]);
	graph->InsertGlobalIndices(row, num_indices, indices);
      }
      graph->FillComplete();
      
      RCP<ParameterList> sg_op_params = rcp(new ParameterList);
      RCP<Stokhos::MatrixFreeOperator> sg_A = 
	rcp(new Stokhos::MatrixFreeOperator(sg_comm, basis, epetraCijk, 
					    x_map, x_map, sg_x_map, sg_x_map, 
					    sg_op_params));
      RCP<Epetra_BlockMap> sg_A_overlap_map =
	rcp(new Epetra_LocalMap(
	      stoch_length, 0, *(sg_parallel_data->getStochasticComm())));
      RCP< Stokhos::EpetraOperatorOrthogPoly > A_sg_blocks = 
	rcp(new Stokhos::EpetraOperatorOrthogPoly(
	      basis, sg_A_overlap_map, x_map, x_map, sg_x_map, sg_comm));
      for (int i=0; i<stoch_length; i++) {
	RCP<Epetra_CrsMatrix> A = rcp(new Epetra_CrsMatrix(Copy, *graph));

	for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < fem_length ; ++iRowFEM ) {
	  for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
	    const int iColFEM = fem_graph[iRowFEM][iRowEntryFEM] ;

	    double v = generate_matrix_coefficient<double>( 
	      fem_length , stoch_length , iRowFEM , iColFEM , i );
	    A->ReplaceGlobalValues(iRowFEM, 1, &v, &iColFEM);
	  }
	}
	A->FillComplete();
	//A->PutScalar(1.0);
	A_sg_blocks->setCoeffPtr(i, A);
      }
      sg_A->setupOperator(A_sg_blocks);

      sg_x = rcp(new Stokhos::EpetraVectorOrthogPoly(
		   basis, stoch_row_map, x_map, sg_x_map, sg_comm));
      sg_y = rcp(new Stokhos::EpetraVectorOrthogPoly(
		   basis, stoch_row_map, x_map, sg_x_map, sg_comm));

      // Initialize vectors
      for (int iColFEM=0; iColFEM < fem_length; ++iColFEM ) {
	for (int iColStoch=0 ; iColStoch < stoch_length; ++iColStoch ) {
	  (*sg_x)[iColStoch][iColFEM] = generate_vector_coefficient<double>( 
	    fem_length , stoch_length , iColFEM , iColStoch );
	}
      }
      sg_y->init(0.0);

      // Apply operator
      sg_A->Apply( *(sg_x->getBlockVector()), *(sg_y->getBlockVector()) );

      // Transpose y to commuted layout
      sg_y_commuted = rcp(new Stokhos::ProductEpetraVector(
			    x_map, stoch_row_map, sg_comm));
      for (int block=0; block<stoch_length; ++block) {
	for (int i=0; i<fem_length; ++i)
	  (*sg_y_commuted)[i][block] = (*sg_y)[block][i];
      }

      const std::vector<unsigned> var_degree( d , p );
      const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
	tensor( var_degree );
      perm.resize(stoch_length);
      inv_perm.resize(stoch_length);
      for (int i=0; i<stoch_length; ++i) {
	const unsigned char *t = tensor[i];
	Stokhos::MultiIndex<int> term(d);
	for (int j=0; j<d; ++j)
	  term[j] = t[j];
	int idx = basis->index(term);
	perm[idx] = i;
	inv_perm[i] = idx;
      }

      typedef KokkosArray::NormalizedLegendrePolynomialBases<8> polynomial ;
      typedef KokkosArray::StochasticProductTensor< value_type , polynomial , KokkosArray::Host , KokkosArray::CrsProductTensor > tensor_type ;

      tensor_type tensor2 = 
	KokkosArray::create_product_tensor< tensor_type >( var_degree );
  
      perm2.resize(stoch_length);
      inv_perm2.resize(stoch_length);
      for (int i=0; i<stoch_length; ++i) {
	Stokhos::MultiIndex<int> term(d);
	for (int j=0; j<d; ++j)
	  term[j] = tensor2.bases_degree(i,j);
	int idx = basis->index(term);
	perm2[idx] = i;
	inv_perm2[i] = idx;
      }
    }

    template <typename vec_type>
    bool test_original(const std::vector<vec_type>& y,
		       Teuchos::FancyOStream& out) const {
      bool success = true;
      for (int block=0; block<stoch_length; ++block) {
	for (int i=0; i<fem_length; ++i) {
	  double diff = std::abs( (*sg_y)[block][i] - y[block][i] );
	  double tol = rel_tol*std::abs((*sg_y)[block][i]) + abs_tol;
	  bool s = diff < tol;
	  if (!s) {
	    out << "y_expected[" << block << "][" << i << "] - "
		<< "y[" << block << "][" << i << "] == "
		<< diff << " < " << tol << " : failed"
		<< std::endl;
	  }
	  success = success && s;
	}
      }

      return success;
    }

    template <typename vec_type>
    bool test_commuted(const vec_type& y,
		       Teuchos::FancyOStream& out) const {
      bool success = true;
      for (int block=0; block<stoch_length; ++block) {
	int b = perm[block];
	for (int i=0; i<fem_length; ++i) {
	  double diff = std::abs( (*sg_y)[block][i] - y(b,i) );
	  double tol = rel_tol*std::abs((*sg_y)[block][i]) + abs_tol;
	  bool s = diff < tol;
	  if (!s) {
	    out << "y_expected[" << block << "][" << i << "] - "
		<< "y(" << b << "," << i << ") == "
		<< diff << " < " << tol << " : failed"
		<< std::endl;
	  }
	  success = success && s;
	}
      }

      return success;
    }

    template <typename vec_type>
    bool test_commuted2(const vec_type& y,
			Teuchos::FancyOStream& out) const {
      bool success = true;
      for (int block=0; block<stoch_length; ++block) {
	int b = perm2[block];
	for (int i=0; i<fem_length; ++i) {
	  double diff = std::abs( (*sg_y)[block][i] - y(b,i) );
	  double tol = rel_tol*std::abs((*sg_y)[block][i]) + abs_tol;
	  bool s = diff < tol;
	  if (!s) {
	    out << "y_expected[" << block << "][" << i << "] - "
		<< "y(" << b << "," << i << ") == "
		<< diff << " < " << tol << " : failed"
		<< std::endl;
	  }
	  success = success && s;
	}
      }

      return success;
    }

    template <typename vec_type>
    bool test_commuted_flat(const vec_type& y,
			    Teuchos::FancyOStream& out) const {
      bool success = true;
      for (int block=0; block<stoch_length; ++block) {
	int b = perm[block];
	for (int i=0; i<fem_length; ++i) {
	  double diff = std::abs( (*sg_y)[block][i] - y(i*stoch_length+b) );
	  double tol = rel_tol*std::abs((*sg_y)[block][i]) + abs_tol;
	  bool s = diff < tol;
	  if (!s) {
	    out << "y_expected[" << block << "][" << i << "] - "
		<< "y(" << b << "," << i << ") == "
		<< diff << " < " << tol << " : failed"
		<< std::endl;
	  }
	  success = success && s;
	}
      }

      return success;
    }

    template <typename vec_type>
    bool test_original_flat(const vec_type& y,
			    Teuchos::FancyOStream& out) const {
      bool success = true;
      for (int block=0; block<stoch_length; ++block) {
	int b = perm[block];
	for (int i=0; i<fem_length; ++i) {
	  double diff = std::abs( (*sg_y)[block][i] - y(b*fem_length+i) );
	  double tol = rel_tol*std::abs((*sg_y)[block][i]) + abs_tol;
	  bool s = diff < tol;
	  if (!s) {
	    out << "y_expected[" << block << "][" << i << "] - "
		<< "y(" << b << "," << i << ") == "
		<< diff << " < " << tol << " : failed"
		<< std::endl;
	  }
	  success = success && s;
	}
      }

      return success;
    }

    template <typename vec_type>
    bool test_commuted_x(const vec_type& x,
			 Teuchos::FancyOStream& out) const {
      bool success = true;
      for (int block=0; block<stoch_length; ++block) {
	int b = perm[block];
	for (int i=0; i<fem_length; ++i) {
	  double diff = std::abs( (*sg_x)[block][i] - x(b,i) );
	  double tol = rel_tol*std::abs((*sg_x)[block][i]) + abs_tol;
	  bool s = diff < tol;
	  if (!s) {
	    out << "x_expected[" << block << "][" << i << "] - "
		<< "x(" << b << "," << i << ") == "
		<< diff << " < " << tol << " : failed"
		<< std::endl;
	  }
	  success = success && s;
	}
      }

      return success;
    }
    
  };

  template <typename value_type, typename Device>
  bool test_crs_matrix_free(const UnitTestSetup& setup,
			    bool test_block,
			    Teuchos::FancyOStream& out) {
    typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
    typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;
    typedef KokkosArray::View<value_type[],Device> vec_type ;

    //------------------------------
    
    std::vector<matrix_type> matrix( setup.stoch_length ) ;
    std::vector<vec_type> x( setup.stoch_length ) ;
    std::vector<vec_type> y( setup.stoch_length ) ;
    std::vector<vec_type> tmp( setup.stoch_length ) ;

    for (int block=0; block<setup.stoch_length; ++block) {
      matrix[block].graph = KokkosArray::create_crsarray<crsarray_type>( 
	std::string("testing") , setup.fem_graph );

      matrix[block].values = vec_type( "matrix" , setup.fem_graph_length );

      x[block]   = vec_type( "x" , setup.fem_length );
      y[block]   = vec_type( "y" , setup.fem_length );
      tmp[block] = vec_type( "tmp" , setup.fem_length );

      typename vec_type::HostMirror hM =
	KokkosArray::create_mirror( matrix[block].values );

      for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
	for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
	  const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

	  hM(iEntryFEM) = generate_matrix_coefficient<value_type>( 
	    setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , block );
	  //hM(iEntryFEM) = 1.0;
	}
      }
      
      KokkosArray::deep_copy( matrix[block].values , hM );

      typename vec_type::HostMirror hx =
	KokkosArray::create_mirror( x[block] );
      
      for ( int i = 0 ; i < setup.fem_length ; ++i ) {
	hx(i) = generate_vector_coefficient<value_type>( 
	  setup.fem_length , setup.stoch_length , i , block ); ;
      }
      
      KokkosArray::deep_copy( x[block] , hx );
    }

    // Original matrix-free multiply algorithm using a block apply
    typename UnitTestSetup::Cijk_type::k_iterator k_begin = 
      setup.Cijk->k_begin();
    typename UnitTestSetup::Cijk_type::k_iterator k_end = 
      setup.Cijk->k_end();
    for (typename UnitTestSetup::Cijk_type::k_iterator k_it=k_begin; 
	 k_it!=k_end; ++k_it) {
      int nj = setup.Cijk->num_j(k_it);
      if (nj > 0) {
	int k = index(k_it);
	typename UnitTestSetup::Cijk_type::kj_iterator j_begin = 
	  setup.Cijk->j_begin(k_it);
	typename UnitTestSetup::Cijk_type::kj_iterator j_end = 
	  setup.Cijk->j_end(k_it);
	std::vector<vec_type> xx(nj), yy(nj);
	int jdx = 0;
	for (typename UnitTestSetup::Cijk_type::kj_iterator j_it = j_begin; 
	     j_it != j_end; 
	     ++j_it) {
	  int j = index(j_it);
	  xx[jdx] = x[j];
	  yy[jdx] = tmp[j];
	  jdx++;
	}
        KokkosArray::multiply( matrix[k] , xx , yy, test_block );
	jdx = 0;
	for (typename UnitTestSetup::Cijk_type::kj_iterator j_it = j_begin; 
	     j_it != j_end; ++j_it) {
	  typename UnitTestSetup::Cijk_type::kji_iterator i_begin = 
	    setup.Cijk->i_begin(j_it);
	  typename UnitTestSetup::Cijk_type::kji_iterator i_end =  
	    setup.Cijk->i_end(j_it);
	  for (typename UnitTestSetup::Cijk_type::kji_iterator i_it = i_begin; 
	       i_it != i_end; 
	       ++i_it) {
	    int i = index(i_it);
	    value_type c = value(i_it);
	    KokkosArray::update( value_type(1.0) , y[i] , c , yy[jdx] );
	  }
	  jdx++;
	}
      }
    }

    std::vector<typename vec_type::HostMirror> hy(setup.stoch_length);
    for (int i=0; i<setup.stoch_length; ++i) {
      hy[i] = KokkosArray::create_mirror( y[i] );
      KokkosArray::deep_copy( hy[i] , y[i] );
    }
    bool success = setup.test_original(hy, out);

    return success;
  }

  template <typename Scalar, typename Device>
  bool test_crs_product_legendre(const UnitTestSetup& setup,
				 Teuchos::FancyOStream& out) {
    typedef Scalar VectorScalar;
    typedef Scalar MatrixScalar;
    typedef Scalar TensorScalar;

    typedef KokkosArray::View< VectorScalar** , KokkosArray::LayoutLeft ,
      Device > vector_type ;

    typedef KokkosArray::CrsProductTensorLegendre< TensorScalar , Device >  tensor_type ;

    typedef KokkosArray::BlockCrsMatrix< tensor_type , MatrixScalar , Device > matrix_type ;

    typedef typename matrix_type::graph_type graph_type ;

    //------------------------------
    // Generate CRS block-tensor matrix:

    const std::vector<unsigned> var_degree( setup.d , setup.p );

    const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
      tensor( var_degree );

    TEUCHOS_TEST_FOR_EXCEPTION(
      setup.stoch_length != static_cast<int>(tensor.bases_count()), 
      std::logic_error,
      "tensor.bases_count() == " << tensor.bases_count() << " != " << 
      "setup.stoch_length == " << setup.stoch_length << std::endl);

    std::vector< std::vector< size_t > > stoch_graph( setup.stoch_length );
    
    for ( int i = 0 ; i < setup.stoch_length ; ++i ) {
      for ( int j = 0 ; j < setup.stoch_length ; ++j ) {
	if ( KokkosArray::matrix_nonzero(tensor,i,j) ) {
	  stoch_graph[i].push_back(j);
	}
      }
    }

    //------------------------------
    // Generate input multivector:
    
    vector_type x = vector_type( "x" , setup.stoch_length , setup.fem_length );
    vector_type y = vector_type( "y" , setup.stoch_length , setup.fem_length );

    typename vector_type::HostMirror hx = 
      KokkosArray::create_mirror( x );

    for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
      for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
	hx(setup.perm[iColStoch],iColFEM) = 
	  generate_vector_coefficient<Scalar>( 
	    setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
      }
    }

    KokkosArray::deep_copy( x , hx );

    //------------------------------

    matrix_type matrix ;

    matrix.block = tensor_type( var_degree );

    matrix.graph = KokkosArray::create_crsarray<graph_type>( 
      std::string("test crs graph") , setup.fem_graph );

    TEUCHOS_TEST_FOR_EXCEPTION(
      setup.stoch_length != static_cast<int>(matrix.block.dimension()), 
      std::logic_error,
      "matrix.block.dimension() == " << matrix.block.dimension() << " != " << 
      "setup.stoch_length == " << setup.stoch_length << std::endl);

    matrix.values = vector_type( "matrix" , setup.stoch_length , setup.fem_graph_length );

    typename vector_type::HostMirror hM = 
      KokkosArray::create_mirror( matrix.values );

    for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
      for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
	const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

	for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
	  hM(setup.perm[k],iEntryFEM) = generate_matrix_coefficient<Scalar>( 
	    setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , k );
	  //hM(k,iEntryFEM) = 1.0;
	}
      }
    }

    KokkosArray::deep_copy( matrix.values , hM );

    //------------------------------

    const KokkosArray::Impl::Multiply< matrix_type , vector_type , vector_type > op( matrix , x , y );
    op.run();

    typename vector_type::HostMirror hy = KokkosArray::create_mirror( y );
    KokkosArray::deep_copy( hy , y );

    bool success = setup.test_commuted(hy, out);
    return success;
  }

  template< typename ScalarType , class Device >
  bool
  test_crs_dense_block(const UnitTestSetup& setup,
		       Teuchos::FancyOStream& out)
  {
    typedef ScalarType value_type ;
    typedef KokkosArray::View< value_type**, KokkosArray::LayoutLeft , Device > block_vector_type ;
    
    //------------------------------

    typedef KokkosArray::BlockCrsMatrix< KokkosArray::SymmetricDiagonalSpec< Device > , value_type , Device > matrix_type ;
    
    typedef typename matrix_type::graph_type  graph_type ;
    
    //------------------------------
    // Generate product tensor from variables' degrees
    
    const std::vector<unsigned> var_degree( setup.d , setup.p );
    
    const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
      tensor( var_degree );
    
    TEUCHOS_TEST_FOR_EXCEPTION(
      setup.stoch_length != static_cast<int>(tensor.bases_count()), 
      std::logic_error,
      "tensor.bases_count() == " << tensor.bases_count() << " != " << 
      "setup.stoch_length == " << setup.stoch_length << std::endl);
    
    
    //------------------------------
    
    block_vector_type x = block_vector_type( "x" , setup.stoch_length , setup.fem_length );
    block_vector_type y = block_vector_type( "y" , setup.stoch_length , setup.fem_length );
    
    typename block_vector_type::HostMirror hx        = KokkosArray::create_mirror( x );

    for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
      for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
	hx(setup.perm[iColStoch],iColFEM) =
	  generate_vector_coefficient<ScalarType>( 
	    setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
      }
    }

    KokkosArray::deep_copy( x , hx );

    //------------------------------
    // Generate CRS matrix of blocks with symmetric diagonal storage

    matrix_type matrix ;

    matrix.block  = KokkosArray::SymmetricDiagonalSpec< Device >( setup.stoch_length );
    matrix.graph  = KokkosArray::create_crsarray<graph_type>( std::string("test product tensor graph") , setup.fem_graph );
    matrix.values = block_vector_type( "matrix" , matrix.block.matrix_size() , setup.fem_graph_length );

    {
      typename block_vector_type::HostMirror hM =
	KokkosArray::create_mirror( matrix.values );

      for ( int iRowStoch = 0 ; iRowStoch < setup.stoch_length ; ++iRowStoch ) {
	for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {

	  for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
	    const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM];

	    for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
	      
	      const size_t offset = 
		matrix.block.matrix_offset( setup.perm[iRowStoch] , 
					    setup.perm[iColStoch] );
	      
	      ScalarType value = 0 ;
	      
	      for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
		value += tensor( setup.perm[iRowStoch] , 
				 setup.perm[iColStoch] , 
				 setup.perm[k] ) *
		  generate_matrix_coefficient<ScalarType>( 
		    setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , k );
	      }

	      hM( offset , iEntryFEM ) = value ;
	    }
	  }

	}
      }

      KokkosArray::deep_copy( matrix.values , hM );
    }

    //------------------------------

    KokkosArray::multiply( matrix , x , y );

    typename block_vector_type::HostMirror hy = KokkosArray::create_mirror( y );
    KokkosArray::deep_copy( hy , y );

    bool success = setup.test_commuted(hy, out);

    return success;
  }

  template< typename ScalarType , class Device >
  bool
  test_crs_flat_commuted(const UnitTestSetup& setup,
			 Teuchos::FancyOStream& out)
  {
    typedef ScalarType value_type ;
    typedef KokkosArray::View< value_type[] , Device > vector_type ;

    //------------------------------

    typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
    typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

    //------------------------------
    // Tensor evaluation:

    const std::vector<unsigned> var_degree( setup.d , setup.p );

    const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
      tensor( var_degree );

    TEUCHOS_TEST_FOR_EXCEPTION(
      setup.stoch_length != static_cast<int>(tensor.bases_count()), 
      std::logic_error,
      "tensor.bases_count() == " << tensor.bases_count() << " != " << 
      "setup.stoch_length == " << setup.stoch_length << std::endl);

    std::vector< std::vector< int > > stoch_graph( setup.stoch_length );

    for ( int i = 0 ; i < setup.stoch_length ; ++i ) {
      for ( int j = 0 ; j < setup.stoch_length ; ++j ) {
	if ( KokkosArray::matrix_nonzero(tensor,i,j) ) {
	  stoch_graph[i].push_back(j);
	}
      }
    }

    //------------------------------
    // Generate flattened graph with FEM outer and stochastic inner

    const int flat_length = setup.fem_length * setup.stoch_length ;

    std::vector< std::vector<int> > flat_graph( flat_length );

    for ( int iOuterRow = 0 ; iOuterRow < setup.fem_length ; ++iOuterRow ) {

      const size_t iOuterRowNZ = setup.fem_graph[iOuterRow].size();

      for ( int iInnerRow = 0 ; iInnerRow < setup.stoch_length ; ++iInnerRow ) {

	const size_t iInnerRowNZ = stoch_graph[ iInnerRow ].size(); ;
	const int iFlatRowNZ  = iOuterRowNZ * iInnerRowNZ ;
	const int iFlatRow    = iInnerRow + iOuterRow * setup.stoch_length ;

	flat_graph[iFlatRow].resize( iFlatRowNZ );

	int iFlatEntry = 0 ;

	for ( size_t iOuterEntry = 0 ; iOuterEntry < iOuterRowNZ ; ++iOuterEntry ) {

	  const int iOuterCol = setup.fem_graph[iOuterRow][iOuterEntry];

	  for ( size_t iInnerEntry = 0 ; iInnerEntry < iInnerRowNZ ; ++iInnerEntry ) {

	    const int iInnerCol   = stoch_graph[iInnerRow][iInnerEntry] ;
	    const int iFlatColumn = iInnerCol + iOuterCol * setup.stoch_length ;

	    flat_graph[iFlatRow][iFlatEntry] = iFlatColumn ;

	    ++iFlatEntry ;
	  }
	}
      }
    }

    //------------------------------

    vector_type x = vector_type( "x" , flat_length );
    vector_type y = vector_type( "y" , flat_length );

    typename vector_type::HostMirror hx = KokkosArray::create_mirror( x );

    for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
      for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
	hx(setup.perm[iColStoch] + iColFEM*setup.stoch_length) =
	  generate_vector_coefficient<ScalarType>( 
	    setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
      }
    }

    KokkosArray::deep_copy( x , hx );

    //------------------------------

    matrix_type matrix ;

    matrix.graph = KokkosArray::create_crsarray<crsarray_type>( 
      std::string("testing") , flat_graph );

    const size_t flat_graph_length = matrix.graph.entries.dimension_0();

    matrix.values = vector_type( "matrix" , flat_graph_length );
    {
      typename vector_type::HostMirror hM =
	KokkosArray::create_mirror( matrix.values );

      for ( int iRow = 0 , iEntry = 0 ; iRow < flat_length ; ++iRow ) {
	const int iRowFEM   = iRow / setup.stoch_length ;
	const int iRowStoch = iRow % setup.stoch_length ;

	for ( size_t iRowEntry = 0 ; iRowEntry < flat_graph[ iRow ].size() ; ++iRowEntry , ++iEntry ) {
	  const int iCol = flat_graph[ iRow ][ iRowEntry ];
	  const int iColFEM   = iCol / setup.stoch_length ;
	  const int iColStoch = iCol % setup.stoch_length ;

	  double value = 0 ;
	  for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
	    const double A_fem_k = 
	      generate_matrix_coefficient<ScalarType>( 
		setup.fem_length , setup.stoch_length , iRowFEM, iColFEM, 
		setup.inv_perm[k] );
	    value += tensor(iRowStoch,iColStoch,k) * A_fem_k ;
	  }
	  hM( iEntry ) = value ;
	}
      }

      KokkosArray::deep_copy( matrix.values , hM );
    }

    //------------------------------

    
    KokkosArray::multiply( matrix , x , y );

    typename vector_type::HostMirror hy = KokkosArray::create_mirror( y );
    KokkosArray::deep_copy( hy , y );

    bool success = setup.test_commuted_flat(hy, out);
    return success;
  }

  template< typename ScalarType , class Device >
  bool
  test_crs_flat_original(const UnitTestSetup& setup,
			 Teuchos::FancyOStream& out)
  {
    typedef ScalarType value_type ;
    typedef KokkosArray::View< value_type[] , Device > vector_type ;

    //------------------------------

    typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
    typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

    //------------------------------
    // Tensor evaluation:

    const std::vector<unsigned> var_degree( setup.d , setup.p );

    const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
      tensor( var_degree );

    TEUCHOS_TEST_FOR_EXCEPTION(
      setup.stoch_length != static_cast<int>(tensor.bases_count()), 
      std::logic_error,
      "tensor.bases_count() == " << tensor.bases_count() << " != " << 
      "setup.stoch_length == " << setup.stoch_length << std::endl);

    std::vector< std::vector< int > > stoch_graph( setup.stoch_length );

    for ( int i = 0 ; i < setup.stoch_length ; ++i ) {
      for ( int j = 0 ; j < setup.stoch_length ; ++j ) {
	if ( KokkosArray::matrix_nonzero(tensor,i,j) ) {
	  stoch_graph[i].push_back(j);
	}
      }
    }

    //------------------------------
    // Generate flattened graph with stochastic outer and FEM inner

    const size_t flat_length = setup.fem_length * setup.stoch_length ;

    std::vector< std::vector<int> > flat_graph( flat_length );
  
    for ( int iOuterRow = 0 ; iOuterRow < setup.stoch_length ; ++iOuterRow ) {

      const size_t iOuterRowNZ = stoch_graph[iOuterRow].size();

      for ( int iInnerRow = 0 ; iInnerRow < setup.fem_length ; ++iInnerRow ) {

	const size_t iInnerRowNZ = setup.fem_graph[iInnerRow].size();
	const int iFlatRowNZ  = iOuterRowNZ * iInnerRowNZ ;
	const int iFlatRow    = iInnerRow + iOuterRow * setup.fem_length ;

	flat_graph[iFlatRow].resize( iFlatRowNZ );

	int iFlatEntry = 0 ;

	for ( size_t iOuterEntry = 0 ; iOuterEntry < iOuterRowNZ ; ++iOuterEntry ) {

	  const int iOuterCol = stoch_graph[ iOuterRow ][ iOuterEntry ];

	  for ( size_t iInnerEntry = 0 ; iInnerEntry < iInnerRowNZ ; ++iInnerEntry ) {

	    const int iInnerCol   = setup.fem_graph[ iInnerRow][iInnerEntry];
	    const int iFlatColumn = iInnerCol + iOuterCol * setup.fem_length ;
	    
	    flat_graph[iFlatRow][iFlatEntry] = iFlatColumn ;
	    ++iFlatEntry ;
	  }
	}
      }
    }

    //------------------------------

    vector_type x = vector_type( "x" , flat_length );
    vector_type y = vector_type( "y" , flat_length );

    typename vector_type::HostMirror hx = KokkosArray::create_mirror( x );

    for ( size_t iCol = 0 ; iCol < flat_length ; ++iCol ) {
      const int iColStoch = iCol / setup.fem_length ;
      const int iColFEM   = iCol % setup.fem_length ;

      hx(iCol) = generate_vector_coefficient<ScalarType>( 
	setup.fem_length , setup.stoch_length ,
	iColFEM , setup.inv_perm[iColStoch] );
    }

    KokkosArray::deep_copy( x , hx );

    //------------------------------

    matrix_type matrix ;

    matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , flat_graph );

    const size_t flat_graph_length = matrix.graph.entries.dimension_0();

    matrix.values = vector_type( "matrix" , flat_graph_length );
    {
      typename vector_type::HostMirror hM =
	KokkosArray::create_mirror( matrix.values );

      for ( size_t iRow = 0 , iEntry = 0 ; iRow < flat_length ; ++iRow ) {
	const int iRowStoch = iRow / setup.fem_length ;
	const int iRowFEM   = iRow % setup.fem_length ;
	
	for ( size_t iRowEntry = 0 ; iRowEntry < flat_graph[ iRow ].size() ; ++iRowEntry , ++iEntry ) {
	  const int iCol = flat_graph[ iRow ][ iRowEntry ];
	  const int iColStoch = iCol / setup.fem_length ;
	  const int iColFEM   = iCol % setup.fem_length ;
	  
	  double value = 0 ;
	  for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
	    const double A_fem_k =
	      generate_matrix_coefficient<ScalarType>( 
		setup.fem_length , setup.stoch_length ,
		iRowFEM , iColFEM , setup.inv_perm[k] );
	    value += tensor(iRowStoch,iColStoch,k) * A_fem_k ;
	  }
	  hM( iEntry ) = value ;
	  
	}
	
      }
      
      KokkosArray::deep_copy( matrix.values , hM );
    }
    
    KokkosArray::multiply( matrix , x , y );

    typename vector_type::HostMirror hy = KokkosArray::create_mirror( y );
    KokkosArray::deep_copy( hy , y );

    bool success = setup.test_original_flat(hy, out);
    return success;
  }

  template< typename ScalarType , class Device ,
	    template< unsigned , typename , class > class TensorType >
  bool test_crs_product_tensor(const UnitTestSetup& setup,
			       Teuchos::FancyOStream& out) {
    typedef ScalarType value_type ;
    typedef KokkosArray::View< value_type** , KokkosArray::LayoutLeft ,
                               Device > block_vector_type ;

    typedef KokkosArray::NormalizedLegendrePolynomialBases<8> polynomial ;

    typedef KokkosArray::StochasticProductTensor< value_type , polynomial , Device , TensorType > tensor_type ;

    typedef KokkosArray::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
    typedef typename matrix_type::graph_type graph_type ;

    const std::vector<unsigned> var_degree( setup.d , setup.p );

    //------------------------------
    // Generate input multivector:

    block_vector_type x = 
      block_vector_type( "x" , setup.stoch_length , setup.fem_length );
    block_vector_type y = 
      block_vector_type( "y" , setup.stoch_length , setup.fem_length );

    typename block_vector_type::HostMirror hx = 
      KokkosArray::create_mirror( x );

    for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
      for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
	hx(setup.perm2[iColStoch],iColFEM) = 
	  generate_vector_coefficient<ScalarType>( 
	    setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
      }
    }

    KokkosArray::deep_copy( x , hx );

    //------------------------------

    matrix_type matrix ;

    matrix.block = 
      KokkosArray::create_product_tensor< tensor_type >( var_degree );

    matrix.graph = KokkosArray::create_crsarray<graph_type>( 
      std::string("test crs graph") , setup.fem_graph );

    matrix.values = block_vector_type( 
      "matrix" , setup.stoch_length , setup.fem_graph_length );

    typename block_vector_type::HostMirror hM = 
      KokkosArray::create_mirror( matrix.values );

    for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
      for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
	const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

	for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
	  hM(setup.perm2[k],iEntryFEM) = 
	    generate_matrix_coefficient<ScalarType>(
	      setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , k );
	  //hM(k,iEntryFEM) = 1.0;
	}
      }
    }
    
    KokkosArray::deep_copy( matrix.values , hM );
    
    //------------------------------
    
    KokkosArray::multiply( matrix , x , y );

    typename block_vector_type::HostMirror hy = KokkosArray::create_mirror( y );
    KokkosArray::deep_copy( hy , y );

    bool success = setup.test_commuted2(hy, out);
    //bool success = true;
    return success;
  }

}

#endif
