// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_KOKKOS_CORE_KERNELS_UNIT_TEST_HPP
#define STOKHOS_KOKKOS_CORE_KERNELS_UNIT_TEST_HPP

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestHelpers.hpp"
#include "Teuchos_TestForException.hpp"

#include "Stokhos_Epetra.hpp"
#include "Stokhos_Sparse3TensorUtilities.hpp"
#include "EpetraExt_BlockUtility.h"
#include "Stokhos_UnitTestHelpers.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Kokkos_Macros.hpp"

#include "Stokhos_Update.hpp"
#include "Stokhos_CrsMatrix.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_StochasticProductTensor.hpp"
#include "Stokhos_CrsProductTensor.hpp"
#include "Stokhos_TiledCrsProductTensor.hpp"
#include "Stokhos_SimpleTiledCrsProductTensor.hpp"
#include "Stokhos_CooProductTensor.hpp"
#include "Stokhos_SymmetricDiagonalSpec.hpp"
#include "Stokhos_FlatSparse3Tensor.hpp"
#include "Stokhos_FlatSparse3Tensor_kji.hpp"
#include "Stokhos_LinearSparse3Tensor.hpp"
#include "Stokhos_LexicographicBlockSparse3Tensor.hpp"

#ifdef HAVE_STOKHOS_KOKKOSLINALG
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_update.hpp"
#endif

namespace KokkosKernelsUnitTest {

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
  //return 1.0;
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
  //return 1.0;
}

template <typename Device>
struct UnitTestSetup {
  typedef double value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::LegendreBasis<int,value_type> basis_type;
  //typedef Stokhos::CompletePolynomialBasis<int,value_type> product_basis_type;
  typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
  typedef Stokhos::TotalOrderBasis<int,value_type,less_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  int p, d, nGrid, fem_length, stoch_length, stoch_length_aligned, fem_graph_length;
  double rel_tol, abs_tol;
  std::vector< std::vector<int> > fem_graph ;
  RCP< product_basis_type> basis;
  RCP<Cijk_type> Cijk;
  RCP<const Epetra_Comm> globalComm;
  RCP<Stokhos::EpetraVectorOrthogPoly> sg_x, sg_y;
  RCP<Stokhos::ProductEpetraVector> sg_y_commuted;
  Teuchos::Array<int> perm, inv_perm;

  // Can't be a constructor because MPI will not be initialized
  void setup(int p_ = 5, int d_ = 2) {

    p = p_;
    d = d_;
    nGrid = 3;
    rel_tol = 1e-12;
    abs_tol = 1e-12;

    // Create a communicator for Epetra objects
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

    typedef typename Kokkos::ViewTraits<value_type,Device,void,void>::host_mirror_space host_device;
    typedef Stokhos::CrsProductTensor< value_type , host_device > tensor_type;
    typedef Stokhos::StochasticProductTensor< value_type , tensor_type ,
      host_device > stoch_tensor_type ;

    stoch_tensor_type tensor =
      Stokhos::create_stochastic_product_tensor< tensor_type >( *basis, *Cijk );
    stoch_length_aligned = tensor.aligned_dimension();

    perm.resize(stoch_length);
    inv_perm.resize(stoch_length);
    for (int i=0; i<stoch_length; ++i) {
      Stokhos::MultiIndex<int> term(d);
      for (int j=0; j<d; ++j)
        term[j] = tensor.bases_degree(i,j);
      int idx = basis->index(term);
      perm[idx] = i;
      inv_perm[i] = idx;
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
              << "y[" << block << "][" << i << "] = " << (*sg_y)[block][i]
              << " - " << y[block][i] << " == "
              << diff << " < " << tol << " : failed"
              << std::endl;
        }
        success = success && s;
      }
    }

    return success;
  }

  template <typename multi_vec_type>
  bool test_original(const multi_vec_type& y,
                     Teuchos::FancyOStream& out) const {
    bool success = true;
    for (int block=0; block<stoch_length; ++block) {
      for (int i=0; i<fem_length; ++i) {
        double diff = std::abs( (*sg_y)[block][i] - y(i,block) );
        double tol = rel_tol*std::abs((*sg_y)[block][i]) + abs_tol;
        bool s = diff < tol;
        if (!s) {
          out << "y_expected[" << block << "][" << i << "] - "
              << "y(" << i << "," << block << ") = " << (*sg_y)[block][i]
              << " - " << y(i,block) << " == "
              << diff << " < " << tol << " : failed"
              << std::endl;
        }
        success = success && s;
      }
    }

    return success;
  }

  template <typename vec_type>
  bool test_commuted_perm(const vec_type& y,
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
              << "y(" << b << "," << i << ") = " << (*sg_y)[block][i]
              << " - " << y(b,i) << " == "
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
      int b = block;
      for (int i=0; i<fem_length; ++i) {
        double diff = std::abs( (*sg_y)[block][i] - y(b,i) );
        double tol = rel_tol*std::abs((*sg_y)[block][i]) + abs_tol;
        bool s = diff < tol;
        if (!s) {
          out << "y_expected[" << block << "][" << i << "] - "
              << "y(" << b << "," << i << ") = " << (*sg_y)[block][i] << " - "
              << y(b,i) << " == "
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
      int b = block;
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
      int b = block;
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

};

template <typename value_type, typename Device, typename SparseMatOps>
bool test_crs_matrix_free(const UnitTestSetup<Device>& setup,
                          Teuchos::FancyOStream& out) {
  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::graph_type matrix_graph_type;
  typedef Kokkos::View<value_type*,Device> vec_type;

  //------------------------------

  std::vector<matrix_type> matrix( setup.stoch_length ) ;
  std::vector<vec_type> x( setup.stoch_length ) ;
  std::vector<vec_type> y( setup.stoch_length ) ;
  std::vector<vec_type> tmp( setup.stoch_length ) ;

  for (int block=0; block<setup.stoch_length; ++block) {
    matrix[block].graph = Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("testing") , setup.fem_graph );

    matrix[block].values =
      matrix_values_type( "matrix" , setup.fem_graph_length );

    x[block]   = vec_type( "x" , setup.fem_length );
    y[block]   = vec_type( "y" , setup.fem_length );
    tmp[block] = vec_type( "tmp" , setup.fem_length );

    typename matrix_values_type::HostMirror hM =
      Kokkos::create_mirror( matrix[block].values );

    for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
      for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
        const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

        hM(iEntryFEM) = generate_matrix_coefficient<value_type>(
          setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , block );
      }
    }

    Kokkos::deep_copy( matrix[block].values , hM );

    typename vec_type::HostMirror hx =
      Kokkos::create_mirror( x[block] );
    typename vec_type::HostMirror hy =
      Kokkos::create_mirror( y[block] );

    for ( int i = 0 ; i < setup.fem_length ; ++i ) {
      hx(i) = generate_vector_coefficient<value_type>(
        setup.fem_length , setup.stoch_length , i , block );
      hy(i) = 0.0;
    }

    Kokkos::deep_copy( x[block] , hx );
    Kokkos::deep_copy( y[block] , hy );
  }

  // Original matrix-free multiply algorithm using a block apply
  SparseMatOps smo;
  typename UnitTestSetup<Device>::Cijk_type::k_iterator k_begin =
    setup.Cijk->k_begin();
  typename UnitTestSetup<Device>::Cijk_type::k_iterator k_end =
    setup.Cijk->k_end();
  for (typename UnitTestSetup<Device>::Cijk_type::k_iterator k_it=k_begin;
       k_it!=k_end; ++k_it) {
    int nj = setup.Cijk->num_j(k_it);
    if (nj > 0) {
      int k = index(k_it);
      typename UnitTestSetup<Device>::Cijk_type::kj_iterator j_begin =
        setup.Cijk->j_begin(k_it);
      typename UnitTestSetup<Device>::Cijk_type::kj_iterator j_end =
        setup.Cijk->j_end(k_it);
      std::vector<vec_type> xx(nj), yy(nj);
      int jdx = 0;
      for (typename UnitTestSetup<Device>::Cijk_type::kj_iterator j_it = j_begin;
           j_it != j_end;
           ++j_it) {
        int j = index(j_it);
        xx[jdx] = x[j];
        yy[jdx] = tmp[j];
        jdx++;
      }
      Stokhos::multiply( matrix[k] , xx , yy, smo );
      jdx = 0;
      for (typename UnitTestSetup<Device>::Cijk_type::kj_iterator j_it = j_begin;
           j_it != j_end; ++j_it) {
        typename UnitTestSetup<Device>::Cijk_type::kji_iterator i_begin =
          setup.Cijk->i_begin(j_it);
        typename UnitTestSetup<Device>::Cijk_type::kji_iterator i_end =
          setup.Cijk->i_end(j_it);
        for (typename UnitTestSetup<Device>::Cijk_type::kji_iterator i_it = i_begin;
             i_it != i_end;
             ++i_it) {
          int i = index(i_it);
          value_type c = value(i_it);
          Stokhos::update( value_type(1.0) , y[i] , c , yy[jdx] );
        }
        jdx++;
      }
    }
  }

  std::vector<typename vec_type::HostMirror> hy(setup.stoch_length);
  for (int i=0; i<setup.stoch_length; ++i) {
    hy[i] = Kokkos::create_mirror( y[i] );
    Kokkos::deep_copy( hy[i] , y[i] );
  }
  bool success = setup.test_original(hy, out);

  return success;
}

template <typename value_type, typename Device, typename SparseMatOps>
bool test_crs_matrix_free_view(const UnitTestSetup<Device>& setup,
                               Teuchos::FancyOStream& out) {
  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::graph_type matrix_graph_type;
  typedef Kokkos::View<value_type*, Kokkos::LayoutLeft, Device, Kokkos::MemoryUnmanaged> vec_type;
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, Device> multi_vec_type;

  //------------------------------

  std::vector<matrix_type> matrix( setup.stoch_length ) ;
  multi_vec_type x( "x", setup.fem_length, setup.stoch_length  ) ;
  multi_vec_type y( "y", setup.fem_length, setup.stoch_length ) ;
  multi_vec_type tmp_x( "tmp_x", setup.fem_length, setup.stoch_length ) ;
  multi_vec_type tmp_y( "tmp_y", setup.fem_length, setup.stoch_length ) ;

  typename multi_vec_type::HostMirror hx = Kokkos::create_mirror( x );
  typename multi_vec_type::HostMirror hy = Kokkos::create_mirror( y );

  for (int block=0; block<setup.stoch_length; ++block) {
    matrix[block].graph = Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("testing") , setup.fem_graph );

    matrix[block].values =
      matrix_values_type( "matrix" , setup.fem_graph_length );

    typename matrix_values_type::HostMirror hM =
      Kokkos::create_mirror( matrix[block].values );

    for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
      for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
        const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

        hM(iEntryFEM) = generate_matrix_coefficient<value_type>(
          setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , block );
      }
    }

    Kokkos::deep_copy( matrix[block].values , hM );

    for ( int i = 0 ; i < setup.fem_length ; ++i ) {
      hx(i, block) = generate_vector_coefficient<value_type>(
        setup.fem_length , setup.stoch_length , i , block );
      hy(i, block) = 0.0;
    }

  }

  Kokkos::deep_copy( x , hx );
  Kokkos::deep_copy( y , hy );

  // Original matrix-free multiply algorithm using a block apply
  typedef typename UnitTestSetup<Device>::Cijk_type::k_iterator k_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kj_iterator kj_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kji_iterator kji_iterator;
  SparseMatOps smo;
  k_iterator k_begin = setup.Cijk->k_begin();
  k_iterator k_end = setup.Cijk->k_end();
  for (k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
    unsigned nj = setup.Cijk->num_j(k_it);
    if (nj > 0) {
      int k = index(k_it);
      kj_iterator j_begin = setup.Cijk->j_begin(k_it);
      kj_iterator j_end = setup.Cijk->j_end(k_it);
      std::vector<int> j_indices(nj);
      unsigned jdx = 0;
      for (kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
        int j = index(j_it);
        vec_type xx = Kokkos::subview( x, Kokkos::ALL(), j );
        vec_type tt = Kokkos::subview( tmp_x, Kokkos::ALL(), jdx++ );
        Kokkos::deep_copy(tt, xx);
      }
      multi_vec_type tmp_x_view =
        Kokkos::subview( tmp_x, Kokkos::ALL(),
                                         std::make_pair(0u,nj));
      multi_vec_type tmp_y_view =
        Kokkos::subview( tmp_y, Kokkos::ALL(),
                                         std::make_pair(0u,nj));
      Stokhos::multiply( matrix[k] , tmp_x_view , tmp_y_view, smo );
      jdx = 0;
      for (kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
        vec_type tmp_y_view =
          Kokkos::subview( tmp_y, Kokkos::ALL(), jdx++ );
        kji_iterator i_begin = setup.Cijk->i_begin(j_it);
        kji_iterator i_end = setup.Cijk->i_end(j_it);
        for (kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
          int i = index(i_it);
          value_type c = value(i_it);
          vec_type y_view = Kokkos::subview( y, Kokkos::ALL(), i );
          Stokhos::update( value_type(1.0) , y_view , c , tmp_y_view );
        }
      }
    }
  }

  Kokkos::deep_copy( hy , y );
  bool success = setup.test_original(hy, out);

  return success;
}

#ifdef HAVE_STOKHOS_KOKKOSLINALG

template <typename value_type, typename Device>
bool test_crs_matrix_free_kokkos(const UnitTestSetup<Device>& setup,
                                 Teuchos::FancyOStream& out) {
  typedef int ordinal_type;
  typedef KokkosSparse::CrsMatrix<value_type,ordinal_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef Kokkos::View<value_type*, Kokkos::LayoutLeft, Device, Kokkos::MemoryUnmanaged> vec_type;
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, Device> multi_vec_type;

  //------------------------------

  std::vector<matrix_type> matrix( setup.stoch_length ) ;
  multi_vec_type x( "x", setup.fem_length, setup.stoch_length  ) ;
  multi_vec_type y( "y", setup.fem_length, setup.stoch_length ) ;
  multi_vec_type tmp_x( "tmp_x", setup.fem_length, setup.stoch_length ) ;
  multi_vec_type tmp_y( "tmp_y", setup.fem_length, setup.stoch_length ) ;

  typename multi_vec_type::HostMirror hx = Kokkos::create_mirror( x );
  typename multi_vec_type::HostMirror hy = Kokkos::create_mirror( y );

  for (int block=0; block<setup.stoch_length; ++block) {
    matrix_graph_type matrix_graph =
      Kokkos::create_staticcrsgraph<matrix_graph_type>(
        std::string("test crs graph"), setup.fem_graph);
    matrix_values_type matrix_values =
      matrix_values_type( "matrix" , setup.fem_graph_length );
    typename matrix_values_type::HostMirror hM =
      Kokkos::create_mirror( matrix_values );

    for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
      for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
        const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

        hM(iEntryFEM) = generate_matrix_coefficient<value_type>(
          setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , block );
      }
    }

    Kokkos::deep_copy( matrix_values , hM );
    matrix[block] = matrix_type("matrix", setup.fem_length, matrix_values,
                                matrix_graph);

    for ( int i = 0 ; i < setup.fem_length ; ++i ) {
      hx(i, block) = generate_vector_coefficient<value_type>(
        setup.fem_length , setup.stoch_length , i , block );
      hy(i, block) = 0.0;
    }

  }

  Kokkos::deep_copy( x , hx );
  Kokkos::deep_copy( y , hy );

  // Original matrix-free multiply algorithm using a block apply
  typedef typename UnitTestSetup<Device>::Cijk_type::k_iterator k_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kj_iterator kj_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kji_iterator kji_iterator;
  k_iterator k_begin = setup.Cijk->k_begin();
  k_iterator k_end = setup.Cijk->k_end();
  for (k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
    int nj = setup.Cijk->num_j(k_it);
    if (nj > 0) {
      int k = index(k_it);
      kj_iterator j_begin = setup.Cijk->j_begin(k_it);
      kj_iterator j_end = setup.Cijk->j_end(k_it);
      unsigned jdx = 0;
      for (kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
        int j = index(j_it);
        vec_type xx = Kokkos::subview( x, Kokkos::ALL(), j );
        vec_type tt = Kokkos::subview( tmp_x, Kokkos::ALL(), jdx++ );
        Kokkos::deep_copy(tt, xx);
      }
      multi_vec_type tmp_x_view =
        Kokkos::subview( tmp_x, Kokkos::ALL(),
                                         std::make_pair(0u,jdx));
      multi_vec_type tmp_y_view =
        Kokkos::subview( tmp_y, Kokkos::ALL(),
                                         std::make_pair(0u,jdx));
      KokkosSparse::spmv(  "N", value_type(1.0), matrix[k] , tmp_x_view , value_type(0.0) , tmp_y_view );
      jdx = 0;
      for (kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
        vec_type tmp_y_view =
          Kokkos::subview( tmp_y, Kokkos::ALL(), jdx++ );
        kji_iterator i_begin = setup.Cijk->i_begin(j_it);
        kji_iterator i_end = setup.Cijk->i_end(j_it);
        for (kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
          int i = index(i_it);
          value_type c = value(i_it);
          vec_type y_view = Kokkos::subview( y, Kokkos::ALL(), i );
          //Stokhos::update( value_type(1.0) , y_view , c , tmp_y_view );
          KokkosBlas::update(c, tmp_y_view, value_type(1.0), y_view, value_type(0.0), y_view);
        }
      }
    }
  }

  Kokkos::deep_copy( hy , y );
  bool success = setup.test_original(hy, out);

  return success;
}

#endif

template< typename ScalarType , class Device >
bool
test_crs_dense_block(const UnitTestSetup<Device>& setup,
                     Teuchos::FancyOStream& out)
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type**, Kokkos::LayoutLeft , Device > block_vector_type ;

  //------------------------------

  typedef Stokhos::BlockCrsMatrix< Stokhos::SymmetricDiagonalSpec< Device > , value_type , Device > matrix_type ;

  typedef typename matrix_type::graph_type  graph_type ;

  // Convert sparse Cijk to dense for faster assembly
  typedef typename UnitTestSetup<Device>::Cijk_type::k_iterator k_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kj_iterator kj_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kji_iterator kji_iterator;
  Stokhos::Dense3Tensor<int,double> dense_cijk(setup.stoch_length);
  for (k_iterator k_it=setup.Cijk->k_begin();
       k_it!=setup.Cijk->k_end(); ++k_it) {
    int k = index(k_it);
    for (kj_iterator j_it = setup.Cijk->j_begin(k_it);
         j_it != setup.Cijk->j_end(k_it); ++j_it) {
      int j = index(j_it);
      for (kji_iterator i_it = setup.Cijk->i_begin(j_it);
           i_it != setup.Cijk->i_end(j_it); ++i_it) {
        int i = index(i_it);
        double c = value(i_it);
        dense_cijk(i,j,k) = c;
      }
    }
  }

  //------------------------------
  // Generate input multivector:

  block_vector_type x =
    block_vector_type( "x" , setup.stoch_length , setup.fem_length );
  block_vector_type y =
    block_vector_type( "y" , setup.stoch_length , setup.fem_length );

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror( x );

  for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
    for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
      hx(iColStoch,iColFEM) =
        generate_vector_coefficient<ScalarType>(
          setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------
  // Generate CRS matrix of blocks with symmetric diagonal storage

  matrix_type matrix ;

  matrix.block =
    Stokhos::SymmetricDiagonalSpec< Device >( setup.stoch_length );
  matrix.graph = Kokkos::create_staticcrsgraph<graph_type>(
    std::string("test crs graph") , setup.fem_graph );
  matrix.values = block_vector_type(
    "matrix" , matrix.block.matrix_size() , setup.fem_graph_length );

  {
    typename block_vector_type::HostMirror hM =
      Kokkos::create_mirror( matrix.values );

    for ( int iRowStoch = 0 ; iRowStoch < setup.stoch_length ; ++iRowStoch ) {
      for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {

        for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
          const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM];

          for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {

            const size_t offset =
              matrix.block.matrix_offset( iRowStoch , iColStoch );

            ScalarType value = 0 ;

            for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
              value += dense_cijk( iRowStoch , iColStoch , k ) *
                generate_matrix_coefficient<ScalarType>(
                  setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , k );
            }

            hM( offset , iEntryFEM ) = value ;
          }
        }

      }
    }

    Kokkos::deep_copy( matrix.values , hM );
  }

  //------------------------------

  Stokhos::multiply( matrix , x , y );

  typename block_vector_type::HostMirror hy = Kokkos::create_mirror( y );
  Kokkos::deep_copy( hy , y );

  bool success = setup.test_commuted(hy, out);

  return success;
}

template< typename ScalarType , class Device >
bool
test_crs_flat_commuted(const UnitTestSetup<Device>& setup,
                       Teuchos::FancyOStream& out)
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type* , Device > vector_type ;

  //------------------------------

  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::graph_type matrix_graph_type;

  // Build stochastic graph
  std::vector< std::vector< int > > stoch_graph( setup.stoch_length );
  Teuchos::RCP<Epetra_CrsGraph> cijk_graph = Stokhos::sparse3Tensor2CrsGraph(
    *setup.basis, *setup.Cijk, *setup.globalComm);
  for ( int i = 0 ; i < setup.stoch_length ; ++i ) {
    int len = cijk_graph->NumGlobalIndices(i);
    stoch_graph[i].resize(len);
    int len2;
    cijk_graph->ExtractGlobalRowCopy(i, len, len2, &stoch_graph[i][0]);
  }

  // Convert sparse Cijk to dense for faster assembly in debug builds
  typedef typename UnitTestSetup<Device>::Cijk_type::k_iterator k_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kj_iterator kj_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kji_iterator kji_iterator;
  Stokhos::Dense3Tensor<int,double> dense_cijk(setup.stoch_length);
  for (k_iterator k_it=setup.Cijk->k_begin();
       k_it!=setup.Cijk->k_end(); ++k_it) {
    int k = index(k_it);
    for (kj_iterator j_it = setup.Cijk->j_begin(k_it);
         j_it != setup.Cijk->j_end(k_it); ++j_it) {
      int j = index(j_it);
      for (kji_iterator i_it = setup.Cijk->i_begin(j_it);
           i_it != setup.Cijk->i_end(j_it); ++i_it) {
        int i = index(i_it);
        double c = value(i_it);
        dense_cijk(i,j,k) = c;
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

  typename vector_type::HostMirror hx = Kokkos::create_mirror( x );

  for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
    for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
      hx(iColStoch + iColFEM*setup.stoch_length) =
        generate_vector_coefficient<ScalarType>(
          setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = Kokkos::create_staticcrsgraph<matrix_graph_type>(
    std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.extent(0);

  matrix.values = matrix_values_type( "matrix" , flat_graph_length );
  {
    typename matrix_values_type::HostMirror hM =
      Kokkos::create_mirror( matrix.values );

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
              setup.fem_length , setup.stoch_length , iRowFEM, iColFEM, k );
          value += dense_cijk(iRowStoch,iColStoch,k) * A_fem_k ;
        }
        hM( iEntry ) = value ;
      }
    }

    Kokkos::deep_copy( matrix.values , hM );
  }

  //------------------------------


  Stokhos::multiply( matrix , x , y );

  typename vector_type::HostMirror hy = Kokkos::create_mirror( y );
  Kokkos::deep_copy( hy , y );

  bool success = setup.test_commuted_flat(hy, out);
  return success;
}

template< typename ScalarType , class Device >
bool
test_crs_flat_original(const UnitTestSetup<Device>& setup,
                       Teuchos::FancyOStream& out)
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type* , Device > vector_type ;

  //------------------------------

  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::graph_type matrix_graph_type;

  // Build stochastic graph
  std::vector< std::vector< int > > stoch_graph( setup.stoch_length );
  Teuchos::RCP<Epetra_CrsGraph> cijk_graph = Stokhos::sparse3Tensor2CrsGraph(
    *setup.basis, *setup.Cijk, *setup.globalComm);
  for ( int i = 0 ; i < setup.stoch_length ; ++i ) {
    int len = cijk_graph->NumGlobalIndices(i);
    stoch_graph[i].resize(len);
    int len2;
    cijk_graph->ExtractGlobalRowCopy(i, len, len2, &stoch_graph[i][0]);
  }

  // Convert sparse Cijk to dense for faster assembly in debug builds
  typedef typename UnitTestSetup<Device>::Cijk_type::k_iterator k_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kj_iterator kj_iterator;
  typedef typename UnitTestSetup<Device>::Cijk_type::kji_iterator kji_iterator;
  Stokhos::Dense3Tensor<int,double> dense_cijk(setup.stoch_length);
  for (k_iterator k_it=setup.Cijk->k_begin();
       k_it!=setup.Cijk->k_end(); ++k_it) {
    int k = index(k_it);
    for (kj_iterator j_it = setup.Cijk->j_begin(k_it);
         j_it != setup.Cijk->j_end(k_it); ++j_it) {
      int j = index(j_it);
      for (kji_iterator i_it = setup.Cijk->i_begin(j_it);
           i_it != setup.Cijk->i_end(j_it); ++i_it) {
        int i = index(i_it);
        double c = value(i_it);
        dense_cijk(i,j,k) = c;
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

  typename vector_type::HostMirror hx = Kokkos::create_mirror( x );

  for ( size_t iCol = 0 ; iCol < flat_length ; ++iCol ) {
    const int iColStoch = iCol / setup.fem_length ;
    const int iColFEM   = iCol % setup.fem_length ;

    hx(iCol) = generate_vector_coefficient<ScalarType>(
      setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = Kokkos::create_staticcrsgraph<matrix_graph_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.extent(0);

  matrix.values = matrix_values_type( "matrix" , flat_graph_length );
  {
    typename matrix_values_type::HostMirror hM =
      Kokkos::create_mirror( matrix.values );

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
              iRowFEM , iColFEM , k );
          value += dense_cijk(iRowStoch,iColStoch,k) * A_fem_k ;
        }
        hM( iEntry ) = value ;

      }

    }

    Kokkos::deep_copy( matrix.values , hM );
  }

  Stokhos::multiply( matrix , x , y );

  typename vector_type::HostMirror hy = Kokkos::create_mirror( y );
  Kokkos::deep_copy( hy , y );

  bool success = setup.test_original_flat(hy, out);
  return success;
}

template< typename ScalarType , typename TensorType, class Device >
bool test_crs_product_tensor(
  const UnitTestSetup<Device>& setup,
  Teuchos::FancyOStream& out,
  const Teuchos::ParameterList& params = Teuchos::ParameterList()) {
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type** , Kokkos::LayoutLeft ,
    Device > block_vector_type ;

  typedef Stokhos::StochasticProductTensor< value_type , TensorType , Device > tensor_type ;

  typedef Stokhos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  //------------------------------
  // Generate input multivector:

  block_vector_type x =
    block_vector_type( "x" , setup.stoch_length_aligned , setup.fem_length );
  block_vector_type y =
    block_vector_type( "y" , setup.stoch_length_aligned , setup.fem_length );

  typename block_vector_type::HostMirror hx =
    Kokkos::create_mirror( x );

  for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
    for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
      hx(setup.perm[iColStoch],iColFEM) =
        generate_vector_coefficient<ScalarType>(
          setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.block =
    Stokhos::create_stochastic_product_tensor< TensorType >( *setup.basis,
                                                             *setup.Cijk,
                                                             params);

  matrix.graph = Kokkos::create_staticcrsgraph<graph_type>(
    std::string("test crs graph") , setup.fem_graph );

  matrix.values = block_vector_type(
    "matrix" , setup.stoch_length_aligned , setup.fem_graph_length );

  typename block_vector_type::HostMirror hM =
    Kokkos::create_mirror( matrix.values );

  for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
    for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
      const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

      for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
        hM(setup.perm[k],iEntryFEM) =
          generate_matrix_coefficient<ScalarType>(
            setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , k );
        //hM(k,iEntryFEM) = 1.0;
      }
    }
  }

  Kokkos::deep_copy( matrix.values , hM );

  //------------------------------

  Stokhos::multiply( matrix , x , y );

  typename block_vector_type::HostMirror hy = Kokkos::create_mirror( y );
  Kokkos::deep_copy( hy , y );

  bool success = setup.test_commuted_perm(hy, out);
  //bool success = true;
  return success;
}

template< typename ScalarType , class Device , int BlockSize >
bool test_linear_tensor(const UnitTestSetup<Device>& setup,
                        Teuchos::FancyOStream& out,
                        const bool symmetric) {
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type** , Kokkos::LayoutLeft ,
    Device > block_vector_type ;
  typedef Stokhos::LinearSparse3Tensor<value_type,Device,BlockSize> TensorType;
  typedef Stokhos::StochasticProductTensor< value_type , TensorType , Device > tensor_type ;

  typedef Stokhos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  // Build tensor
  matrix_type matrix ;

  Teuchos::ParameterList params;
  params.set("Symmetric",symmetric);
  matrix.block =
    Stokhos::create_stochastic_product_tensor< TensorType >( *setup.basis,
                                                             *setup.Cijk,
                                                             params );
  int aligned_stoch_length = matrix.block.tensor().aligned_dimension();

  //------------------------------
  // Generate input multivector:

  block_vector_type x =
    block_vector_type( "x" , aligned_stoch_length , setup.fem_length );
  block_vector_type y =
    block_vector_type( "y" , aligned_stoch_length , setup.fem_length );

  typename block_vector_type::HostMirror hx =
    Kokkos::create_mirror( x );

  for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
    for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
      hx(iColStoch,iColFEM) =
        generate_vector_coefficient<ScalarType>(
          setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
      //hx(iColStoch,iColFEM) = 1.0;
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  matrix.graph = Kokkos::create_staticcrsgraph<graph_type>(
    std::string("test crs graph") , setup.fem_graph );

  matrix.values = block_vector_type(
    "matrix" , aligned_stoch_length , setup.fem_graph_length );

  typename block_vector_type::HostMirror hM =
    Kokkos::create_mirror( matrix.values );

  for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
    for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
      const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

      for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
        hM(k,iEntryFEM) =
          generate_matrix_coefficient<ScalarType>(
            setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , k );
        //hM(k,iEntryFEM) = 1.0;
      }
    }
  }

  Kokkos::deep_copy( matrix.values , hM );

  //------------------------------

  Stokhos::multiply( matrix , x , y );

  typename block_vector_type::HostMirror hy = Kokkos::create_mirror( y );
  Kokkos::deep_copy( hy , y );

  bool success = setup.test_commuted(hy, out);
  //bool success = true;
  return success;
}

template< typename ScalarType , class Device >
bool test_lexo_block_tensor(const UnitTestSetup<Device>& setup,
                            Teuchos::FancyOStream& out) {
  typedef ScalarType value_type ;
  typedef int ordinal_type;
  typedef Kokkos::View< value_type** , Kokkos::LayoutLeft ,
    Device > block_vector_type ;
  typedef Stokhos::LexicographicBlockSparse3Tensor<value_type,Device> TensorType;
  typedef Stokhos::StochasticProductTensor< value_type , TensorType , Device > tensor_type ;

  typedef Stokhos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  //------------------------------
  // Generate input multivector:

  block_vector_type x =
    block_vector_type( "x" , setup.stoch_length , setup.fem_length );
  block_vector_type y =
    block_vector_type( "y" , setup.stoch_length , setup.fem_length );

  typename block_vector_type::HostMirror hx =
    Kokkos::create_mirror( x );

  for ( int iColFEM = 0 ;   iColFEM < setup.fem_length ;   ++iColFEM ) {
    for ( int iColStoch = 0 ; iColStoch < setup.stoch_length ; ++iColStoch ) {
      hx(iColStoch,iColFEM) =
        generate_vector_coefficient<ScalarType>(
          setup.fem_length , setup.stoch_length , iColFEM , iColStoch );
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  /*
    typedef UnitTestSetup<Device>::abstract_basis_type abstract_basis_type;
    typedef UnitTestSetup<Device>::basis_type basis_type;
    typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
    typedef Stokhos::TotalOrderBasis<ordinal_type,value_type,less_type> product_basis_type;
    Teuchos::Array< Teuchos::RCP<const abstract_basis_type> > bases(setup.d);
    for (int i=0; i<setup.d; i++)
      bases[i] = rcp(new basis_type(setup.p,true));
    product_basis_type basis(bases, 1e-12);
  */
  const bool symmetric = true;
  Teuchos::RCP< Stokhos::LTBSparse3Tensor<ordinal_type, value_type> > Cijk =
    Stokhos::computeTripleProductTensorLTBBlockLeaf(*setup.basis, symmetric);

  matrix.block =
    Stokhos::create_stochastic_product_tensor< TensorType >( *setup.basis,
                                                             *Cijk );

  matrix.graph = Kokkos::create_staticcrsgraph<graph_type>(
    std::string("test crs graph") , setup.fem_graph );

  matrix.values = block_vector_type(
    "matrix" , setup.stoch_length , setup.fem_graph_length );

  typename block_vector_type::HostMirror hM =
    Kokkos::create_mirror( matrix.values );

  for ( int iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < setup.fem_length ; ++iRowFEM ) {
    for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < setup.fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
      const int iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM] ;

      for ( int k = 0 ; k < setup.stoch_length ; ++k ) {
        hM(k,iEntryFEM) =
          generate_matrix_coefficient<ScalarType>(
            setup.fem_length , setup.stoch_length , iRowFEM , iColFEM , k );
      }
    }
  }

  Kokkos::deep_copy( matrix.values , hM );

  //------------------------------

  Stokhos::multiply( matrix , x , y );

  typename block_vector_type::HostMirror hy = Kokkos::create_mirror( y );
  Kokkos::deep_copy( hy , y );

  bool success = setup.test_commuted(hy, out);
  return success;
}

}

#endif
