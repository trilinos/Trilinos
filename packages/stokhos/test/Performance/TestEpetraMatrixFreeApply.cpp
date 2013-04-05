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

#include <string>
#include <iostream>
#include <cstdlib>

#ifdef EPETRA_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"
#include "EpetraExt_BlockUtility.h"

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

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
			   Teuchos::Array< Teuchos::Array<ordinal> > & graph )
{
  graph.resize( N * N * N , Teuchos::Array<ordinal>() );

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

Teuchos::Array<double>
test_original_matrix_free_epetra(
  const Teuchos::Array<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag ,
  const bool test_block ,
  const bool check )
{
  typedef double value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::LegendreBasis<int,value_type> basis_type;
  typedef Stokhos::CompletePolynomialBasis<int,value_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;
  using Teuchos::ParameterList;

   // Create a communicator for Epetra objects
  RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
  RCP<const Epetra_MpiComm> globalMpiComm = 
    rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  //globalMpiComm->Print(std::cout);
  globalComm = globalMpiComm;
#else
   RCP<const Epetra_SerialComm> globalSerialComm =
     rcp(new Epetra_SerialComm);
   //globalSerialComm->Print(std::cout);
   globalComm = globalSerialComm;
#endif

  //int MyPID = globalComm->MyPID();

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL); 
  for (size_t i=0; i<num_KL; i++)
    bases[i] = rcp(new basis_type(var_degree[i],true));
  RCP<const product_basis_type> basis = 
    rcp(new product_basis_type(bases, 1e-12));
  const size_t stoch_length = basis->size();
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor();

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
  // Generate FEM graph:

  Teuchos::Array< Teuchos::Array<int> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;
  generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  // Generate Epetra objects
  RCP<const Epetra_Map> x_map = rcp(new Epetra_Map(static_cast<int>(fem_length), 0, *app_comm));
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
  int nnz = graph->NumGlobalNonzeros();
  
  RCP<ParameterList> sg_op_params = rcp(new ParameterList);
  RCP<Stokhos::MatrixFreeOperator> sg_A = 
    rcp(new Stokhos::MatrixFreeOperator(sg_comm, basis, epetraCijk, 
					x_map, x_map, sg_x_map, sg_x_map, 
					sg_op_params));
  RCP<Epetra_BlockMap> sg_A_overlap_map =
    rcp(new Epetra_LocalMap(
	  static_cast<int>(stoch_length), 0, *(sg_parallel_data->getStochasticComm())));
  RCP< Stokhos::EpetraOperatorOrthogPoly > A_sg_blocks = 
    rcp(new Stokhos::EpetraOperatorOrthogPoly(
	  basis, sg_A_overlap_map, x_map, x_map, sg_x_map, sg_comm));
  for (unsigned int i=0; i<stoch_length; i++) {
    RCP<Epetra_CrsMatrix> A = rcp(new Epetra_CrsMatrix(Copy, *graph));
    A->FillComplete();
    A->PutScalar(1.0);
    A_sg_blocks->setCoeffPtr(i, A);
  }
  sg_A->setupOperator(A_sg_blocks);

  RCP<Stokhos::EpetraVectorOrthogPoly> sg_x =
    rcp(new Stokhos::EpetraVectorOrthogPoly(
	  basis, stoch_row_map, x_map, sg_x_map, sg_comm));
  RCP<Stokhos::EpetraVectorOrthogPoly> sg_y =
    rcp(new Stokhos::EpetraVectorOrthogPoly(
	  basis, stoch_row_map, x_map, sg_x_map, sg_comm));
  sg_x->init(1.0);
  sg_y->init(0.0);


  // Time apply
  Teuchos::Time clock("apply") ;
  clock.start();
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    sg_A->Apply( *(sg_x->getBlockVector()), *(sg_y->getBlockVector()) );
  }
  clock.stop();
  double seconds_per_iter = 
    clock.totalElapsedTime() / ((double) iterCount );

  // Averge time across proc's
  double average_seconds_per_iter;
  globalComm->SumAll(&seconds_per_iter, &average_seconds_per_iter, 1);
  average_seconds_per_iter /= globalComm->NumProc();

  // Compute number of fem mat-vec's
  int n_apply = 0;
  int n_add = 0;
  for (Cijk_type::k_iterator k_it=Cijk->k_begin(); 
       k_it!=Cijk->k_end(); ++k_it) {
    for (Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	 j_it != Cijk->j_end(k_it); ++j_it) {
      ++n_apply;
      for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it); 
	   i_it != Cijk->i_end(j_it); ++i_it) {
	++n_add;
      }
    }
  }

  const double flops = 1.0e-9*(2.0*static_cast<double>(n_apply)*nnz + 
			       static_cast<double>(n_add)*fem_length);

  //------------------------------

  Teuchos::Array<double> perf(8);
  perf[0] = stoch_length;
  perf[1] = fem_length;
  perf[2] = stoch_length * fem_length;
  perf[3] = Cijk->num_entries();
  perf[4] = nnz;
  perf[5] = average_seconds_per_iter ;
  perf[6] = flops/average_seconds_per_iter;
  perf[7] = flops;

  return perf;
}

void performance_test_driver_epetra( const int pdeg ,
				     const int minvar ,
				     const int maxvar ,
				     const int nGrid ,
				     const int nIter ,
				     const bool print ,
				     const bool test_block ,
				     const bool check ,
				     Teuchos::FancyOStream& out)
{
  out.precision(8);

  //------------------------------

  out << std::endl
      << "\"#nGrid\" , \"NNZ\" "
      << "\"#Variable\" , \"PolyDegree\" , \"#Bases\" , "
      << "\"#TensorEntry\" , "
      << "\"VectorSize\" , "
      << "\"Original-Matrix-Free-Block-MXV-Time\" , "
      << "\"Original-Matrix-Free-Block-MXV-GFLOPS\" , "
      << std::endl ;
  
  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {
    Teuchos::Array<int> var_degree( nvar , pdeg );

    const Teuchos::Array<double> perf_original_mat_free_epetra =
      test_original_matrix_free_epetra( var_degree , nGrid , nIter , print , test_block , check );

    out << nGrid << " , "
	<< perf_original_mat_free_epetra[4] << " , "
	<< nvar << " , " << pdeg << " , "
	<< perf_original_mat_free_epetra[0] << " , "
	<< perf_original_mat_free_epetra[3] << " , "
	<< perf_original_mat_free_epetra[2] << " , "
	<< perf_original_mat_free_epetra[5] << " , "
	<< perf_original_mat_free_epetra[6] << " , "
	<< std::endl ;
    
  }

    
}

int main(int argc, char *argv[])
{
  bool success = true;

  try {
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
    Teuchos::RCP< Teuchos::FancyOStream > out = 
      Teuchos::VerboseObjectBase::getDefaultOStream();
    
    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    int p = 3;
    CLP.setOption("p", &p, "Polynomial order");
    int d_min = 1;
    CLP.setOption("dmin", &d_min, "Starting stochastic dimension");
    int d_max = 12;
    CLP.setOption("dmax", &d_max, "Ending stochastic dimension");
    int nGrid = 64;
    CLP.setOption("n", &nGrid, "Number of spatial grid points in each dimension");
    int nIter = 1;
    CLP.setOption("niter", &nIter, "Number of iterations");
    bool test_block = true;
    CLP.setOption("block", "no-block", &test_block, "Use block algorithm");
    CLP.parse( argc, argv );
    
    bool print = false ;
    bool check = false;
    performance_test_driver_epetra( 
      p , d_min , d_max , nGrid , nIter , print , test_block , check , *out );
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    
  if (!success)
    return -1;
  return 0;
}

