// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include <string>
#include <iostream>
#include <cstdlib>

#include "Teuchos_StandardCatchMacros.hpp"

#include "Stokhos_Epetra.hpp"
#include "Stokhos_Sparse3TensorUtilities.hpp"
#include "EpetraExt_BlockUtility.h"

#include "Kokkos_Timer.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

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

void
run_test(const int p, const int d, const int nGrid, const int nIter,
         const RCP<const Epetra_Comm>& globalComm,
         const RCP<const Epetra_Map>& map,
         const RCP<Epetra_CrsGraph>& graph)
{
  typedef double value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::LegendreBasis<int,value_type> basis_type;
  typedef Stokhos::LexographicLess< Stokhos::MultiIndex<int> > less_type;
  typedef Stokhos::TotalOrderBasis<int,value_type,less_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  // Create Stochastic Galerkin basis and expansion
  Array< RCP<const abstract_basis_type> > bases(d);
  for (int i=0; i<d; i++)
    bases[i] = rcp(new basis_type(p,true));
  RCP< product_basis_type> basis = rcp(new product_basis_type(bases, 1e-12));
  int stoch_length = basis->size();
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

  // Generate Epetra objects
  RCP<const Epetra_Map> sg_map =
    rcp(EpetraExt::BlockUtility::GenerateBlockMap(
          *map, *stoch_row_map, *sg_comm));
  RCP<ParameterList> sg_op_params = rcp(new ParameterList);
  RCP<Stokhos::MatrixFreeOperator> sg_A =
    rcp(new Stokhos::MatrixFreeOperator(sg_comm, basis, epetraCijk,
                                        map, map, sg_map, sg_map,
                                        sg_op_params));
  RCP<Epetra_BlockMap> sg_A_overlap_map =
    rcp(new Epetra_LocalMap(
          stoch_length, 0, *(sg_parallel_data->getStochasticComm())));
  RCP< Stokhos::EpetraOperatorOrthogPoly > A_sg_blocks =
    rcp(new Stokhos::EpetraOperatorOrthogPoly(
          basis, sg_A_overlap_map, map, map, sg_map, sg_comm));
  for (int i=0; i<stoch_length; i++) {
    RCP<Epetra_CrsMatrix> A = rcp(new Epetra_CrsMatrix(Copy, *graph));
    A->PutScalar(1.0);
    A->FillComplete();
    A_sg_blocks->setCoeffPtr(i, A);
  }
  sg_A->setupOperator(A_sg_blocks);

  RCP<Stokhos::EpetraVectorOrthogPoly> sg_x =
    rcp(new Stokhos::EpetraVectorOrthogPoly(
          basis, stoch_row_map, map, sg_map, sg_comm));
  RCP<Stokhos::EpetraVectorOrthogPoly> sg_y =
    rcp(new Stokhos::EpetraVectorOrthogPoly(
          basis, stoch_row_map, map, sg_map, sg_comm));
  sg_x->init(1.0);
  sg_y->init(0.0);

  // Apply operator
  Kokkos::Timer clock;
  for (int iter=0; iter<nIter; ++iter)
    sg_A->Apply( *(sg_x->getBlockVector()), *(sg_y->getBlockVector()) );

  const double t = clock.seconds() / ((double) nIter );
  const double flops = sg_A->countApplyFlops();
  const double gflops = 1.0e-9 * flops / t;

  if (globalComm->MyPID() == 0)
    std::cout << nGrid << " , "
              << d << " , "
              << p << " , "
              << t << " , "
              << gflops << " , "
              << std::endl;
}

int main(int argc, char *argv[])
{
  bool success = true;

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  try {

// Create a communicator for Epetra objects
#ifdef HAVE_MPI
    RCP<const Epetra_Comm> globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    RCP<const Epetra_Comm> globalComm = rcp(new Epetra_SerialComm);
#endif

    // Print header
    if (globalComm->MyPID() == 0)
      std::cout << std::endl
                << "\"#nGrid\" , "
                << "\"#Variable\" , "
                << "\"PolyDegree\" , "
                << "\"MXV Time\" , "
                << "\"MXV GFLOPS\" , "
                << std::endl;

    const int nIter = 1;
    const int nGrid = 32;

    // Generate FEM graph:
    const int fem_length = nGrid * nGrid * nGrid;
    RCP<const Epetra_Map> map = rcp(new Epetra_Map(fem_length, 0, *globalComm));
    std::vector< std::vector<int> > fem_graph;
    generate_fem_graph(nGrid, fem_graph);
    RCP<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy, *map, 27));
    int *my_GIDs = map->MyGlobalElements();
    int num_my_GIDs = map->NumMyElements();
    for (int i=0; i<num_my_GIDs; ++i) {
      int row = my_GIDs[i];
      int num_indices = fem_graph[row].size();
      int *indices = &(fem_graph[row][0]);
      graph->InsertGlobalIndices(row, num_indices, indices);
    }
    graph->FillComplete();

    {
      const int p = 3;
      for (int d=1; d<=12; ++d)
        run_test(p, d, nGrid, nIter, globalComm, map, graph);
    }

    {
      const int p = 5;
      for (int d=1; d<=6; ++d)
        run_test(p, d, nGrid, nIter, globalComm, map, graph);
    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (!success)
    return -1;
  return 0 ;
}
