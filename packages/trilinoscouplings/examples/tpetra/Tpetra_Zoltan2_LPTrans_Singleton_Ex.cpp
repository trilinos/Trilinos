// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//--------------------------------------------------------------------
// This file is a self-contained example of rebalancing the contents of
// a matrix and a vector among different MPI nodes. The matrix is formed,
// on purpose, in an unbalanced way.
//
// Originally, this example used EpetraExt + Isorropia for rebalance. It
// was a modified version of the matrix_1.cpp example in Isorropia.
//
// This example now uses Tpetra + Zoltan2, combining some aspects of the
// original example using Isorropia with aspects from example/graph/graph.cpp
// in Zoltan2.
//--------------------------------------------------------------------

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <TrilinosCouplings_Rebalance_LinearProblem.hpp>
#include <Tpetra_CrsSingletonFilter_LinearProblem.hpp>
#include <Tpetra_SolverMap_LinearProblem.hpp>
#include <Tpetra_Reindex_LinearProblem.hpp>

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <MatrixMarket_Tpetra.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

bool runExample(const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  int numProcs( comm->getSize() );
  int localProc( comm->getRank() );

  using Map_t      = Tpetra::Map<>;
  using LocalId_t  = Map_t::local_ordinal_type;
  using GlobalId_t = Map_t::global_ordinal_type;
  using Node_t     = Map_t::node_type;
  using Scalar_t   = Tpetra::Details::DefaultTypes::scalar_type;
  using MultiV_t   = Tpetra::MultiVector<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Matrix_t   = Tpetra::CrsMatrix<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Problem_t  = Tpetra::LinearProblem<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using SolverMap_t= Tpetra::SolverMap_LinearProblem<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Reindex_t  = Tpetra::Reindex_LinearProblem<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using CrsFilter_t= Tpetra::CrsSingletonFilter_LinearProblem<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Rebalance_t= TrilinosCouplings::Rebalance_LinearProblem<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Reader_t   = Tpetra::MatrixMarket::Reader<Matrix_t>;

  // ****************************************************************
  // Step 1/7: load a linear system to redistribute
  //            over processors.
  // ****************************************************************
  RCP<Matrix_t> A_Original;
  RCP<MultiV_t> LHS_Original, RHS_Original;

  A_Original                      = Reader_t::readSparseFile("SF2_Matrix_Original.mm", comm);
  RCP<const Map_t> A_Original_Map = A_Original->getRangeMap();
  LHS_Original                    = Reader_t::readDenseFile("SF2_LHS_Original.mm", comm, A_Original_Map);
  RHS_Original                    = Reader_t::readDenseFile("SF2_RHS_Original.mm", comm, A_Original_Map);

  bool verbose                             = true;
  bool run_on_host                         = false;
  RCP<MultiV_t> x                     = rcp(new MultiV_t(A_Original_Map, LHS_Original->getNumVectors()));
  RCP<Problem_t> preSingletonProblem  = rcp(new Problem_t(A_Original, x, RHS_Original));

  // ****************************************************************
  // Step 2/7: create and apply singleton filter
  // ****************************************************************
  CrsFilter_t SingletonTransform(run_on_host, verbose);
  RCP<Problem_t> postSingletonProblem = SingletonTransform(preSingletonProblem);
  SingletonTransform.fwd();

  // ****************************************************************
  // Step 3/7: create and apply solver map transform
  // ****************************************************************
  SolverMap_t SolverMapTransform;
  RCP<Problem_t> postSolverMapProblem = SolverMapTransform(postSingletonProblem);
  SolverMapTransform.fwd();

  // ****************************************************************
  // Step 4/7: create and apply reindex transform
  // ****************************************************************
  RCP<Map_t const> newRowMap(Teuchos::null);
  Reindex_t ReindexTransform(newRowMap);
  RCP<Problem_t> postReindexProblem = ReindexTransform(postSolverMapProblem);
  ReindexTransform.fwd();

  // ****************************************************************
  // Step 5/7: create and apply the rebalance transform
  // ****************************************************************
  std::string partitioningMethod("block"); // "scotch" "parmetis"
  RCP<Teuchos::ParameterList> paramListForZoltan2PartitioningProblem = rcp<Teuchos::ParameterList>( new Teuchos::ParameterList() );
  paramListForZoltan2PartitioningProblem->set("partitioning_approach", "partition");
  paramListForZoltan2PartitioningProblem->set("algorithm", partitioningMethod);

  Rebalance_t RebalanceTransform( paramListForZoltan2PartitioningProblem );
  RCP<Problem_t> transformedLP = RebalanceTransform(postReindexProblem);
  RebalanceTransform.fwd();

  RCP<Matrix_t> rebalancedMatrix = rcp<Matrix_t>( dynamic_cast<Matrix_t *>(transformedLP->getMatrix().get()), false );
  RCP<MultiV_t> rebalancedLhs    = transformedLP->getLHS();
  RCP<MultiV_t> rebalancedRhs    = transformedLP->getRHS();

  // ****************************************************************
  // Step 6/7: Copy solution to reduced system and reverse transform
  // ****************************************************************
  auto reducedRowMap                 = SingletonTransform.ReducedMatrix()->getRowMap();
  RCP<MultiV_t> reducedSolution = Reader_t::readDenseFile("SF2_LHS_Reduced.mm", comm, reducedRowMap);
  
  rebalancedLhs->update(1.0, *reducedSolution, 0.0);

  // Rebalance Transform Reverse
  RebalanceTransform.rvs();
  
  // Reindex Transform Reverse
  ReindexTransform.rvs();
  
  // Solver Map Transform Reverse
  SolverMapTransform.rvs();
  
  // Singleton Transform Reverse
  SingletonTransform.rvs();

  // ****************************************************************
  // Step 7/7: check that originalRhs == originalMat * originalLhs
  // ****************************************************************
  int numVectors = RHS_Original->getNumVectors(); 
  std::vector<Scalar_t> originalRhsNorms2_vec( numVectors );
  Teuchos::ArrayView<Scalar_t> originalRhsNorms2_array = Teuchos::arrayViewFromVector( originalRhsNorms2_vec );
  RHS_Original->norm2(originalRhsNorms2_array);

  Teuchos::RCP<MultiV_t> originalMatLhs = Tpetra::createMultiVector<Scalar_t,LocalId_t,GlobalId_t>( A_Original_Map, numVectors );
  A_Original->apply(*x, *originalMatLhs);
  MultiV_t diff(*RHS_Original);
  diff.update(-1., *originalMatLhs, 1.); 

  std::vector<Scalar_t> diffNorms2_vec( numVectors );
  Teuchos::ArrayView<Scalar_t> diffNorms2_array = Teuchos::arrayViewFromVector( diffNorms2_vec );
  diff.norm2(diffNorms2_array);

  bool isOkay = true;
  for (size_t v(0); v < numVectors; ++v) {
    Scalar_t diffNorm2       = diffNorms2_array[v];
    Scalar_t originalNorm2   = originalRhsNorms2_array[v];

    Scalar_t ratio = diffNorm2 / originalNorm2;

    if (localProc == 0) {
      std::cout << "||originalRhs[" << v << "] - originalMat * originalLhs[" << v << "]||_2 = " << diffNorms2_array[v]
                << ", ||originalRhs[" << v << "]||_2 = " << originalNorm2
                << ", ratio = " << ratio
                << ", Scalar_t = " << typeid(Scalar_t).name()
                << ", eps = " << Teuchos::ScalarTraits<Scalar_t>::eps()
                << std::endl;
      std::cout.flush();
    }
    isOkay &= ( std::fabs(ratio) < 10. * Teuchos::ScalarTraits<Scalar_t>::eps() );
  }

  return isOkay;
}

int main(int argc, char** argv) {

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  { // Avoiding objects to persist after MPI_Finalize or Kokkos::finalize
    Teuchos::RCP< const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    if ( runExample(comm) ) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
    else {
      std::cout << "End Result: TEST FAILED!" << std::endl;
    }

  }

  return 0;
}

