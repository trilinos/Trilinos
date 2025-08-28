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

#include <Tpetra_CrsMatrix.hpp>

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

bool runExample(const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  int numProcs( comm->getSize() );
  int localProc( comm->getRank() );

  using Map_t      = Tpetra::Map<>;
  using localId_t  = Map_t::local_ordinal_type;
  using globalId_t = Map_t::global_ordinal_type;
  using scalar_t   = Tpetra::Details::DefaultTypes::scalar_type;
  using Vector_t   = Tpetra::Vector<scalar_t, localId_t, globalId_t>;
  using Matrix_t   = Tpetra::CrsMatrix<scalar_t, localId_t, globalId_t>;

  // ****************************************************************
  // Step 1/6: create a Tpetra::CrsMatrix with rows spread un-evenly
  //           over processors.
  // ****************************************************************
  localId_t  local_num_rows  = 200;
  localId_t  nnz_per_row     = local_num_rows/4+1;
  globalId_t global_num_rows = numProcs*local_num_rows;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;
  localId_t adjustment = local_num_rows/2;

  // Adjust local_num_rows so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_num_rows -= adjustment;
  }
  else {
    local_num_rows += adjustment;
  }

  // If numProcs is not an even number, undo the local_num_rows adjustment
  // on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_num_rows -= adjustment;
    }
  }

  // Now we're ready to create a row-map.
  Teuchos::RCP< Map_t > rowmap = Teuchos::rcp( new Map_t(global_num_rows, local_num_rows, 0, comm) );
  Teuchos::RCP< Matrix_t > originalMatrix = Teuchos::rcp( new Matrix_t(rowmap, nnz_per_row) );

  std::vector<globalId_t> indices(nnz_per_row);
  std::vector<double> coefs(nnz_per_row);

  for (localId_t i(0); i < local_num_rows; ++i) {
    globalId_t global_row = rowmap->getGlobalElement(i);
    globalId_t first_col = global_row - nnz_per_row/2;

    if (first_col < 0) {
      first_col = 0;
    }
    else if (first_col > (global_num_rows - nnz_per_row)) {
      first_col = global_num_rows - nnz_per_row;
    }

    for (localId_t j(0); j < nnz_per_row; ++j) {
      indices[j] = first_col + j;
      coefs[j] = 1.0;
    }

    originalMatrix->insertGlobalValues(global_row,
                                       nnz_per_row,
                                       &coefs[0],
                                       &indices[0]);
  }

  originalMatrix->fillComplete();

  // ****************************************************************
  // Step 2/6: instantiate a partitioning problem and solve it.
  // ****************************************************************
  using MatrixAdapter_t = Zoltan2::XpetraCrsMatrixAdapter<Matrix_t>;
  MatrixAdapter_t matrixAdapter(originalMatrix);

  std::string partitioningMethod("block"); // "scotch" "parmetis"
  Teuchos::ParameterList param;
  param.set("partitioning_approach", "partition");
  param.set("algorithm", partitioningMethod);

  Zoltan2::PartitioningProblem<MatrixAdapter_t> partitioningProblem(&matrixAdapter, &param);
  partitioningProblem.solve();

  // ****************************************************************
  // Step 3/6: create a random vector compatible (distribution wise)
  //           with the original matrix. The original matrix will
  //           multiply this vector in the step 6/6 below.
  // ****************************************************************
  Teuchos::RCP<Vector_t> originalVector;
  originalVector  = Tpetra::createVector<scalar_t,localId_t,globalId_t>( originalMatrix->getDomainMap() );
  originalVector->randomize();

  // ****************************************************************
  // Step 4/6: balance the matrix and the vector
  // ****************************************************************
  Teuchos::RCP<Matrix_t> balancedMatrix;
  matrixAdapter.applyPartitioningSolution(*originalMatrix,
                                          balancedMatrix,
                                          partitioningProblem.getSolution());

  using MultiVectorAdapter_t = Zoltan2::XpetraMultiVectorAdapter<Vector_t>;
  MultiVectorAdapter_t vectorAdapter(originalVector);

  Teuchos::RCP<Vector_t> balancedVector;
  vectorAdapter.applyPartitioningSolution(*originalVector,
                                          balancedVector,
                                          partitioningProblem.getSolution());

  // ****************************************************************
  // Step 5/6: compare the matrix balancing among MPI nodes, before
  //           and after the redistribution.
  // ****************************************************************
  localId_t originalLocalNumRows = originalMatrix->getLocalNumRows();
  localId_t balancedLocalNumRows = balancedMatrix->getLocalNumRows();
  localId_t overallMin_originalLocalNumRows(0);
  localId_t overallMax_originalLocalNumRows(0);
  localId_t overallMin_balancedLocalNumRows(0);
  localId_t overallMax_balancedLocalNumRows(0);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &originalLocalNumRows, &overallMin_originalLocalNumRows);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &originalLocalNumRows, &overallMax_originalLocalNumRows);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &balancedLocalNumRows, &overallMin_balancedLocalNumRows);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &balancedLocalNumRows, &overallMax_balancedLocalNumRows);

  localId_t originalLocalNNZ     = originalMatrix->getLocalNumEntries();
  localId_t balancedLocalNNZ     = balancedMatrix->getLocalNumEntries();
  localId_t overallMin_originalLocalNNZ(0);
  localId_t overallMax_originalLocalNNZ(0);
  localId_t overallMin_balancedLocalNNZ(0);
  localId_t overallMax_balancedLocalNNZ(0);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &originalLocalNNZ, &overallMin_originalLocalNNZ);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &originalLocalNNZ, &overallMax_originalLocalNNZ);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &balancedLocalNNZ, &overallMin_balancedLocalNNZ);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &balancedLocalNNZ, &overallMax_balancedLocalNNZ);

  comm->barrier();
  for (int p(0); p < numProcs; ++p) {
    comm->barrier();

    if (p != localProc) continue;

    std::cout << "proc " << p
              << ": original matrix local num rows = " << originalLocalNumRows
              << "( min = " << overallMin_originalLocalNumRows
              << ", max = " << overallMax_originalLocalNumRows
              << "); original local NNZ = " << originalLocalNNZ
              << "( min = " << overallMin_originalLocalNNZ
              << ", max = " << overallMax_originalLocalNNZ
              << ")"
              << std::endl;
    std::cout << "proc " << p
              << ": balanced matrix local num rows = " << balancedLocalNumRows
              << "( min = " << overallMin_balancedLocalNumRows
              << ", max = " << overallMax_balancedLocalNumRows
              << "); balanced local NNZ = " << balancedLocalNNZ
              << "( min = " << overallMin_balancedLocalNNZ
              << ", max = " << overallMax_balancedLocalNNZ
              << ")"
              << std::endl;
  }

  bool allOk1(true);
  if ((overallMin_originalLocalNumRows <= overallMin_balancedLocalNumRows) &&
      (overallMax_balancedLocalNumRows <= overallMax_originalLocalNumRows) &&
      (overallMin_originalLocalNNZ     <= overallMin_balancedLocalNNZ    ) &&
      (overallMax_balancedLocalNNZ     <= overallMax_originalLocalNNZ    )) {
    // Ok
  }
  else {
    allOk1 = false;
  }

  if (localProc == 0) {
    std::cout << "allOk1 = " << allOk1 << std::endl;
  }
  
  comm->barrier();
  
  // ****************************************************************
  // Step 6/6: check sparse matvec
  // ****************************************************************
  Teuchos::RCP<Vector_t> originalMatVec = Tpetra::createVector<scalar_t,localId_t,globalId_t>( originalMatrix->getRangeMap() );
  originalMatrix->apply(*originalVector, *originalMatVec);
  scalar_t originalNorm = originalMatVec->norm2();

  Teuchos::RCP<Vector_t> balancedMatVec = Tpetra::createVector<scalar_t,localId_t,globalId_t>( balancedMatrix->getRangeMap() );
  balancedMatrix->apply(*balancedVector, *balancedMatVec);
  scalar_t balancedNorm = balancedMatVec->norm2();

  scalar_t relativeDiff = (balancedNorm - originalNorm) / originalNorm;
  
  if (localProc == 0) {
    std::cout << "||originalMatVec||_2 = " << originalNorm
              << ", ||balancedMatVec||_2 = " << balancedNorm
              << ", relativeDiff = " << relativeDiff
              << ", scalar_t = " << typeid(scalar_t).name()
              << ", Teuchos::ScalarTraits<scalar_t>::eps() = " << Teuchos::ScalarTraits<scalar_t>::eps()
              << std::endl;
  }

  bool allOk2( std::fabs(relativeDiff) < 10. * Teuchos::ScalarTraits<scalar_t>::eps() );

  if (localProc == 0) {
    std::cout << "allOk2 = " << allOk2 << std::endl;
  }

  comm->barrier();

  return (allOk1 && allOk2);
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

