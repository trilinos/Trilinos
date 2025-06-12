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
// This example now uses Tpetra + Zoltan2. It is a modified version of
// the zoltan2/example/graph/graph.cpp example in Zoltan2.
//--------------------------------------------------------------------

#ifdef HAVE_MPI

#include <mpi.h>
#include <Teuchos_DefaultMpiComm.hpp>

#else

#include <Teuchos_DefaultSerialComm.hpp>

#endif // ifdef HAVE_MPI

//#include <Tpetra_Map_decl.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
//#include <Tpetra_MultiVector_decl.hpp>

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#ifdef HAVE_MPI
void runExample() {
  int numProcs(1);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int localProc(0);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);

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
  Teuchos::RCP< Teuchos::MpiComm<int> > comm = Teuchos::rcp( new Teuchos::MpiComm<int>(MPI_COMM_WORLD) );
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

  std::string partitioningMethod("scotch");
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
  //           and after solving the partitionaing problem.
  // ****************************************************************
  localId_t originalLocalNumRows = originalMatrix->getLocalNumRows();
  localId_t balancedLocalNumRows = balancedMatrix->getLocalNumRows();
  localId_t originalLocalNNZ     = originalMatrix->getLocalNumEntries();
  localId_t balancedLocalNNZ     = balancedMatrix->getLocalNumEntries();

  MPI_Barrier(MPI_COMM_WORLD);
  for (int p(0); p < numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p
              << ": original matrix local num rows = " << originalLocalNumRows
              << ", local NNZ = " << originalLocalNNZ
              << std::endl;
    std::cout << "proc " << p
              << ": balanced matrix local num rows = " << balancedLocalNumRows
              << ", local NNZ = " << balancedLocalNNZ
              << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // ****************************************************************
  // Step 6/6: check sparse matvec
  // ****************************************************************
  Teuchos::RCP<Vector_t> originalMatVec = Tpetra::createVector<scalar_t,localId_t,globalId_t>( originalMatrix->getRangeMap() );
  originalMatrix->apply(*originalVector, *originalMatVec);
  scalar_t originalNorm = originalMatVec->norm2();

  Teuchos::RCP<Vector_t> balancedMatVec = Tpetra::createVector<scalar_t,localId_t,globalId_t>( balancedMatrix->getRangeMap() );
  balancedMatrix->apply(*balancedVector, *balancedMatVec);
  scalar_t balancedNorm = balancedMatVec->norm2();

  if (localProc == 0) {
    std::cout << "||originalMatVec||_2 = " << originalNorm
              << ", ||balancedMatVec||_2 = " << balancedNorm
              << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
}
#endif // ifdef HAVE_MPI

int main(int argc, char** argv) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  runExample();
  MPI_Finalize();
#else
  std::cout << "ERROR:  This MUST be a MPI build with Tpetra enabled!" << std::endl;
#endif

  return(0);
}

