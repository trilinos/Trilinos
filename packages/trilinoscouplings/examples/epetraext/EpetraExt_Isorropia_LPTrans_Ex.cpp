// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//--------------------------------------------------------------------
//This file is a self-contained example of creating an Epetra_CrsGraph
//and Epetra_CrsMatrix object, and using Isorropia through the EpetraExt
//transform interface to rebalance them.
// NOTE:  This is a modified version of the matrix_1.cpp example in Isorropia.
//--------------------------------------------------------------------

#include <Isorropia_ConfigDefs.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_LinearProblem.h>
#include <EpetraExt_Isorropia_CrsGraph.h>
#include <EpetraExt_LPTrans_From_GraphTrans.h>
#include <Isorropia_Exception.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

//Declaration for helper-function that creates epetra objects. These
//functions are implemented at the bottom of this file.
Teuchos::RCP<Epetra_CrsMatrix>
  create_epetra_matrix(int numProcs, int localProc);


int main(int argc, char** argv) {
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

  int p, numProcs = 1;
  int localProc = 0;

  //first, set up our MPI environment...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  //Now create a balanced copy of the matrix using the EpetraExt
  //transform interface to Isorropia.
  //NOTE: By default, Isorropia will use Zoltan for the
  //repartitioning, if Isorropia was configured with Zoltan support.
  //(i.e., --enable-isorropia-zoltan flag to configure, plus Zoltan include
  //paths and library directives)
  //If Isorropia was not configured with Zoltan support, then a simple
  //built-in linear partitioner will be used to make sure the number
  //of nonzeros on each processor is equal or close to equal.

  Teuchos::ParameterList paramlist;

#ifdef HAVE_ISORROPIA_ZOLTAN

  // If Zoltan is available, we'll specify that the Zoltan package use
  // graph-partitioning for the partitioning operation and specifically
  // 1D hypergraph partitioning, by creating a parameter sublist named
  // "Zoltan" and setting the appropriate values.
  // (See Zoltan documentation for other valid parameters...)

  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_METHOD", "HYPERGRAPH");

#else
  //If Zoltan is not available, we don't need to set any parameters.
#endif

  // Create the object to generate the balanced graph.
  Teuchos::RCP<EpetraExt::Isorropia_CrsGraph> Trans = 
    Teuchos::rcp( new EpetraExt::Isorropia_CrsGraph( paramlist ) );

  // Create a linear problem transform interface.
  Teuchos::RCP<EpetraExt::LinearProblem_GraphTrans> LPTrans =
    Teuchos::rcp( new EpetraExt::LinearProblem_GraphTrans(
    *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(Trans.get())) ) );

  // Create an Epetra_CrsMatrix object.
  Teuchos::RCP<Epetra_CrsMatrix> crsmatrix;
  try {
    crsmatrix = create_epetra_matrix(numProcs, localProc);
  }
  catch(std::exception& e) {
    std::cout << "ERROR: create_epetra_matrix threw exception '"
          << e.what() << "' on proc " << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  // Create Epetra_MultiVectors for the lhs and rhs of the linear problem.
  Epetra_MultiVector lhs(crsmatrix->Map(), 1);
  Epetra_MultiVector rhs(crsmatrix->Map(), 1);
  rhs.PutScalar( 1.0 );

  // Create the linear problem with the original partitioning.
  Epetra_LinearProblem problem( &*crsmatrix, &lhs, &rhs );

  // Create the new linear problem and perform the balanced partitioning.
  // NOTE:  The balanced linear system will be in tProblem after fwd() is called.
  //        It is not necessary for the RCP to manage the transformed problem.
  Teuchos::RCP<Epetra_LinearProblem> tProblem = Teuchos::rcp( &((*LPTrans)( problem )), false );
  LPTrans->fwd();

  int graphrows1 = crsmatrix->NumMyRows();
  int bal_graph_rows = tProblem->GetMatrix()->NumMyRows();
  int graphnnz1 = crsmatrix->NumMyNonzeros();
  int bal_graph_nnz = tProblem->GetMatrix()->NumMyNonzeros();

  for(p=0; p<numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p << ": input matrix local rows: " << graphrows1
       << ", local NNZ: " << graphnnz1 << std::endl;
    std::cout << "proc " << p << ": balanced matrix local rows: "
       << bal_graph_rows << ", local NNZ: " << bal_graph_nnz << std::endl;
  }

  MPI_Finalize();

#else
  std::cout << "ERROR:  This MUST be and MPI build with Epetra enabled!" << std::endl;
#endif

  return(0);
}

//Below are implementations of the helper-functions that create the
//poorly-balanced epetra objects for use in the above example program.

#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

Teuchos::RCP<Epetra_CrsMatrix>
  create_epetra_matrix(int numProcs, int localProc)
{
  if (localProc == 0) {
    std::cout << " creating Epetra_CrsMatrix with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsMatrix with rows spread un-evenly over
  //processors.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int local_num_rows = 200;
  int nnz_per_row = local_num_rows/4+1;
  int global_num_rows = numProcs*local_num_rows;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;

  int adjustment = local_num_rows/2;

  //adjust local_num_rows so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_num_rows -= adjustment;
  }
  else {
    local_num_rows += adjustment;
  }

  //if numProcs is not an even number, undo the local_num_rows adjustment
  //on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_num_rows -= adjustment;
    }
  }

  //now we're ready to create a row-map.
  Epetra_Map rowmap(global_num_rows, local_num_rows, 0, comm);

  //create a matrix
  Teuchos::RCP<Epetra_CrsMatrix> matrix =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row));

  std::vector<int> indices(nnz_per_row);
  std::vector<double> coefs(nnz_per_row);

  int err = 0;

  for(int i=0; i<local_num_rows; ++i) {
    int global_row = rowmap.GID(i);
    int first_col = global_row - nnz_per_row/2;

    if (first_col < 0) {
      first_col = 0;
    }
    else if (first_col > (global_num_rows - nnz_per_row)) {
      first_col = global_num_rows - nnz_per_row;
    }

    for(int j=0; j<nnz_per_row; ++j) {
      indices[j] = first_col + j;
      coefs[j] = 1.0;
    }

    int err = matrix->InsertGlobalValues(global_row, nnz_per_row,
                                         &coefs[0], &indices[0]);
    if (err < 0) {
      err = matrix->ReplaceGlobalValues(global_row, nnz_per_row,
                                        &coefs[0], &indices[0]);
      if (err < 0) {
        throw Isorropia::Exception("create_epetra_matrix: error inserting matrix values.");
      }
    }
  }

  err = matrix->FillComplete();
  if (err != 0) {
    throw Isorropia::Exception("create_epetra_matrix: error in matrix.FillComplete()");
  }

  return(matrix);
}

#endif //HAVE_MPI && HAVE_EPETRA

