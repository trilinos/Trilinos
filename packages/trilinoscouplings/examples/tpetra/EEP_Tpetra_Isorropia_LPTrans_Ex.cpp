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

//#include <Isorropia_ConfigDefs.hpp>

#ifdef HAVE_MPI

#include <mpi.h>
#include <Teuchos_DefaultMpiComm.hpp>

#else

#include <Teuchos_DefaultSerialComm.hpp>

#endif

#include <Tpetra_Map_decl.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_LinearProblem_decl.hpp>
#include <EEP_Tpetra_Isorropia_CrsGraph.hpp> // In 'packages/trilinoscouplings/src/tpetra/'
#include <EEP_Tpetra_LPTrans_From_GraphTrans.hpp> // In 'packages/tpetra/core/src/'

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#ifdef HAVE_MPI
void runExample() {
  std::cout << "Entering runExample()" << std::endl;

  int numProcs = 1;
  int localProc = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  //Now create a balanced copy of the matrix using the Tpetra
  //transform interface to Isorropia.
  //NOTE: By default, Isorropia will use Zoltan for the
  //repartitioning, if Isorropia was configured with Zoltan support.
  //(i.e., --enable-isorropia-zoltan flag to configure, plus Zoltan include
  //paths and library directives)
  //If Isorropia was not configured with Zoltan support, then a simple
  //built-in linear partitioner will be used to make sure the number
  //of nonzeros on each processor is equal or close to equal.

  Teuchos::ParameterList paramlist;

  std::cout << "In main(), pos 000" << std::endl;

#if 1 // def HAVE_ISORROPIA_ZOLTAN // AquiToDo

  // If Zoltan is available, we'll specify that the Zoltan package use
  // graph-partitioning for the partitioning operation and specifically
  // 1D hypergraph partitioning, by creating a parameter sublist named
  // "Zoltan" and setting the appropriate values.
  // (See Zoltan documentation for other valid parameters...)

  //std::string partitioning_method_str("PARTITIONING METHOD"); // Aqui
  //std::string partitioning_method =
  //paramlist_.get(partitioning_method_str, "UNSPECIFIED");

  std::cout << "In main(), pos 001" << std::endl;
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_METHOD", "HYPERGRAPH"); // AquiAqui

#else
  //If Zoltan is not available, we don't need to set any parameters.
#endif

  std::cout << "In main(), pos 002" << std::endl;

  // Create the object to generate the balanced graph.
  Teuchos::RCP< Tpetra::Isorropia_CrsGraph< Tpetra::CrsMatrix<double>::local_ordinal_type
                                          , Tpetra::CrsMatrix<double>::global_ordinal_type
                                          , Tpetra::CrsMatrix<double>::node_type
                                          > > Trans = Teuchos::rcp( new Tpetra::Isorropia_CrsGraph< Tpetra::CrsMatrix<double>::local_ordinal_type
                                                                                                  , Tpetra::CrsMatrix<double>::global_ordinal_type
                                                                                                  , Tpetra::CrsMatrix<double>::node_type
                                                                                                  >( paramlist ) );

  std::cout << "In main(), pos 003" << std::endl;
  
  // Create a linear problem transform interface.
  Teuchos::RCP< Tpetra::LinearProblem_GraphTrans< double
                                                , Tpetra::CrsMatrix<double>::local_ordinal_type
                                                , Tpetra::CrsMatrix<double>::global_ordinal_type
                                                , Tpetra::CrsMatrix<double>::node_type
                                                > > LPTrans;

  std::cout << "In main(), pos 004" << std::endl;

#if 0 // EEP
  LPTrans = Teuchos::rcp( new Tpetra::LinearProblem_GraphTrans( *(dynamic_cast<Tpetra::StructuralSameTypeTransform<Tpetra_CrsGraph>*>(Trans.get())) ) );
#endif // EEP

  // Create a Tpetra::CrsMatrix with rows spread un-evenly over processors.
  if (localProc == 0) {
    std::cout << "Creating Tpetra::CrsMatrix with un-even distribution..."
              << std::endl;
  }

  Teuchos::RCP< Teuchos::MpiComm<int> > comm = Teuchos::rcp( new Teuchos::MpiComm<int>(MPI_COMM_WORLD) );
  Tpetra::CrsMatrix<double>::local_ordinal_type local_num_rows = 200;
  Tpetra::CrsMatrix<double>::local_ordinal_type nnz_per_row = local_num_rows/4+1;
  Tpetra::CrsMatrix<double>::global_ordinal_type global_num_rows = numProcs*local_num_rows;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;

  Tpetra::CrsMatrix<double>::local_ordinal_type adjustment = local_num_rows/2;

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
  Teuchos::RCP< Tpetra::Map<> > rowmap = Teuchos::rcp( new Tpetra::Map<>(global_num_rows, local_num_rows, 0, comm) );
  Teuchos::RCP< Tpetra::CrsMatrix<double> > crsmatrix = Teuchos::rcp( new Tpetra::CrsMatrix<double>(rowmap, nnz_per_row) );

  std::cout << "crsmatrix = " << crsmatrix << std::endl;

  std::vector<Tpetra::CrsMatrix<double>::global_ordinal_type> indices(nnz_per_row);
  std::vector<double> coefs(nnz_per_row);

  for (Tpetra::CrsMatrix<double>::local_ordinal_type i(0); i < local_num_rows; ++i) {
    Tpetra::CrsMatrix<double>::global_ordinal_type global_row = rowmap->getGlobalElement(i);
    Tpetra::CrsMatrix<double>::global_ordinal_type first_col = global_row - nnz_per_row/2;

    if (first_col < 0) {
      first_col = 0;
    }
    else if (first_col > (global_num_rows - nnz_per_row)) {
      first_col = global_num_rows - nnz_per_row;
    }

    for (Tpetra::CrsMatrix<double>::local_ordinal_type j(0); j < nnz_per_row; ++j) {
      indices[j] = first_col + j;
      coefs[j] = 1.0;
    }

    crsmatrix->insertGlobalValues(global_row, nnz_per_row,
                                  &coefs[0], &indices[0]);
  }

  crsmatrix->fillComplete();

  std::cout << "In main(), pos 005" << std::endl;

  // Create Tpetra::MultiVectors for the lhs and rhs of the linear problem.
  Teuchos::RCP< Tpetra::MultiVector<double> > lhs = Teuchos::rcp( new Tpetra::MultiVector<double>(crsmatrix->getMap(), 1) );
  Teuchos::RCP< Tpetra::MultiVector<double> > rhs = Teuchos::rcp( new Tpetra::MultiVector<double>(crsmatrix->getMap(), 1) );
  rhs->putScalar( 1.0 );

  std::cout << "In main(), pos 006" << std::endl;

  // Create the linear problem with the original partitioning.
  Tpetra::LinearProblem< double
                       , Tpetra::CrsMatrix<double>::local_ordinal_type
                       , Tpetra::CrsMatrix<double>::global_ordinal_type
                       , Tpetra::CrsMatrix<double>::node_type
                       > problem( crsmatrix, lhs, rhs );

  std::cout << "In main(), pos 007" << std::endl;

  // Create the new linear problem and perform the balanced partitioning.
  // NOTE:  The balanced linear system will be in transformedProblem after fwd() is called.
  //        It is not necessary for the RCP to manage the transformed problem.
#if 0 // EEP
  Teuchos::RCP<Epetra_LinearProblem> transformedProblem = Teuchos::rcp( &((*LPTrans)( problem )), false );
  std::cout << "In main(), pos 008" << std::endl;
  LPTrans->fwd();
  std::cout << "In main(), pos 009" << std::endl;
#endif // EEP
  
  Tpetra::CrsMatrix<double>::local_ordinal_type graphrows1 = crsmatrix->getLocalNumRows();
  Tpetra::CrsMatrix<double>::local_ordinal_type bal_graph_rows = 9999; // transformedProblem->GetMatrix()->getLocalNumRows();
  Tpetra::CrsMatrix<double>::local_ordinal_type graphnnz1 = crsmatrix->getLocalNumEntries();
  Tpetra::CrsMatrix<double>::local_ordinal_type bal_graph_nnz = 9999; // transformedProblem->GetMatrix()->getLocalNumEntries();

  std::cout << "In main(), pos 010" << std::endl;

  for (int p(0); p < numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p
              << ": input matrix local rows = " << graphrows1
              << ", local NNZ = " << graphnnz1
              << std::endl;
    std::cout << "proc " << p
              << ": balanced matrix local rows = " << bal_graph_rows
              << ", local NNZ = " << bal_graph_nnz
              << std::endl;
  }

  std::cout << "Leaving runExample()" << std::endl;
}
#endif

int main(int argc, char** argv) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  runExample();
  MPI_Finalize();
#else
  std::cout << "ERROR:  This MUST be and MPI build with Tpetra enabled!" << std::endl;
#endif

  return(0);
}

