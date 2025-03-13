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
//#include <Epetra_LinearProblem.h>
//#include <EpetraExt_Isorropia_CrsGraph.h>
//#include <EpetraExt_LPTrans_From_GraphTrans.h>
//#include <Isorropia_Exception.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

//Declaration for helper-function that creates tpetra objects. These
//functions are implemented at the bottom of this file.
Teuchos::RCP< Tpetra::CrsMatrix<double/*, int, int*//*, std::int64_t, std::int64_t*//*, Tpetra::KokkosCompat::KokkosSerialWrapperNode*/> >
  create_tpetra_matrix(int numProcs, int localProc);

int main(int argc, char** argv) {
#ifdef HAVE_MPI
  //#if defined(HAVE_MPI) && defined(HAVE_TPETRA)

  int numProcs = 1;
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

  std::cout << "In main(), pos 000" << std::endl;

#if 0 // EEP
  
#ifdef HAVE_ISORROPIA_ZOLTAN

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

#endif // EEP

  std::cout << "In main(), pos 002" << std::endl;

#if 0 // EEP

  // Create the object to generate the balanced graph.
  Teuchos::RCP<EpetraExt::Isorropia_CrsGraph> Trans = 
    Teuchos::rcp( new EpetraExt::Isorropia_CrsGraph( paramlist ) );

  std::cout << "In main(), pos 003" << std::endl;
  
  // Create a linear problem transform interface.
  Teuchos::RCP<EpetraExt::LinearProblem_GraphTrans> LPTrans =
    Teuchos::rcp( new EpetraExt::LinearProblem_GraphTrans(
    *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(Trans.get())) ) );

#endif // EEP

  std::cout << "In main(), pos 004" << std::endl;

  // Create an Epetra_CrsMatrix object.
  Teuchos::RCP< Tpetra::CrsMatrix<double/*, int, int*//*, std::int64_t, std::int64_t*//*, Tpetra::KokkosCompat::KokkosSerialWrapperNode*/> > crsmatrix;
  try {
    crsmatrix = create_tpetra_matrix(numProcs, localProc);
  }
  catch(std::exception& e) {
    std::cout << "ERROR: create_tpetra_matrix threw exception '"
          << e.what() << "' on proc " << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  std::cout << "In main(), pos 005" << std::endl;

  // Create Epetra_MultiVectors for the lhs and rhs of the linear problem.
  Tpetra::MultiVector<double/*, int, int*//*, std::int64_t, std::int64_t*//*, Tpetra::KokkosCompat::KokkosSerialWrapperNode*/> lhs(crsmatrix->getMap(), 1);
  Tpetra::MultiVector<double/*, int, int*//*, std::int64_t, std::int64_t*//*, Tpetra::KokkosCompat::KokkosSerialWrapperNode*/> rhs(crsmatrix->getMap(), 1);
  rhs.putScalar( 1.0 );

  std::cout << "In main(), pos 006" << std::endl;

#if 0 // EEP
  
  // Create the linear problem with the original partitioning.
  Epetra_LinearProblem problem( &*crsmatrix, &lhs, &rhs );

  std::cout << "In main(), pos 007" << std::endl;

  // Create the new linear problem and perform the balanced partitioning.
  // NOTE:  The balanced linear system will be in tProblem after fwd() is called.
  //        It is not necessary for the RCP to manage the transformed problem.
  Teuchos::RCP<Epetra_LinearProblem> tProblem = Teuchos::rcp( &((*LPTrans)( problem )), false );
  std::cout << "In main(), pos 008" << std::endl;
  LPTrans->fwd();
  std::cout << "In main(), pos 009" << std::endl;

#endif // EEP
  
  int graphrows1 = crsmatrix->getLocalNumRows();
  //int bal_graph_rows = tProblem->GetMatrix()->getLocalNumRows();
  int graphnnz1 = crsmatrix->getLocalNumEntries();
  //int bal_graph_nnz = tProblem->GetMatrix()->getLocalNumEntries();

  std::cout << "In main(), pos 010" << std::endl;

  for(int p(0); p < numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p << ": input matrix local rows: " << graphrows1
       << ", local NNZ: " << graphnnz1 << std::endl;
    //std::cout << "proc " << p << ": balanced matrix local rows: "
    //   << bal_graph_rows << ", local NNZ: " << bal_graph_nnz << std::endl;
  }

  MPI_Finalize();

#else
  std::cout << "ERROR:  This MUST be and MPI build with Tpetra enabled!" << std::endl;
#endif

  return(0);
}

//Below are implementations of the helper-functions that create the
//poorly-balanced epetra objects for use in the above example program.

#ifdef HAVE_MPI
//#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

Teuchos::RCP< Tpetra::CrsMatrix<double/*, int, int*//*, std::int64_t, std::int64_t*//*, Tpetra::KokkosCompat::KokkosSerialWrapperNode*/> >
  create_tpetra_matrix(int numProcs, int localProc)
{
  if (localProc == 0) {
    std::cout << " creating Epetra_CrsMatrix with un-even distribution..."
            << std::endl;
  }
  //create an Epetra_CrsMatrix with rows spread un-evenly over
  //processors.
  Teuchos::MpiComm<int> comm(MPI_COMM_WORLD);
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
  Tpetra::Map</*int, int*//*std::int64_t, std::int64_t*//*, Tpetra::KokkosCompat::KokkosSerialWrapperNode*/> rowmap(global_num_rows, local_num_rows, 0, Teuchos::rcp(&comm));
  //create a matrix
  Teuchos::RCP< Tpetra::CrsMatrix<double/*, int, int*//*, std::int64_t, std::int64_t*//*, Tpetra::KokkosCompat::KokkosSerialWrapperNode*/> > matrix =
    Teuchos::rcp( new Tpetra::CrsMatrix<double/*, int, int*//*, std::int64_t, std::int64_t*//*, Tpetra::KokkosCompat::KokkosSerialWrapperNode*/>(Teuchos::rcp(&rowmap), nnz_per_row) );

#if 0 // EEP
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
        throw Isorropia::Exception("create_tpetra_matrix: error inserting matrix values.");
      }
    }
  }

  err = matrix->FillComplete();
  if (err != 0) {
    throw Isorropia::Exception("create_tpetra_matrix: error in matrix.FillComplete()");
  }
#endif // EEP
  return matrix;
}

#endif //HAVE_MPI && HAVE_EPETRA

