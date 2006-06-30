//@HEADER
// ************************************************************************
//
//               Isorropia: Partitioning and Load Balancing Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Alan Williams (william@sandia.gov)
//                 or Erik Boman    (egboman@sandia.gov)
//
// ************************************************************************
//@HEADER

//--------------------------------------------------------------------
//This file is a self-contained example of creating an Epetra_CrsGraph
//and Epetra_CrsMatrix object, and using Isorropia's create_balanced_copy
//function to rebalance them.
//--------------------------------------------------------------------

//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>

//The Isorropia user-interface functions being demonstrated are declared
//in Isorropia_Epetra.hpp.
#include <Isorropia_Epetra.hpp>

#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#endif

//Declarations for helper-functions that create epetra objects. These
//functions are implemented at the bottom of this file.
#ifdef HAVE_EPETRA
Teuchos::RefCountPtr<Epetra_CrsGraph>
  create_epetra_graph(int numProcs, int localProc);

Teuchos::RefCountPtr<Epetra_CrsMatrix>
  create_epetra_matrix(int numProcs, int localProc);
#endif

/** matrix_1 example program demonstrating Isorropia usage.
*/
int main(int argc, char** argv) {
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

  int p, numProcs = 1;
  int localProc = 0;

  //first, set up our MPI environment...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  //Create a Epetra_CrsGraph object. This graph will be the input to
  //the Isorropia rebalancing function...

  Teuchos::RefCountPtr<Epetra_CrsGraph> crsgraph;
  try {
    crsgraph = create_epetra_graph(numProcs, localProc);
  }
  catch(std::exception& exc) {
    std::cout << "matrix_1 example: create_epetra_graph threw exception '"
          << exc.what() << "' on proc " << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  //Now call Isorropia::create_balanced_copy to create a balanced
  //copy of crsgraph. By default, the result graph is balanced so that
  //the number of nonzeros on each processor is equal or close to equal.
  //Use a try-catch block because Isorropia will throw an exception
  //if it encounters a fatal error.

  if (localProc == 0) {
    std::cout << " calling Isorropia::create_balanced_copy..." << std::endl;
  }

  Teuchos::ParameterList paramlist;
  // If Zoltan is available, we'll specify that the Zoltan package be
  // used for the partitioning operation, by creating a parameter
  // sublist named "Zoltan".
  // In the sublist, we'll set parameters that we want sent to Zoltan.
#ifdef HAVE_ISORROPIA_ZOLTAN
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_METHOD", "GRAPH");
  sublist.set("PARMETIS_METHOD", "PARTKWAY");
#else
  // If Zoltan is not available, a simple linear partitioner will be
  // used to partition such that the number of nonzeros is equal (or
  // close to equal) on each processor. No parameter is necessary to
  // specify this.
#endif

  Teuchos::RefCountPtr<Epetra_CrsGraph> balanced_graph;
  try {
    balanced_graph =
      Isorropia::Epetra::create_balanced_copy(*crsgraph, paramlist);
  }
  catch(std::exception& exc) {
    std::cout << "matrix_1 example: Isorropia::create_balanced_copy threw "
         << "exception '" << exc.what() << "' on proc "
         << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  //Now query and print out information regarding the local sizes
  //of the input graph and the resulting balanced graph.

  int graphrows1 = crsgraph->NumMyRows();
  int bal_graph_rows = balanced_graph->NumMyRows();
  int graphnnz1 = crsgraph->NumMyNonzeros();
  int bal_graph_nnz = balanced_graph->NumMyNonzeros();

  for(p=0; p<numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p << ": input graph local rows: " << graphrows1
       << ", local NNZ: " << graphnnz1 << std::endl;
  }

  for(p=0; p<numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p << ": balanced graph local rows: "
       << bal_graph_rows << ", local NNZ: " << bal_graph_nnz << std::endl;
  }

  if (localProc == 0) {
    std::cout << std::endl;
  }

  //Next, do a similar exercise with a Epetra_CrsMatrix. Like the
  //Epetra_CrsGraph example above, we'll create a matrix to use as input,
  //and then have Isorropia::create_balanced_copy create a copy which is
  //balanced so that the number of nonzeros are equal on each processor.

  Teuchos::RefCountPtr<Epetra_CrsMatrix> crsmatrix;
  try {
    crsmatrix = create_epetra_matrix(numProcs, localProc);
  }
  catch(std::exception& exc) {
    std::cout << "matrix_1 example: create_epetra_matrix threw exception '"
          << exc.what() << "' on proc " << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

#ifdef HAVE_ISORROPIA_ZOLTAN
  //This time, we'll try hypergraph partitioning.
  sublist.set("LB_METHOD", "HYPERGRAPH");
#endif

  if (localProc == 0) {
    std::cout << " calling Isorropia::Epetra::create_balanced_copy..."
            << std::endl;
  }

  Teuchos::RefCountPtr<Epetra_CrsMatrix> balanced_matrix;
  try {
    balanced_matrix =
      Isorropia::Epetra::create_balanced_copy(*crsmatrix, paramlist);
  }
  catch(std::exception& exc) {
    std::cout << "matrix_1 example: Isorropia::create_balanced_copy(matrix)"
        << " threw exception '" << exc.what() << "' on proc "
         << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  int matrows1 = crsmatrix->NumMyRows();
  int bal_mat_rows = balanced_matrix->NumMyRows();
  int matnnz1 = crsmatrix->NumMyNonzeros();
  int bal_mat_nnz = balanced_matrix->NumMyNonzeros();

  for(p=0; p<numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p << ": input matrix local rows: " << matrows1
       << ", local NNZ: " << matnnz1 << std::endl;
  }

  for(p=0; p<numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p << ": balanced matrix local rows: "
        << bal_mat_rows << ", local NNZ: " << bal_mat_nnz << std::endl;
  }

  MPI_Finalize();

#else
  std::cout << "matrix_1: must have both MPI and EPETRA. Make sure Trilinos "
    << "is configured with --enable-mpi and --enable-epetra." << std::endl;
#endif

  return(0);
}

//Below are implementations of the helper-functions that create the
//poorly-balanced epetra objects for use in the above example program.

#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

Teuchos::RefCountPtr<Epetra_CrsMatrix>
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
  Teuchos::RefCountPtr<Epetra_CrsMatrix> matrix =
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

Teuchos::RefCountPtr<Epetra_CrsGraph>
  create_epetra_graph(int numProcs, int localProc)
{
  if (localProc == 0) {
    std::cout << " creating Epetra_CrsGraph with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsGraph with rows spread un-evenly over
  //processors.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int local_num_rows = 800;
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

  //create a graph
  Teuchos::RefCountPtr<Epetra_CrsGraph> graph =
    Teuchos::rcp(new Epetra_CrsGraph(Copy, rowmap, nnz_per_row));

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

    err = graph->InsertGlobalIndices(global_row, nnz_per_row,
                                         &indices[0]);
    if (err < 0) {
      throw Isorropia::Exception("create_epetra_graph: error inserting indices in graph");
    }
  }

  err = graph->FillComplete();
  if (err != 0) {
    throw Isorropia::Exception("create_epetra_graph: error in graph.FillComplete()");
  }

  return(graph);
}

#endif //HAVE_MPI && HAVE_EPETRA

