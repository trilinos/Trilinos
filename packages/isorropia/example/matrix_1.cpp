//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

//--------------------------------------------------------------------
//This file is a self-contained example of creating an Epetra_CrsGraph
//and Epetra_CrsMatrix object, and using Isorropia's createBalancedCopy
//function to rebalance them.
//--------------------------------------------------------------------

//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

//The Isorropia user-interface functions being demonstrated are declared
//in Isorropia_Epetra.hpp.
#include <Isorropia_Epetra.hpp>

#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#endif

#ifdef HAVE_ISPATEST
#include "ispatest_lbeval_utils.hpp"
#endif

//Declarations for helper-functions that create epetra objects. These
//functions are implemented at the bottom of this file.
#ifdef HAVE_EPETRA
Teuchos::RCP<Epetra_CrsGraph>
  create_epetra_graph(int numProcs, int localProc);

Teuchos::RCP<Epetra_CrsMatrix>
  create_epetra_matrix(int numProcs, int localProc);
#endif

/** matrix_1 example program demonstrating simple Isorropia usage.
*/
int main(int argc, char** argv) {
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

  int numProcs = 1;
  int localProc = 0;

  //first, set up our MPI environment...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  //Create a Epetra_CrsGraph object. This graph will be the input to
  //the Isorropia rebalancing function...

  Teuchos::RCP<Epetra_CrsGraph> crsgraph;
  try {
    crsgraph = create_epetra_graph(numProcs, localProc);
  }
  catch(std::exception& exc) {
    std::cout << "matrix_1 example: create_epetra_graph threw exception '"
          << exc.what() << "' on proc " << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  //Now call Isorropia::createBalancedCopyto create a balanced
  //copy of crsgraph.

  if (localProc == 0) {
    std::cout << "Hypergraph partitioning" << std::endl;
  }

  Teuchos::ParameterList paramlist;
  //No parameters. By default, Isorropia will use Zoltan hypergraph
  //partitioning, treating the graph columns as hyperedges and the
  //graph rows as vertices.

  Epetra_CrsGraph *balanced_graph = NULL;
  try {
    balanced_graph =
      Isorropia::Epetra::createBalancedCopy(*crsgraph, paramlist);

  }
  catch(std::exception& exc) {
    std::cout << "matrix_1 example: Isorropia::createBalancedCopy threw "
         << "exception '" << exc.what() << "' on proc "
         << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  // Results

#ifdef HAVE_ISPATEST
  Isorropia::Epetra::CostDescriber emptyCostObject;
  double bal0, bal1, cutn0, cutn1, cutl0, cutl1;

  // Balance and cut quality before partitioning

  double goalWeight = 1.0 / (double)numProcs;
  ispatest::compute_hypergraph_metrics(*crsgraph, emptyCostObject, goalWeight,
                     bal0, cutn0, cutl0);

  // Balance and cut quality after partitioning

  ispatest::compute_hypergraph_metrics(*balanced_graph, emptyCostObject, goalWeight,
                     bal1, cutn1, cutl1);

  if (localProc == 0){
    std::cout << "Before partitioning hypergraph: ";
    std::cout << "Balance " << bal0 << " cutN " << cutn0 << " cutL " << cutl0;
    std::cout << std::endl;

    std::cout << "After partitioning hypergraph:  ";
    std::cout << "Balance " << bal1 << " cutN " << cutn1 << " cutL " << cutl1;
    std::cout << std::endl;
    std::cout << std::endl;
  }
#else
  (void) balanced_graph; // unused
#endif


  //Next, do a similar exercise with a Epetra_CrsMatrix. Like the
  //Epetra_CrsGraph example above, we'll create a matrix to use as input,
  //and then have Isorropia::createBalancedCopy create a copy which is
  //balanced so that the number of nonzeros are equal on each processor.

  Teuchos::RCP<Epetra_CrsMatrix> crsmatrix;
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
  //This time, we'll try graph partitioning..
  //Create a parameter sublist for Zoltan parameters.
  paramlist.set("PARTITIONING METHOD", "GRAPH");

  //We could also call ParMetis, if Zoltan was built with ParMetis.
  //Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  //sublist.set("GRAPH_PACKAGE", "PARMETIS");
  //sublist.set("PARMETIS_METHOD", "PARTKWAY");

#else
  //This should never happen, since Zoltan is now required by Isorropia.
#endif


  if (localProc == 0) {
    std::cout << "Specifying GRAPH partitioning." << std::endl;
  }

  Epetra_CrsMatrix *balanced_matrix;
  try {
    balanced_matrix =
      Isorropia::Epetra::createBalancedCopy(*crsmatrix, paramlist);
  }
  catch(std::exception& exc) {
    std::cout << "matrix_1 example: Isorropia::createBalancedCopy(matrix)"
        << " threw exception '" << exc.what() << "' on proc "
         << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  // Results

#ifdef HAVE_ISPATEST
  double cutWgt0, cutWgt1;
  int numCuts0, numCuts1;

  // Balance and cut quality before partitioning

  ispatest::compute_graph_metrics(*crsmatrix, emptyCostObject, goalWeight,
                     bal0, numCuts0, cutWgt0, cutn0, cutl0);

  // Balance and cut quality after partitioning

  ispatest::compute_graph_metrics(*balanced_matrix, emptyCostObject, goalWeight,
                     bal1, numCuts1, cutWgt1, cutn1, cutl1);

  if (localProc == 0){
    std::cout << "Before partitioning graph: Number of cuts " << numCuts0 << " Cut weight " << cutWgt0 << std::endl;
    std::cout << "                     Balance " << bal0 << " cutN " << cutn0 << " cutL " << cutl0;
    std::cout << std::endl;

    std::cout << "After partitioning graph:  Number of cuts " << numCuts1 << " Cut weight " << cutWgt1 << std::endl;
    std::cout << "                     Balance " << bal1 << " cutN " << cutn1 << " cutL " << cutl1;
    std::cout << std::endl;
  }
#else
  (void) balanced_matrix; // unused
#endif

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

Teuchos::RCP<Epetra_CrsMatrix>
  create_epetra_matrix(int numProcs, int localProc)
{
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

    err = matrix->InsertGlobalValues(global_row, nnz_per_row,
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

Teuchos::RCP<Epetra_CrsGraph>
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
  Teuchos::RCP<Epetra_CrsGraph> graph =
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

