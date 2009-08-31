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
//
// ************************************************************************
//@HEADER

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

#include <ispatest_utils.hpp>
#include <ispatest_epetra_utils.hpp>
#include <ispatest_lbeval_utils.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#endif

#ifdef HAVE_EPETRA
Epetra_CrsGraph* create_epetra_test_graph_1(int numProcs,
                                            int localProc,
                                            bool verbose);

Epetra_CrsMatrix* create_epetra_test_matrix_1(int numProcs,
                                              int localProc,
                                              bool verbose);

bool test_rebalance_epetra_linproblem(int numProcs, int localProc, bool verbose);
bool test_rebalance_epetra_crsmatrix(int numProcs, int localProc, bool verbose);
bool test_rebalance_epetra_graph(int numProcs, int localProc, bool verbose);
#endif

int main(int argc, char** argv) {
#ifdef HAVE_MPI

  bool verbose = false;
  int numProcs = 1;
  int localProc = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  try {
    verbose = ispatest::set_verbose(localProc, argc, argv);
  }
  catch(std::exception& exc) {
    std::cout << "err, setting verbosity: " << exc.what() << std::endl;
    std::cout << "End Result: TEST FAILED" << std::endl;
    MPI_Finalize();
    return(-1);
  }

  // Only print on proc 0
  if (localProc != 0) {
    verbose = false;
  }

  bool test_passed = true;
#ifdef HAVE_EPETRA
  test_passed = test_passed &&
    test_rebalance_epetra_graph(numProcs, localProc, verbose);

  test_passed = test_passed &&
    test_rebalance_epetra_crsmatrix(numProcs, localProc, verbose);

  test_passed = test_passed &&
    test_rebalance_epetra_linproblem(numProcs, localProc, verbose);
#else
  std::cout << "rebalance_zoltan main: currently can only test "
         << "rebalancing with Epetra enabled." << std::endl;
  test_passed = false;
#endif

  if (test_passed && verbose) {
    std::cout << "rebalance_zoltan main: tests passed."<<std::endl;
  }

#else
  std::cout << "rebalance_zoltan: don't have MPI, can't run test."
            << std::endl;
#endif

  MPI_Finalize();

  return(0);
}

#ifdef HAVE_EPETRA
//-------------------------------------------------------------------
bool test_rebalance_epetra_crsmatrix(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif

  double bal1, bal2, cutl1, cutl2, cutn1, cutn2;
  double myShare = 1.0 / numProcs;

  Epetra_CrsMatrix* input_matrix =
    create_epetra_test_matrix_1(numProcs, localProc, verbose);

  Isorropia::Epetra::CostDescriber costs;   // default costs

  ispatest::compute_hypergraph_metrics(*input_matrix, costs, myShare, bal1, cutn1, cutl1);

  //We'll use Zoltan for the rebalancing:

  Teuchos::ParameterList paramlist;
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_METHOD", "HYPERGRAPH");

  Epetra_CrsMatrix *balanced_matrix;
  try {
    if (verbose) {
      std::cout << " calling Isorropia::createBalancedCopy(Epetra_CrsMatrix)..."
                << std::endl;
    }

    balanced_matrix =
      Isorropia::Epetra::createBalancedCopy(*input_matrix, paramlist);
  }
  catch(std::exception& exc) {
    std::cout << "caught exception: " << exc.what() << std::endl;
    return(false);
  }
  ispatest::compute_hypergraph_metrics(*balanced_matrix, costs, myShare, bal2, cutn2, cutl2);

  if (verbose) {
    std::cout << "before balancing hypergraph" << std::endl;
    std::cout << "balance " << bal1 << " cutn " << cutn1 << " cutl " << cutl1 << std::endl;
    std::cout << "after balancing hypergraph" << std::endl;
    std::cout << "balance " << bal2 << " cutn " << cutn2 << " cutl " << cutl2 << std::endl;
  }

  test_passed = (cutl2 <= cutl1);

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }

  return(test_passed);
}

//-------------------------------------------------------------------
bool test_rebalance_epetra_linproblem(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif
  double bal1, bal2, cutl1, cutl2, cutn1, cutn2;
  double myShare = 1.0 / numProcs;

  Isorropia::Epetra::CostDescriber costs;   // default costs

  Epetra_CrsMatrix* input_matrix =
    create_epetra_test_matrix_1(numProcs, localProc, verbose);

  ispatest::compute_hypergraph_metrics(*input_matrix, costs, myShare, bal1, cutn1, cutl1);

  Epetra_Vector* x = new Epetra_Vector(input_matrix->RowMap());
  Epetra_Vector* b = new Epetra_Vector(input_matrix->RowMap());
  Epetra_LinearProblem problem(input_matrix, x, b);

  Teuchos::ParameterList paramlist;
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  // No parameters will test default Zoltan method (hypergraph)
  //sublist.set("LB_METHOD", "GRAPH");
  //sublist.set("PARMETIS_METHOD", "PARTKWAY");

  //Wrap a RefCountPtr around the matrix graph, and specify 'false', meaning
  //that the RefCountPtr will not take ownership of the graph (will not
  //delete it).
  Teuchos::RefCountPtr<const Epetra_CrsGraph> graph =
    Teuchos::rcp( &(input_matrix->Graph()), false);

  Teuchos::RefCountPtr<Isorropia::Epetra::Partitioner> partitioner =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner(graph, paramlist));

  Isorropia::Epetra::Redistributor rd(partitioner);

  Teuchos::RefCountPtr<Epetra_CrsMatrix> bal_matrix =
    rd.redistribute(*(problem.GetMatrix()));

  ispatest::compute_hypergraph_metrics(*bal_matrix, costs, myShare, bal2, cutn2, cutl2);

  Teuchos::RefCountPtr<Epetra_MultiVector> bal_x =
    rd.redistribute(*(problem.GetLHS()));

  Teuchos::RefCountPtr<Epetra_MultiVector> bal_b =
    rd.redistribute(*(problem.GetRHS()));

  Teuchos::RefCountPtr<Epetra_LinearProblem> balanced_problem =
    Teuchos::rcp(new Epetra_LinearProblem(bal_matrix.get(),
					  bal_x.get(), bal_b.get()));

  if (verbose) {
    std::cout << "before balancing hypergraph" << std::endl;
    std::cout << "balance " << bal1 << " cutn " << cutn1 << " cutl " << cutl1 << std::endl;
    std::cout << "after balancing hypergraph" << std::endl;
    std::cout << "balance " << bal2 << " cutn " << cutn2 << " cutl " << cutl2 << std::endl;
  }

  test_passed = (cutl2 <= cutl1);

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }

  delete input_matrix;
  delete x;
  delete b;

  return(test_passed);
}

//-------------------------------------------------------------------
bool test_rebalance_epetra_graph(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif
  Isorropia::Epetra::CostDescriber costs;   // default costs
  double bal1, bal2, cutl1, cutl2, cutn1, cutn2;
  double myShare = 1.0 / numProcs;

  Epetra_CrsGraph* input_graph =
    create_epetra_test_graph_1(numProcs, localProc, verbose);

  //We'll specify that the package Zoltan should be used to perform the
  //repartitioning, and call the Isorropia::createBalancedCopy function.
  //Default behavior should be to balance the graph so that the number of
  //nonzeros on each processor is roughly equal. i.e., by default, weights
  //for each row are assumed to be the number of nonzeros in that row.

  Teuchos::ParameterList paramlist;
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  // No parameters will test default Zoltan method (hypergraph)
  //sublist.set("LB_METHOD", "GRAPH");
  //sublist.set("PARMETIS_METHOD", "PARTKWAY");

  Epetra_CrsGraph *balanced_graph;
  try {
    if (verbose) {
      std::cout << " calling Isorropia::createBalancedCopy(Epetra_CrsGraph)..."
                << std::endl;
    }

    balanced_graph =
      Isorropia::Epetra::createBalancedCopy(*input_graph, paramlist);
  }
  catch(std::exception& exc) {
    std::cout << "caught exception: " << exc.what() << std::endl;
    delete input_graph;
    return(false);
  }

  ispatest::compute_hypergraph_metrics(*input_graph, costs, myShare, bal1, cutn1, cutl1);
  ispatest::compute_hypergraph_metrics(*balanced_graph, costs, myShare, bal2, cutn2, cutl2);

  if (verbose) {
    std::cout << "before balancing hypergraph" << std::endl;
    std::cout << "balance " << bal1 << " cutn " << cutn1 << " cutl " << cutl1 << std::endl;
    std::cout << "after balancing hypergraph" << std::endl;
    std::cout << "balance " << bal2 << " cutn " << cutn2 << " cutl " << cutl2 << std::endl;
  }

  test_passed = (cutl2 <= cutl1);

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }
  return(test_passed);
}

Epetra_CrsMatrix* create_epetra_test_matrix_1(int numProcs,
                                              int localProc,
                                              bool verbose)
{
  if (verbose) {
    std::cout << " creating Epetra_CrsMatrix with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsMatrix with rows spread un-evenly over
  //processors.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int local_num_rows = 360;
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
  Epetra_CrsMatrix* input_matrix =
    new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  int err = ispatest::fill_matrix(*input_matrix, nnz_per_row, verbose);
  if (err != 0) {
    delete input_matrix;
    return(0);
  }

  return(input_matrix);
}

Epetra_CrsGraph* create_epetra_test_graph_1(int numProcs,
                                              int localProc,
                                              bool verbose)
{
  if (verbose) {
    std::cout << " creating Epetra_CrsGraph with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsGraph with rows spread un-evenly over
  //processors.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int local_num_rows = 360;
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
  Epetra_CrsGraph* input_graph =
    new Epetra_CrsGraph(Copy, rowmap, nnz_per_row);

  int err = ispatest::fill_graph(*input_graph, nnz_per_row, verbose);
  if (err != 0) {
    delete input_graph;
    return(0);
  }

  return(input_graph);
}

#endif //HAVE_EPETRA

