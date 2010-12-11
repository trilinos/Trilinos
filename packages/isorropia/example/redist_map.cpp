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

//--------------------------------------------------------------------
// This file is a self-contained example of creating an Epetra_LinearProblem
// object, and using Isorropia to create a rebalanced copy of it by providing a
// target map (not the partitioner). This example also demonstrates that it
// is possible to use different maps for the vectors and the matrix of a 
// linear problem and redistribute them correctly.
//--------------------------------------------------------------------

//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>

//The Isorropia symbols being demonstrated are declared
//in these headers:
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

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
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

#include "ispatest_lbeval_utils.hpp"

//Declaration for helper-function that creates epetra objects. This
//function is implemented at the bottom of this file.
#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
Epetra_LinearProblem* create_epetra_problem(int numProcs,
                                            int localProc,
                                            int local_n,
                                            Epetra_MpiComm comm);
#endif
#endif

int main(int argc, char** argv) {
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

  int numProcs = 1;
  int localProc = 0;

  //first, set up our MPI environment...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  int local_n = 600;

  //Create a Epetra_LinearProblem object.

  Epetra_LinearProblem* linprob = 0;
  try {
    linprob = create_epetra_problem(numProcs, localProc, local_n, comm);
  }
  catch(std::exception& exc) {
    std::cout << "redist_map example: create_epetra_problem threw exception '"
          << exc.what() << "' on proc " << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  // Create a target map to redistribute to (localize the matrix to proc 0)
  int global_n = local_n * numProcs;
  int target_local_n ;
  if (localProc == 0)
      target_local_n = global_n;
  else
      target_local_n = 0;

  Epetra_Map target_map(global_n, target_local_n, 0, comm);

  //Create a Redistributor object and use it to create redistributed
  //copies of the objects in linprob.
  Isorropia::Epetra::Redistributor rd(Teuchos::rcpFromRef(target_map));

  Teuchos::RCP<Epetra_CrsMatrix> bal_matrix;
  Teuchos::RCP<Epetra_MultiVector> bal_x;
  Teuchos::RCP<Epetra_MultiVector> bal_b;

  if (localProc == 0) {
    std::cout << " calling Isorropia::Epetra::Redistributor::redistribute..."
        << std::endl;
  }

  try {
    bal_matrix = rd.redistribute(*linprob->GetMatrix());
    bal_x = rd.redistribute(*linprob->GetLHS());
    bal_b = rd.redistribute(*linprob->GetRHS());
  }
  catch(std::exception& exc) {
    std::cout << "redist_map example: Isorropia::Epetra::Redistributor threw "
         << "exception '" << exc.what() << "' on proc "
         << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  Epetra_LinearProblem balanced_problem(bal_matrix.get(),
                                        bal_x.get(), bal_b.get());

  // Results
  double bal0, bal1, cutn0, cutn1, cutl0, cutl1;
  Isorropia::Epetra::CostDescriber default_costs;

  // Balance and cut quality before partitioning
  double goalWeight = 1.0 / (double)numProcs;
  ispatest::compute_hypergraph_metrics(*(linprob->GetMatrix()), default_costs,
        goalWeight, bal0, cutn0, cutl0);

  // Balance and cut quality after partitioning
  ispatest::compute_hypergraph_metrics(*bal_matrix, default_costs, goalWeight,
                     bal1, cutn1, cutl1);

  if (localProc == 1){
    std::cout << "Before redistribute: ";
    std::cout << "Balance " << bal0 << " cutN " << cutn0 << " cutL " << cutl0;
    std::cout << std::endl;

    std::cout << "After redistribute (localize to proc 0):  ";
    std::cout << "Balance " << bal1 << " cutN " << cutn1 << " cutL " << cutl1;
    std::cout << std::endl;
  }

  // Check the result of the localize
  assert (bal1 == numProcs);
  assert (cutn1 == 0);
  assert (cutl1 == 0);

  //Finally, delete the pointer objects that we asked to be created.
  delete linprob->GetMatrix();
  delete linprob->GetLHS();
  delete linprob->GetRHS();
  delete linprob;

  if (localProc == 0) {
    std::cout << std::endl;
  }

  MPI_Finalize();

#else
  std::cout << "redist_map: must have both MPI and EPETRA. Make sure Trilinos "
    << "is configured with --enable-mpi and --enable-epetra." << std::endl;
#endif

  return(0);
}

//Below is the implementation of the helper-function that creates the
//poorly-balanced epetra objects for use in the above example program.

#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

Epetra_LinearProblem* create_epetra_problem(int numProcs,
                                            int localProc,
                                            int local_n,
                                            Epetra_MpiComm comm)
{
  if (localProc == 0) {
    std::cout << " creating Epetra_CrsMatrix with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsMatrix with rows spread un-evenly over
  //processors.
  //Epetra_MpiComm comm(MPI_COMM_WORLD);
  int global_num_rows = numProcs*local_n;

  // Make sure the vectors use a different distribution than the matrix
  Epetra_Map vector_rowmap(global_num_rows, local_n, 0, comm);
  Epetra_Vector* x = new Epetra_Vector(vector_rowmap);
  Epetra_Vector* b = new Epetra_Vector(vector_rowmap);

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;

  int adjustment = local_n/2;

  //adjust local_n so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_n -= adjustment;
  }
  else {
    local_n += adjustment;
  }

  //if numProcs is not an even number, undo the local_n adjustment
  //on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_n -= adjustment;
    }
  }

  //now we're ready to create a row-map.
  Epetra_Map rowmap(global_num_rows, local_n, 0, comm);

  //create a matrix
  int nnz_per_row = 9;
  Epetra_CrsMatrix* matrix =
    new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  // Add  rows one-at-a-time
  double negOne = -1.0;
  double posTwo = 4.0;
  for (int i=0; i<local_n; i++) {
    int GlobalRow = matrix->GRID(i);
    int RowLess1 = GlobalRow - 1;
    int RowPlus1 = GlobalRow + 1;
    int RowLess2 = GlobalRow - 2;
    int RowPlus2 = GlobalRow + 2;
    int RowLess3 = GlobalRow - 3;
    int RowPlus3 = GlobalRow + 3;
    int RowLess4 = GlobalRow - 4;
    int RowPlus4 = GlobalRow + 4;

    if (RowLess4>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess4);
    }
    if (RowLess3>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess3);
    }
    if (RowLess2>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess2);
    }
    if (RowLess1>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess1);
    }
    if (RowPlus1<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus1);
    }
    if (RowPlus2<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus2);
    }
    if (RowPlus3<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus3);
    }
    if (RowPlus4<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus4);
    }

    matrix->InsertGlobalValues(GlobalRow, 1, &posTwo, &GlobalRow);
  }

  int err = matrix->FillComplete();
  if (err != 0) {
    throw Isorropia::Exception("create_epetra_matrix: error in matrix.FillComplete()");
  }

  return(new Epetra_LinearProblem(matrix, x, b));
}

#endif //HAVE_MPI && HAVE_EPETRA
