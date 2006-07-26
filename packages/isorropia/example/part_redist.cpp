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
//This file is a self-contained example of creating an Epetra_LinearProblem
//object, and using Isorropia to create a rebalanced copy of it.
//--------------------------------------------------------------------

//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>

//The Isorropia symbols being demonstrated are declared
//in these headers:
#include <Isorropia_Epetra.hpp>
#include <Isorropia_Redistributor.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

//Declaration for helper-function that creates epetra objects. This
//function is implemented at the bottom of this file.
#ifdef HAVE_EPETRA
Epetra_LinearProblem* create_epetra_problem(int numProcs,
                                            int localProc,
                                            int local_n);
#endif

int main(int argc, char** argv) {
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

  int p, numProcs = 1;
  int localProc = 0;

  //first, set up our MPI environment...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int local_n = 600;

  //Create a Epetra_LinearProblem object.

  Epetra_LinearProblem* linprob = 0;
  try {
    linprob = create_epetra_problem(numProcs, localProc, local_n);
  }
  catch(std::exception& exc) {
    std::cout << "linsys example: create_epetra_problem threw exception '"
          << exc.what() << "' on proc " << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  //We'll need a Teuchos::ParameterList object to pass to the
  //Isorropia::Partitioner class.
  Teuchos::ParameterList paramlist;

  // If Zoltan is available, the Zoltan package will be used for
  // the partitioning operation. By default, Isorropia selects Zoltan's
  // Hypergraph partitioner. If a method other than Hypergraph is
  // desired, it can be specified by first creating a parameter sublist
  // named "Zoltan", and then setting appropriate Zoltan parameters in
  // that sublist. A sublist is created like this:
  //      Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  //

  // If Zoltan is not available, a simple linear partitioner will be
  // used to partition such that the number of nonzeros is equal (or
  // close to equal) on each processor.


  Epetra_RowMatrix* rowmatrix = linprob->GetMatrix();
  Teuchos::RefCountPtr<const Epetra_RowMatrix> rowmat =
    Teuchos::rcp(rowmatrix, false);


  //Now create the partitioner object using an Isorropia factory-like
  //function...
  Teuchos::RefCountPtr<Isorropia::Partitioner> partitioner =
    Isorropia::Epetra::create_partitioner(rowmat, paramlist);


  //Next create a Redistributor object and use it to create balanced
  //copies of the objects in linprob.

  Isorropia::Redistributor rd(partitioner);

  Teuchos::RefCountPtr<Epetra_CrsMatrix> bal_matrix;
  Teuchos::RefCountPtr<Epetra_MultiVector> bal_x;
  Teuchos::RefCountPtr<Epetra_MultiVector> bal_b;

  //Use a try-catch block because Isorropia will throw an exception
  //if it encounters an error.

  if (localProc == 0) {
    std::cout << " calling Isorropia::Redistributor::redistribute..."
        << std::endl;
  }

  try {
    bal_matrix = rd.redistribute(*linprob->GetMatrix());
    bal_x = rd.redistribute(*linprob->GetLHS());
    bal_b = rd.redistribute(*linprob->GetRHS());
  }
  catch(std::exception& exc) {
    std::cout << "linsys example: Isorropia::Redistributor threw "
         << "exception '" << exc.what() << "' on proc "
         << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  Epetra_LinearProblem balanced_problem(bal_matrix.get(),
                                        bal_x.get(), bal_b.get());


  //Now query and print out information regarding the local sizes
  //of the original problem and the resulting balanced problem.

  int rows1 = linprob->GetMatrix()->NumMyRows();
  int bal_rows = balanced_problem.GetMatrix()->NumMyRows();
  int nnz1 = linprob->GetMatrix()->NumMyNonzeros();
  int bal_nnz = balanced_problem.GetMatrix()->NumMyNonzeros();

  for(p=0; p<numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p << ": original local rows: " << rows1
       << ", local NNZ: " << nnz1 << std::endl;
  }

  for(p=0; p<numProcs; ++p) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (p != localProc) continue;

    std::cout << "proc " << p << ": balanced prob local rows: "
       << bal_rows << ", local NNZ: " << bal_nnz << std::endl;
  }

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
  std::cout << "part_redist: must have both MPI and EPETRA. Make sure Trilinos "
    << "is configured with --enable-mpi and --enable-epetra." << std::endl;
#endif

  return(0);
}

//Below is the implementation of the helper-function that creates the
//poorly-balanced epetra objects for use in the above example program.

#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

Epetra_LinearProblem* create_epetra_problem(int numProcs,
                                            int localProc,
                                            int local_n)
{
  if (localProc == 0) {
    std::cout << " creating Epetra_CrsMatrix with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsMatrix with rows spread un-evenly over
  //processors.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int global_num_rows = numProcs*local_n;

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

  Epetra_Vector* x = new Epetra_Vector(rowmap);
  Epetra_Vector* b = new Epetra_Vector(rowmap);
  return(new Epetra_LinearProblem(matrix, x, b));
}

#endif //HAVE_MPI && HAVE_EPETRA

