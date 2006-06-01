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
//object, and using Isorropia to rebalance it.
//--------------------------------------------------------------------

//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>

//The Isorropia symbols being demonstrated are declared
//in these headers:
#include <Isorropia_Partitioner.hpp>
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

  int local_n = 5000;

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
  Teuchos::RefCountPtr<Teuchos::ParameterList> paramlist =
    Teuchos::rcp(new Teuchos::ParameterList);

  //(We won't pass any parameters in this example. But if we wanted
  // the Zoltan package to be used for the partitioning operation, we
  // could specify that like this:
  //       paramlist->set("Balancing package", "Zoltan");

  //We also need to pass a Epetra_CrsGraph object to the
  //Isorropia::Partitioner class. (This requirement may be removed in
  //the future...) So next we'll get the row-matrix from the linear-
  //problem, cast it to a crs-matrix, then obtain the crs-graph and
  //wrap that in a Teuchos::RefCountPtr.

  Epetra_RowMatrix* rowmatrix = linprob->GetMatrix();
  Epetra_CrsMatrix* crsmatrix = dynamic_cast<Epetra_CrsMatrix*>(rowmatrix);

  Teuchos::RefCountPtr<const Epetra_CrsGraph> graph =
    Teuchos::rcp(&(crsmatrix->Graph()), false);
  //The 'false' parameter specifies that the RefCountPtr should not take
  //ownership of the graph object.

  //Now create the partitioner object...
  Teuchos::RefCountPtr<Isorropia::Partitioner> partitioner =
    Teuchos::rcp(new Isorropia::Partitioner(graph, paramlist));

  partitioner->compute_partitioning();

  //Next create a Redistributor object and use it to create balanced
  //copies of the objects in linprob. By default, the result matrix is
  //balanced so that the number of nonzeros on each processor is equal
  //or close to equal.

  Isorropia::Redistributor rd(partitioner);

  Teuchos::RefCountPtr<Epetra_CrsMatrix> bal_matrix;
  Teuchos::RefCountPtr<Epetra_MultiVector> bal_x;
  Teuchos::RefCountPtr<Epetra_MultiVector> bal_b;

  //Use a try-catch block because Isorropia will throw an exception
  //if it encounters a fatal error.

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


  //Now simply query and print out information regarding the local sizes
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
  std::cout << "linsys: must have both MPI and EPETRA. Make sure Trilinos "
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
  int nnz_per_row = local_n/4+1;
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
  Epetra_CrsMatrix* matrix =
    new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  // Add  rows one-at-a-time
  double negOne = -1.0;
  double posTwo = 2.0;
  for (int i=0; i<local_n; i++) {
    int GlobalRow = matrix->GRID(i);
    int RowLess1 = GlobalRow - 1; int RowPlus1 = GlobalRow + 1;

    if (RowLess1!=-1) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess1);
    }
    if (RowPlus1!=global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus1);
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

