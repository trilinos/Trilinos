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
// This is an example of partitioning a hypergraph, while requesting
// partitions of unequal size.
//--------------------------------------------------------------------

//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>

//The Isorropia symbols being demonstrated are declared
//in these headers:
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>

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

//Declaration for helper-function that creates epetra rowmatrix objects. This
//function is implemented at the bottom of this file.
#ifdef HAVE_EPETRA
Teuchos::RCP<const Epetra_RowMatrix>
  create_epetra_rowmatrix(int numProcs,
                          int localProc,
                          int local_n);
#endif

int main(int argc, char** argv) {
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

  int numProcs = 1;
  int localProc = 0;

  //first, set up our MPI environment...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int local_n = 1200;

  //Create a Epetra_RowMatrix object.

  Teuchos::RCP<const Epetra_RowMatrix> rowmatrix;
  try {
    rowmatrix = create_epetra_rowmatrix(numProcs, localProc, local_n);
  }
  catch(std::exception& exc) {
    std::cout << "vert_weights example: create_epetra_rowmatrix threw"
         << " exception '" << exc.what() << "' on proc "
         << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  // Zoltan parameters

  Teuchos::ParameterList paramlist;

  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_METHOD", "HYPERGRAPH");

  // Now create the partitioner.  By default the partitioning occurs
  // in the constructor.  Include a flag to omit the partitioning.

  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner(
             rowmatrix,  // object to partition
             paramlist,  // parameters
             false));    // NO - do not do the partitioning in the constructor

  //Specify the proportion of rows to be assigned to each partition.  We'll
  //have the even-numbered partitions be twice the size of the odd-numbered
  //partitions.

  float *partSize = new float [numProcs];
  int *partGlobalId = new int [numProcs];

  for (int i=0; i < numProcs; i++){
    partGlobalId[i] = i;
    partSize[i] = (i % 2) ? 1.0 : 2.0;
  }

  partitioner->set_partition_sizes(numProcs, partGlobalId, NULL, partSize);

  delete [] partGlobalId;

  // NOW perform the partitioning

  partitioner->partition();

  //Next create a Redistributor object and use it to create a repartitioned
  //copy of the matrix.

  Isorropia::Epetra::Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsMatrix> bal_matrix;

  //Use a try-catch block because Isorropia will throw an exception
  //if it encounters an error.

  if (localProc == 0) {
    std::cout << " calling Isorropia::Epetra::Redistributor::redistribute..."
        << std::endl;
  }

  try {
    bal_matrix = rd.redistribute(*rowmatrix);
  }
  catch(std::exception& exc) {
    std::cout << "linsys example: Isorropia::Epetra::Redistributor threw "
         << "exception '" << exc.what() << "' on proc "
         << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  // How many rows does each process have?

  std::cout << "Partition " << localProc << ", Size " << partSize[localProc] << ", Number of rows " << bal_matrix->NumMyRows() << std::endl;


  delete [] partSize;

  MPI_Finalize();

#else
  std::cout << "part_sizes: must have both MPI and EPETRA. Make sure "
    << "Trilinos is configured with --enable-mpi and --enable-epetra."
     << std::endl;
#endif

  return(0);
}

//Below is the implementation of the helper-function that creates the
//epetra rowmatrix for use in the above example program.

#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

Teuchos::RCP<const Epetra_RowMatrix>
 create_epetra_rowmatrix(int numProcs,
                         int localProc,
                         int local_n)
{
  if (localProc == 0) {
    std::cout << " creating Epetra_CrsMatrix with even row-distribution..."
            << std::endl;
  }

  //create an Epetra_CrsMatrix with rows spread evenly over
  //processors.

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int global_num_rows = numProcs*local_n;

  Epetra_Map rowmap(global_num_rows, local_n, 0, comm);

  int nnz_per_row = 9;
  Teuchos::RCP<Epetra_CrsMatrix> matrix =
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row));

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

  return(matrix);
}

#endif //HAVE_MPI && HAVE_EPETRA

