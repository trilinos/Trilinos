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
//This file is a self-contained example of creating an Epetra_RowMatrix
//object, and using Isorropia to create a rebalanced copy of it.
//Furthermore graph-edge weights are used to influence the repartitioning.
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
#include <Epetra_FECrsMatrix.h>
#endif

int main(int argc, char** argv) {
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

  int p, numProcs = 1;
  int localProc = 0;

  //first, set up our MPI environment...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

// This program can only run on 3 processors, to make things
// work out easy.
//
  if (numProcs != 3) {
    std::cout << "num-procs="<<numProcs<<". This program can only "
      << "run on 3 procs. Exiting."<<std::endl;
    MPI_Finalize();
    return(0);
  }

//Consider the following mesh of 4 2-D quad elements:
//
//  *-------*-------*
// 8|      7|      6|
//  |  E2   |  E3   |
//  *-------*-------*
// 3|      2|      5|
//  |  E0   |  E1   |
//  *-------*-------*
// 0       1       4
//
// Node-ids are to the lower-left of each node (*).
//
// Mimicing a finite-element application, we will say that
// each node has 1 scalar degree-of-freedom, and assemble
// a matrix which would have 9 global rows and columns.
//
// Each processor will have 3 rows. We'll set up a strange
// initial map, where nodes are distributed as follows:
//
// proc 0: nodes 0,3,8,
// proc 1: nodes 1,2,7
// proc 2: nodes 4,5,6.
//
// After we assemble our matrix, we'll create another matrix
// and populate it with graph edge weights such that the
// partitioner repartitions the problem so that nodes are
// laid out as follows:
//
// proc 0: nodes 0, 1, 4
// proc 1: nodes 3, 2, 5
// proc 2: nodes 8, 7, 6
//

  int nodesPerElem = 4;
  int global_n = 9;

  //First, set up the initial map:

  std::vector<int> mynodes(3);
  if (localProc == 0) {
    mynodes[0] = 0; mynodes[1] = 3; mynodes[2] = 8;
  }
  if (localProc == 1) {
    mynodes[0] = 1; mynodes[1] = 2; mynodes[2] = 7;
  }
  if (localProc == 2) {
    mynodes[0] = 4; mynodes[1] = 5; mynodes[2] = 6;
  }

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map origmap(global_n, 3, &mynodes[0], 0, comm);

  Teuchos::RCP<Epetra_FECrsMatrix> matrix =
    Teuchos::rcp(new Epetra_FECrsMatrix(Copy, origmap, 0));

  //We'll assemble elements E0 and E1 on proc 0,
  //               element E2 or proc 1,
  //               element E3 on proc 2.

  std::vector<int> indices(nodesPerElem);
  std::vector<double> coefs(nodesPerElem*nodesPerElem,2.0);

  if (localProc == 0) {
    //element E0:
    indices[0] = 0; indices[1] = 1; indices[2] = 2; indices[3] = 3;
    matrix->InsertGlobalValues(nodesPerElem, &indices[0], &coefs[0]);

    //element E1:
    indices[0] = 1; indices[1] = 4; indices[2] = 5; indices[3] = 2;
    matrix->InsertGlobalValues(nodesPerElem, &indices[0], &coefs[0]);
  }
  else if (localProc == 1) {
    //element E2:
    indices[0] = 3; indices[1] = 2; indices[2] = 7; indices[3] = 8;
    matrix->InsertGlobalValues(nodesPerElem, &indices[0], &coefs[0]);
  }
  else { //localProc==2
    //element E3:
    indices[0] = 2; indices[1] = 5; indices[2] = 6; indices[3] = 7;
    matrix->InsertGlobalValues(nodesPerElem, &indices[0], &coefs[0]);
  }

  int err = matrix->GlobalAssemble();
  if (err != 0) {
    std::cout << "err="<<err<<" returned from matrix->GlobalAssemble()"
      << std::endl;
  }

//  std::cout << "matrix: " << std::endl;
//  std::cout << *matrix << std::endl;

  //We'll need a Teuchos::ParameterList object to pass to the
  //Isorropia::Epetra::Partitioner class.
  Teuchos::ParameterList paramlist;

#ifdef HAVE_ISORROPIA_ZOLTAN
  // If Zoltan is available, we'll specify that the Zoltan package be
  // used for the partitioning operation, by creating a parameter
  // sublist named "Zoltan".
  // In the sublist, we'll set parameters that we want sent to Zoltan.

  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_METHOD", "GRAPH");
  sublist.set("PARMETIS_METHOD", "PARTKWAY");
  //sublist.set("PARMETIS_ITR", "1000");
#else
  // If Zoltan is not available, a simple linear partitioner will be
  // used to partition such that the number of nonzeros is equal (or
  // close to equal) on each processor. No parameter is necessary to
  // specify this.
#endif


  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs =
    Teuchos::rcp(new Isorropia::Epetra::CostDescriber);

  //Next create a matrix which is a copy of the matrix we just
  //assembled, but we'll replace the values with graph edge weights.

  Teuchos::RCP<Epetra_FECrsMatrix> ge_weights =
    Teuchos::rcp(new Epetra_FECrsMatrix(*matrix));

  Teuchos::RCP<Epetra_CrsMatrix> crs_ge_weights;
  crs_ge_weights = ge_weights;

  //Fill the matrix with a "default" weight of 1.0.
  crs_ge_weights->PutScalar(1.0);

  //Now we'll put a "large" weight on edges that connect nodes
  //0 and 1, 1 and 4,
  //3 and 2, 2 and 5,
  //8 and 7, 7 and 6.

  double weight = 500.0;

  if (localProc == 0) {
    //row 0, edge 1
    indices[0] = 1;
    coefs[0] = weight;
    crs_ge_weights->ReplaceGlobalValues(0, 1, &coefs[0], &indices[0]);

    //row 3, edge 2
    indices[0] = 2;
    coefs[0] = weight;
    crs_ge_weights->ReplaceGlobalValues(3, 1, &coefs[0], &indices[0]);

    //row 8, edge 7
    indices[0] = 7;
    coefs[0] = weight;
    crs_ge_weights->ReplaceGlobalValues(8, 1, &coefs[0], &indices[0]);
  }

  if (localProc == 1) {
    //row 1, edges 0 and 4
    indices[0] = 0; indices[1] = 4;
    coefs[0] = weight; coefs[1] = weight;
    crs_ge_weights->ReplaceGlobalValues(1, 2, &coefs[0], &indices[0]);

    //row 2, edges 3 and 5
    indices[0] = 3; indices[1] = 5;
    coefs[0] = weight; coefs[1] = weight;
    crs_ge_weights->ReplaceGlobalValues(2, 2, &coefs[0], &indices[0]);

    //row 7, edges 6 and 8
    indices[0] = 6; indices[1] = 8;
    coefs[0] = weight;
    crs_ge_weights->ReplaceGlobalValues(7, 2, &coefs[0], &indices[0]);
  }

  if (localProc == 2) {
    //row 4, edge 1
    indices[0] = 1;
    coefs[0] = weight;
    crs_ge_weights->ReplaceGlobalValues(4, 1, &coefs[0], &indices[0]);

    //row 5, edge 2
    indices[0] = 2;
    coefs[0] = weight;
    crs_ge_weights->ReplaceGlobalValues(5, 1, &coefs[0], &indices[0]);

    //row 6, edge 7
    indices[0] = 7;
    coefs[0] = weight;
    crs_ge_weights->ReplaceGlobalValues(6, 1, &coefs[0], &indices[0]);
  }

// std::cout << "crs_ge_weights: " << std::endl
//       << *crs_ge_weights << std::endl;

  //Now give the graph edge weights to the CostDescriber:
  costs->setGraphEdgeWeights(crs_ge_weights);

  Teuchos::RCP<const Epetra_RowMatrix> rowmatrix;
  rowmatrix = matrix;

  //Now create the partitioner object using an Isorropia factory-like
  //function...
  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
    Isorropia::Epetra::create_partitioner(rowmatrix, costs, paramlist);


  //Next create a Redistributor object and use it to create a
  //repartitioned copy of the matrix

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


  //Now print out the original map and then the map of the balanced
  //matrix, to show that the rows have been migrated to different
  //processors.

  std::cout << "origmap: " << std::endl;
  std::cout << origmap << std::endl;

  comm.Barrier();

  std::cout << std::endl;

  std::cout << "bal_matrix->RowMap(): " << std::endl;
  std::cout << bal_matrix->RowMap() << std::endl;

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

