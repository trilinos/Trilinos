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




#ifdef HAVE_MPI



// Need to add a function to utilities to compute cutl or some
// balance measure so we can compare quality with Zoltan w/o Isorropia
//
// Read in matrix market file like MapColoring/cxx_main.cpp in EpetraExt
// Partition with Zoltan hg
// Calculate cuts
//
// For out of source builds, Makefile needs to copy simple.mtx file from
// source to build directory.
//

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_Redistributor.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <ispatest_utils.hpp>
#include <ispatest_epetra_utils.hpp>

#include <Epetra_SerialDenseVector.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_BlockMapIn.h>
#endif

static void show_matrix(const char *txt, Epetra_CrsMatrix *matrix, const Epetra_Comm &Comm);

int main(int argc, char** argv) {

  int rc=0;  
#ifndef HAVE_MPI
  std::cout << "test_simple: don't have MPI, can't run test."
            << std::endl;
  return 1;
#endif
#ifndef HAVE_EPETRA
  std::cout << "test_simple main: currently can only test "
         << "with Epetra enabled." << std::endl;
  return 1;
#endif

  bool verbose = false;
  int numProcs = 1;
  int localProc = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  const Epetra_MpiComm Comm(MPI_COMM_WORLD);

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

  Epetra_CrsMatrix *matrix_ptr;
  rc = EpetraExt::MatrixMarketFileToCrsMatrix("simple.mtx", Comm, matrix_ptr);
  if (rc < 0){
    std::cout << "error reading in simple.mtx" << std::endl;
    exit(1);
  }

  matrix_ptr->FillComplete();

  show_matrix("Before load balancing", matrix_ptr, Comm);

  // By convention, the rows of this matrix are hyperedges, and we
  // wish to balance the vertices, represented by the columns.

  // Create hyperedge weights, equal to the number of vertices in
  // the edge (or equivalently the number of non-zeroes in the row.
  // 
  // For vertex weights, we'll take the default of equal weight
  // vertices.  TODO: the CostDescriber needs an interface to specify
  // vertex weights for hypergraphs.  The existing interface just
  // asks for an Epetra_Vector of weights.  For graphs this is OK, but
  // for hypergraphs we need to say which vertices those weights are for.

  Teuchos::ParameterList params;
  Teuchos::ParameterList &sublist = params.sublist("Zoltan");
  sublist.set("LB_METHOD", "HYPERGRAPH");
  sublist.set("PHG_CUT_OBJECTIVE", "CONNECTIVITY");  // "cutl"
  sublist.set("EDGE_WEIGHT_DIM", "1");  

  const Epetra_Map &myRowMap = matrix_ptr->RowMap();
  int numMyRows = matrix_ptr->NumMyRows();
  double *ewgts = new double[numMyRows];
  for (int i=0; i<numMyRows; i++){
    //ewgts[i] = myRowMap.ElementSize(i);
    ewgts[i] = i+1;
  }

  Teuchos::RCP<Epetra_Vector> edge_weights = 
    Teuchos::rcp(new Epetra_Vector(Copy, myRowMap, ewgts));

  delete [] ewgts;

  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs = 
    Teuchos::rcp(new Isorropia::Epetra::CostDescriber);
  costs->setHypergraphEdgeWeights(edge_weights);

  // Compute balance and cuts in current partitioning

  // Perform hyperedge partitioning with Zoltan

  Teuchos::RCP<const Epetra_RowMatrix> rowmatrix(matrix_ptr);

  Teuchos::RCP<Isorropia::Partitioner> partitioner =
    Isorropia::Epetra::create_partitioner(rowmatrix, costs, params);

  Isorropia::Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsMatrix> new_matrix = rd.redistribute(*matrix_ptr);

  show_matrix("After load balancing", new_matrix.get(), Comm);

  // What is the difference between using Redistributor and doing this:
  //
  // Teuchos::RCP<Epetra_Map> pmap = 
  //       Isorropia::Epetra::create_target_map(Comm, partitioner);
  // Teuchos::RCP<Epetra_CrsMatrix> new_matrix = 
  //       Isorropia::Epetra::redistribute_rows(row_matrix, pmap);
  //

  // Recompute balance and cuts 

  MPI_Finalize();

  return(0);
}
static void show_matrix(const char *txt, Epetra_CrsMatrix *matrix, const Epetra_Comm &Comm)
{
  int localProc = Comm.MyPID();
  int numProcs = Comm.NumProc();

  int base = matrix->IndexBase();

  if (localProc == 0){
    std::cout << txt << std::endl;
    std::cout << matrix->NumGlobalNonzeros() << " non zeroes" << std::endl;
    std::cout << matrix->NumGlobalRows() << " rows" << std::endl;
    std::cout << matrix->NumGlobalCols() << " columns" << std::endl;
    std::cout << matrix->IndexBase() << " base" << std::endl;
  }
  Comm.Barrier();

  for (int p=0; p < numProcs; p++){
    if (p == localProc){
      std::cout << "Process " << p << " has rows: ";
      for (int i=0; i< matrix->NumMyRows(); i++){
        int gid = matrix->GRID(i);
        std::cout << gid << " ";
      }
      std::cout << std::endl;
    }
    Comm.Barrier();
  }
  if (localProc == numProcs-1){
    std::cout << "==============" << std::endl;
  }
}


#else // HAVE_MPI


int main()
{
  return 0;
}


#endif // HAVE_MPI

