
// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// Trilinos Tutorial
// -----------------
// Generate a matrix using triutils and redistribute with Zoltan.

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "EpetraExt_Zoltan_CrsGraph.h"
// #include "EpetraExt_Transform.h"

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  int n = 4; /* n by n grid */

  // Generate Laplacian2d gallery matrix
  Trilinos_Util::CrsMatrixGallery G("laplace_2d", Comm);
  G.Set("problem_size", n*n);
  G.Set("map_type", "linear"); // Linear map initially

  // Copy the CrsMatrix to A. 
  Epetra_CrsMatrix A (*(G.GetMatrix()));
  // Epetra_CrsGraph Graph (A.Graph());

  cout << "Matrix A: " << A << endl;
  // cout << "Graph of A: " << A.Graph() << endl;

  // Repartition graph using Zoltan
  cout << "Calling Zoltan via EpetraExt to repartition the graph." << endl;
  // First create a transform 
  EpetraExt::Zoltan_CrsGraph ZoltanTrans;
  // Then apply the transform to the graph.
  Epetra_CrsGraph & BalGraph = ZoltanTrans(const_cast<Epetra_CrsGraph&>(A.Graph()));
  // Epetra_CrsGraph & BalGraph = ZoltanTrans(Graph);

  // cout << "Graph of A (not changed after balancing) : " << A.Graph() << endl;
  // cout << "Balanced Graph of A: " << BalGraph << endl;

  // Create Exporter to rebalance the matrix from the row maps
  Epetra_Export exporter(A.Graph().RowMap(), BalGraph.RowMap());

  cout << "Old rowmap for A: " << A.Graph().RowMap() << endl;
  cout << "New rowmap for A: " << BalGraph.RowMap() << endl;

  cout << "Proc " << MyPID << ": NumSend = " << exporter.NumSend()
       << ", NumRecv = " << exporter.NumRecv() << endl;

  Epetra_Time Timer(Comm);
  Comm.Barrier();
  double startTime = Timer.ElapsedTime();
  // Export A to new distribution. 
  Epetra_CrsMatrix BalA(Copy, BalGraph);
  BalA.Export(A, exporter, Insert);
  Comm.Barrier();
  double matrixRedistributeTime = Timer.ElapsedTime() - startTime;

  cout << "Rebalanced matrix: " << BalA << endl ;
  if( MyPID==0 ) {
    cout << "Matrix redistribute time (sec) = "
	 << matrixRedistributeTime << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
