
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
// Generate a linearproblem using triutils and redistribute with Zoltan.
// This version uses transforms.
//
// *********** Under construction! Not working yet! ********************** 

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
#include "EpetraExt_Transform.h"
#include "EpetraExt_LPTrans_From_GraphTrans.h"

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  int n=4;

  // Epetra_Map  Map;
  // Epetra_CrsMatrix A;
  // Epetra_Vector x, b, xexact;
   
  // Generate Laplacian2d gallery matrix
  Trilinos_Util::CrsMatrixGallery G("laplace_2d", Comm);
  G.Set("problem_size", n*n);
  G.Set("map_type", "linear"); // Linear map initially
  //G.Set("exact_solution", "random");

  // Get the LinearProblem. 
  Epetra_LinearProblem *Prob = G.GetLinearProblem();

  // Get the exact solution.
  Epetra_MultiVector *sol = G.GetExactSolution();

  // Get the rhs (b) and lhs (x)
  Epetra_MultiVector *b = Prob->GetRHS();
  Epetra_MultiVector *x = Prob->GetLHS();

  cout << "RHS b: " << *b << endl;
  cout << "LHS x: " << *x << endl;
  cout << "Solution: " << *sol << endl;

  // Repartition graph using Zoltan
  // Epetra_CrsGraph NewGraph = EpetraExt::Zoltan_CrsGraph(Graph);
  EpetraExt::Zoltan_CrsGraph * ZoltanTrans = new EpetraExt::Zoltan_CrsGraph();
  EpetraExt::LinearProblem_GraphTrans * ZoltanLPTrans =
    new EpetraExt::LinearProblem_GraphTrans(
         *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(ZoltanTrans)) );
 
  cout << "Creating Load Balanced Linear Problem\n";
  Epetra_LinearProblem &BalancedProb = (*ZoltanLPTrans)(*Prob);

  // Get the rhs (b) and lhs (x)
  Epetra_MultiVector *Balancedb = Prob->GetRHS();
  Epetra_MultiVector *Balancedx = Prob->GetLHS();
  cout << "Balanced b: " << *Balancedb << endl;
  cout << "Balanced x: " << *Balancedx << endl;


/****************** OLD STUFF not being used **********************

  int NumGlobalElements = readMap->NumGlobalElements();

  // Create uniform distributed map
  Epetra_Map map(NumGlobalElements, 0, Comm);

  // Create Exporter to distribute read-in matrix and vectors

  Epetra_Export exporter(*readMap, map);
  Epetra_CrsMatrix A(Copy, map, 0);
  Epetra_Vector x(map);
  Epetra_Vector b(map);
  Epetra_Vector xexact(map);

  Epetra_Time FillTimer(Comm);
  x.Export(*readx, exporter, Add);
  b.Export(*readb, exporter, Add);
  xexact.Export(*readxexact, exporter, Add);
  Comm.Barrier();
  double vectorRedistributeTime = FillTimer.ElapsedTime();
  A.Export(*readA, exporter, Add);
  Comm.Barrier();
  double matrixRedistributeTime = FillTimer.ElapsedTime() - vectorRedistributeTime;
  assert(A.TransformToLocal()==0);    
  Comm.Barrier();
  double fillCompleteTime = FillTimer.ElapsedTime() - matrixRedistributeTime;

  if( MyPID==0 ) {
    cout << "Vector redistribute  time (sec) = "
	 << vectorRedistributeTime<< endl;
    cout << "Matrix redistribute time (sec) = "
	 << matrixRedistributeTime << endl;
    cout << "Transform to Local  time (sec) = "
	 << fillCompleteTime << endl<< endl;
  }

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;

***********/

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
