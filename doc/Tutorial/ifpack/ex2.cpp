
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
// Using IFPACK factorizations as Aztec's preconditioners

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
#include "AztecOO.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"

// function for fancy output

string toString(const int& x) {
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

string toString(const double& x) {
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

// main driver

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  bool verbose = false; 
  if (MyPID==0) verbose = true;

  // B E G I N   O F   M A T R I X   C O N S T R U C T I O N
  
  // matrix downloaded from MatrixMarket
  char FileName[] = "../HBMatrices/bcsstk14.rsa";

  Epetra_Map * readMap; // Pointers because of Trilinos_Util_ReadHb2Epetra
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
   
  // Call routine to read in HB problem
  Trilinos_Util_ReadHb2Epetra(FileName, Comm, readMap, readA, readx, 
			      readb, readxexact);

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
  A.Export(*readA, exporter, Add);
  x.Export(*readx, exporter, Add);
  b.Export(*readb, exporter, Add);
  xexact.Export(*readxexact, exporter, Add);

  A.TransformToLocal();
  
  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;

  // E N D   O F   M A T R I X   C O N S T R U C T I O N  

  // ============================= //
  // Construct RILU preconditioner //
  // ---------------=------------- //

  //  modify those parameters 
  int    LevelFill = 0;
  int    Overlap = 2;
  double Athresh = 0.0;
  double Rthresh = 1.0;

  Ifpack_IlukGraph * Graph = 0;
  Ifpack_CrsRiluk * RILU = 0;

  Graph = new Ifpack_IlukGraph(A.Graph(), LevelFill, Overlap);
  assert(Graph->ConstructFilledGraph()==0);

  RILU = new Ifpack_CrsRiluk(*Graph);
  int initerr = RILU->InitValues(A);
  if (initerr!=0) cout << Comm << "*ERR* InitValues = " << initerr;

  assert(RILU->Factor()==0);
  
  // Define label for printing out during the solve phase
  string label = "Ifpack_CrsRiluk Preconditioner: LevelFill = " + toString(LevelFill) + 
                                                 " Overlap = " + toString(Overlap) + 
                                                 " Athresh = " + toString(Athresh) + 
                                                 " Rthresh = " + toString(Rthresh); 
  RILU->SetLabel(label.c_str());

  // Here we create an AztecOO object
  AztecOO solver;
  solver.SetUserMatrix(&A);
  solver.SetLHS(&x);
  solver.SetRHS(&b);

  // Here we set the IFPACK preconditioner and specify few parameters
  
  solver.SetPrecOperator(RILU);

  int Niters = 1200;
  solver.SetAztecOption(AZ_kspace, Niters); 
  solver.Iterate(Niters, 5.0e-10);

  if (RILU!=0) delete RILU;
  if (Graph!=0) delete Graph;
				       
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
