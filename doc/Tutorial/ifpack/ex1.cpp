
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
//
// NOTE: this example implemenets minor modifications to one of the
// examples included in the AztecOO package. Please give a look
// to file ${TRILINOS_HOME}/packages/aztecoo/examples/IfpackAztecOO/cxx_main.cpp
// for more details.
//
// (output reported at the end of the file)
//
// Marzio Sala, SNL, 9214, 19-Nov-2003

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "AztecOO.h"
#include "Trilinos_Util.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_CrsIct.h"

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
  
  // ============================ //
  // Construct ILU preconditioner //
  // ---------------------------- //

  //  modify those parameters 
  int    LevelFill = 1;
  double DropTol = 0.0;
  double Condest;
  
  Ifpack_CrsIct * ICT = NULL;
  ICT = new Ifpack_CrsIct(A,DropTol,LevelFill);
  // Init values from A
  ICT->InitValues(A);
  // compute the factors
  ICT->Factor();
  // and now estimate the condition number
  ICT->Condest(false,Condest);
  
  cout << Condest << endl;
    
  if( Comm.MyPID() == 0 ) {
    cout << "Condition number estimate (level-of-fill = "
	 << LevelFill <<  ") = " << Condest << endl;
  }

  // Define label for printing out during the solve phase
  string label = "Ifpack_CrsIct Preconditioner: LevelFill = " + toString(LevelFill) + 
                                                 " Overlap = 0"; 
  ICT->SetLabel(label.c_str());
  
  // Here we create an AztecOO object
  AztecOO solver;
  solver.SetUserMatrix(&A);
  solver.SetLHS(&x);
  solver.SetRHS(&b);
  solver.SetAztecOption(AZ_solver,AZ_cg);
  
  // Here we set the IFPACK preconditioner and specify few parameters
  
  solver.SetPrecOperator(ICT);

  int Niters = 1200;
  solver.SetAztecOption(AZ_kspace, Niters);
  solver.SetAztecOption(AZ_output, 20); 
  solver.Iterate(Niters, 5.0e-5);

  if (ICT!=0) delete ICT;
				       
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:ifpack]> mpirun -np 2 ./ex1.exe

*/
