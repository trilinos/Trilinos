
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
// Read a matrix in HB format
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
#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

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

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
