
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
// Use of Epetra_FECrsMatrix
// This code should be run with at least two processes
//
// (output reported at the end of the file)
//
// Marzio Sala, SNL, 9214, 20-Nov-2003

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"

// =========== //
// main driver //
// ----------- //

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // create a linear map of size 5 (could be any number > 1)
  int NumGlobalElements = 5;
  Epetra_Map Map(NumGlobalElements,0,Comm);

  // create a diagonal FE crs matrix (one nonzero per row)
  Epetra_FECrsMatrix A(Copy,Map,2);
  
  // set the entries
  if( Comm.MyPID() == 0 ) {
    for( int i=0 ; i<NumGlobalElements ; ++i ) {
      int indices[2];
      indices[0] = i; indices[1] = i;
      double value = 1.0*i;
      A.SumIntoGlobalValues(1,indices,&value);
    }
  }
  
  A.GlobalAssemble();

  cout << A;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}

