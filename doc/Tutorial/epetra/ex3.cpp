
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
// define a vector of global dimension 5, split among 2 processes p0 and p1
// int the following way:
//  0 -- 1 -- 2 -- 3 -- 4
//  <-p0->    <----p1--->
// This code must be run with two processes
//
// (output reported at the end of the file)

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumGlobalElements=5;  // global dimension of the problem
  int MyElements;           // local dimension of the problem
  int *MyGlobalElements;    // local-to-global map
  int MyPID = Comm.MyPID(); // ID of this process

  if( Comm.NumProc() != 2 ) {
    cerr << "This code must be run with 2 processes\n";
    exit( EXIT_FAILURE );
  }
  
  // hardwired number of local elements and local-to-global map
  switch( MyPID ) {
  case 0:
    MyElements = 2;
    MyGlobalElements = new int[MyElements];
    MyGlobalElements[0] = 0;
    MyGlobalElements[1] = 1;
    break;
  case 1:
    MyElements = 3;
    MyGlobalElements = new int[MyElements];
    MyGlobalElements[0] = 2;
    MyGlobalElements[1] = 3;
    MyGlobalElements[2] = 4;
    break;
  }

  Epetra_Map Map(NumGlobalElements,MyElements,MyGlobalElements,0,Comm);
  
  cout << Map;

  delete MyGlobalElements;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0 ;

} /* main */

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:epetra]> mpirun -np 2 ./ex3.exe
Epetra::Map

Number of Global Elements  = 5
Number of Global Points = 5
Maximum of all GIDs        = 4
Minimum of all GIDs        = 0
Index Base                 = 0
Constant Element Size      = 1

Number of Local Elements   = 2
Number of Local Points  = 2
Maximum of my GIDs         = 1
Minimum of my GIDs         = 0

         MyPID           Local Index        Global Index
             0                 0                 0
             0                 1                 1
Epetra::Map

Number of Local Elements   = 3
Number of Local Points  = 3
Maximum of my GIDs         = 4
Minimum of my GIDs         = 2

         MyPID           Local Index        Global Index
             1                 0                 2
             1                 1                 3
             1                 2                 4

*/
