
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
// Definition of Epetra_Map
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
#include "Epetra_Map.h"

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // create a map given the global number of points (any positive number)  
  int NumGlobalPoint = 4;

  Epetra_Map Map1(NumGlobalPoint,0,Comm);

  // Epetra_Map overloads the << operator
  cout << Map1;

  // now create a map given the local number of points

  int NumMyPoints = Comm.MyPID();
  
  Epetra_Map Map2(-1,NumMyPoints,0,Comm);

  cout << Map2;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

} /* main */

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[sala:epetra]> mpirun -np 2 ./ex2.exe

Epetra::Map

Number of Global Elements  = 4
Number of Global Points = 4
Maximum of all GIDs        = 3
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

Number of Local Elements   = 2
Number of Local Points  = 2
Maximum of my GIDs         = 3
Minimum of my GIDs         = 2

         MyPID           Local Index        Global Index
             1                 0                 2
             1                 1                 3
Epetra::Map

Number of Global Elements  = 3
Number of Global Points = 3
Maximum of all GIDs        = 2
Minimum of all GIDs        = 0
Index Base                 = 0
Constant Element Size      = 1

Number of Local Elements   = 1
Number of Local Points  = 1
Maximum of my GIDs         = 0
Minimum of my GIDs         = 0

         MyPID           Local Index        Global Index
             0                 0                 0
Epetra::Map

Number of Local Elements   = 2
Number of Local Points  = 2
Maximum of my GIDs         = 2
Minimum of my GIDs         = 1

         MyPID           Local Index        Global Index
             1                 0                 1
             1                 1                 2
*/
