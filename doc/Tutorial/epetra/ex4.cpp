
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
// Define serial dense vectors
// This code should be run with at least two processes
//
// (output reported at the end of the file)

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_SerialDenseVector.h"

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Total number of elements in vectors
  int NumRows = 5;

  Epetra_SerialDenseVector x;

  // shape this object
  x.Size( NumRows );

  // set the elements of the vector
  for( int i=0 ; i<NumRows ; ++i )  x[i] = 1.0*i;
  // NOTE: x[i] performs the same operations of x(i); 
  // however, the former checks for the bounds, while the latter
  // do such only if Epetra is compiled with -DEPETRA_ARRAY_BOUNDS_CHECK


  // now try to write out of memory
  //int Length = x.Length();
  
  //x[Length] = 12;
  
  cout << x;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( 0 );

} /* main */

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:epetra]> mpirun -np  1 ./ex4.exe
Epetra::SerialDenseVector
Data access mode: Copy
A_Copied: yes
Length(M): 5
0 1 2 3 4

*/
