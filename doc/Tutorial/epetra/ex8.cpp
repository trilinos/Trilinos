
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
// Epetra_Vectors in View mode; use of ResetView.
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
#include "Epetra_Vector.h"

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

  // Total number of elements in the vector
  int NumLocalElements = 10;

  // Construct a Map with NumElements and index base of 0
  Epetra_Map Map(-1, NumLocalElements, 0, Comm);

  // here it is defined a double vector of size NumLocalElements
  // (that is, a size compatible with Map), and it is filled with
  // some values

  double values[NumLocalElements];
  for( int i=0 ;i<NumLocalElements ; i++ )
    values[i] = 1.0*i;
  
  // Create x as an Epetra_vector with the View mode, using `values'
  // as data
  
  Epetra_Vector x(View, Map, values);

  // now we can change x by modifying values...
  values[0] = 123;
  
  // this affects the object x
  cout << x;

  // now we can reset the view, and let it point to an other douoble
  // vector (having the same size)

  double values2[NumLocalElements];
  for( int i=0 ;i<NumLocalElements ; i++ )
    values2[i] = -1.0*i;

  x.ResetView( values2 );

  cout << x;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
} /* main */

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[sala:epetra]> mpirun -np 1 ./ex8.exe
Epetra::Vector
     MyPID           GID               Value
         0             0                     123
         0             1                       1
         0             2                       2
         0             3                       3
         0             4                       4
         0             5                       5
         0             6                       6
         0             7                       7
         0             8                       8
         0             9                       9
Epetra::Vector
     MyPID           GID               Value
         0             0                      -0
         0             1                      -1
         0             2                      -2
         0             3                      -3
         0             4                      -4
         0             5                      -5
         0             6                      -6
         0             7                      -7
         0             8                      -8
         0             9                      -9
*/


