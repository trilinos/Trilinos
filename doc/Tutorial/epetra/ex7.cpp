
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
// Example of Epetra_MultiVector; use of ExtractView on MultiVector.
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
#include "Epetra_MultiVector.h"

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

  // Total number of elements in the vector
  int NumElements = 10;

  // Construct a Map with NumElements and index base of 0
  Epetra_Map Map(NumElements, 0, Comm);

  // Create x as a 2-component multi-vector
  Epetra_MultiVector x(Map,2);
 
  // get the local size of the vector
  int MyLength = x.MyLength();

  /* First way to define the vector:   */
  /* use the [] operator on the object */

  for( int c=0 ; c<x.NumVectors() ; ++c ) 
    for( int i=0 ; i<MyLength ; ++i ) x[c][i] = 1.0*i+1000*c;

  // need a double pointer because this works with multi-vectors
  double ** pointer;
  
  x.ExtractView( &pointer );

  for( int c=0 ; c<x.NumVectors() ;++c ) 
    for( int i=0 ; i<MyLength ; ++i )
      cout << "on proc " << Comm.MyPID() << ", x["
	   << i << "] = " << pointer[c][i] << endl;

  // now modify the values
  for( int c=0 ; c<x.NumVectors() ;++c ) 
    for( int i=0 ; i<MyLength ; ++i )
      pointer[c][i] *= 10;

  cout << x;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
} /* main */

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[sala:epetra]> mpirun -np 2 ./ex7.exe
on proc 0, x[0] = 0
on proc 0, x[1] = 1
on proc 0, x[2] = 2
on proc 0, x[3] = 3
on proc 0, x[4] = 4
on proc 0, x[0] = 1000
on proc 0, x[1] = 1001
on proc 0, x[2] = 1002
on proc 0, x[3] = 1003
on proc 0, x[4] = 1004
Epetra::MultiVector
     MyPID           GID               Value               Value
         0             0                       0               10000
         0             1                      10               10010
         0             2                      20               10020
         0             3                      30               10030
         0             4                      40               10040
on proc 1, x[0] = 0
on proc 1, x[1] = 1
on proc 1, x[2] = 2
on proc 1, x[3] = 3
on proc 1, x[4] = 4
on proc 1, x[0] = 1000
on proc 1, x[1] = 1001
on proc 1, x[2] = 1002
on proc 1, x[3] = 1003
on proc 1, x[4] = 1004
Epetra::MultiVector
         1             5                       0               10000
         1             6                      10               10010
         1             7                      20               10020
         1             8                      30               10030
         1             9                      40               10040

*/


