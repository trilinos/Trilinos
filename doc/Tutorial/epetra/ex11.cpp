
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
// Basic definition of communicator.
// This code should be run with one process
//
// (output reported at the end of the file)


#include <iostream>
#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Total number of elements in vectors, can be any positive number
  int NumRows = 5;

  Epetra_SerialDenseVector x, b;
  x.Size( NumRows );
  b.Size( NumRows );

  // set the elements of the vector
  for( int i=0 ; i<NumRows ; ++i ) b[i] = 1.0, x[i]=0.0;
  
  Epetra_SerialDenseMatrix A, A2;
  A.Shape( NumRows, NumRows );
  A2.Shape( NumRows, NumRows ); // A2 is a copy of A

  // Hilbert matrix (ill-conditioned)
  for( int i=0 ; i<NumRows ; ++i )
    for( int j=0 ; j<NumRows ; ++j ) 
      A(i,j) = 1.0/(i+j+2);

  cout<< A;

  // set up the solver
  Epetra_SerialDenseSolver Problem;
  Problem.SetMatrix( A );
  Problem.SetVectors( x, b );

  A2 = A; 
  // we make a copy of A because Problem.Solve() will
  // overwrite A with its LU decomposition. Try with
  // cout << A after the following invocation

  b.Multiply('N','N',1.0, A2, x, 0.0);

  cout << "A * x = \n" << b;

  double rcond;
  Problem.ReciprocalConditionEstimate(rcond);
  cout << "The (estimated) condition number of A is " << 1/rcond << endl;

  Problem.SetMatrix( A2 );
  Problem.Invert();
  cout << "The inverse of A is\n";
  cout << A2;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

} /* main */

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:epetra]> mpirun -np 2 ./ex1.exe
Epetra::SerialDenseMatrix

Data access mode: Copy
A_Copied: yes
Rows(M): 5
Columns(N): 5
LDA: 5
0.5 0.333333 0.25 0.2 0.166667
0.333333 0.25 0.2 0.166667 0.142857
0.25 0.2 0.166667 0.142857 0.125
0.2 0.166667 0.142857 0.125 0.111111
0.166667 0.142857 0.125 0.111111 0.1
A * x =
Epetra::SerialDenseVector
Data access mode: Copy
A_Copied: yes
Length(M): 5
0 0 0 0 0
The (estimated) condition number of A is 2.81723e+06
The inverse of A is
Epetra::SerialDenseMatrix

Data access mode: Copy
A_Copied: yes
Rows(M): 5
Columns(N): 5
LDA: 5
450 -4200 12600 -15120 6300
-4200 44100 -141120 176400 -75600
12600 -141120 470400 -604800 264600
-15120 176400 -604800 793800 -352800
6300 -75600 264600 -352800 158760

*/
