
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
// Define dense serial matrices
// This code should be run with one process
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
#include "Epetra_SerialDenseMatrix.h"

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

  // Total number of elements in matrix A
  int NumRowsA = 2, NumColsA = 2;

  Epetra_SerialDenseMatrix A, B;
  A.Shape( NumRowsA, NumColsA );

  // set the element of A using the () operator.
  // Note that i is the row-index, and j the column-index
  for( int i=0 ; i<NumRowsA ; ++i ) 
    for( int j=0 ; j<NumColsA ; ++j ) 
      A(i,j) = i+100*j;

  cout << A;

  cout << "Inf norm of A = " << A.OneNorm() << endl;
  cout << "One norm of A = " << A.InfNorm() << endl;

  // now define an other matrix, B, for matrix multiplication
  int NumRowsB = 2, NumColsB=1;
  B.Shape(NumRowsB, NumColsB);

  // enter the values of B
  for( int i=0 ; i<NumRowsB ; ++i ) 
    for( int j=0 ; j<NumColsB ; ++j ) 
      B(i,j) = 11.0+i+100*j;

  cout << B;

  // define the matrix which will hold A * B
  Epetra_SerialDenseMatrix AtimesB;

  // same number of rows than A, same columns than B
  AtimesB.Shape(NumRowsA,NumColsB);  

  // A * B
  AtimesB.Multiply('N','N',1.0, A, B, 0.0);
  cout << AtimesB;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
} /* main */

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:epetra]> mpirun -np 1 ./ex10.exe
Epetra::SerialDenseMatrix

Data access mode: Copy
A_Copied: yes
Rows(M): 2
Columns(N): 2
LDA: 2
0 100
1 101
Inf norm of A = 201
One norm of A = 102
Epetra::SerialDenseMatrix

Data access mode: Copy
A_Copied: yes
Rows(M): 2
Columns(N): 1
LDA: 2
11
12
Epetra::SerialDenseMatrix

Data access mode: Copy
A_Copied: yes
Rows(M): 2
Columns(N): 1
LDA: 2
1200
1223

*/

     
