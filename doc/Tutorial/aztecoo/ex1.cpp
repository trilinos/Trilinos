
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
// Solve a 2D Laplacian problem
//
// This example builds the matrix and solves it with AztecOO. 

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "AztecOO.h"

// external function
void  get_neighbours( const int i, const int nx, const int ny,
		      int & left, int & right, 
		      int & lower, int & upper);

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // number of nodes along the x- and y-axis
  int nx = 5;
  int ny = 6;
  int NumGlobalElements = nx * ny;

  // create a linear map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  
  // local number of rows
  int NumMyElements = Map.NumMyElements();
  // get update list
  int * MyGlobalElements = new int [NumMyElements];
  Map.MyGlobalElements( MyGlobalElements );

  // Create a Epetra_Matrix with 5 nonzero per rows
  
  Epetra_CrsMatrix A(Copy,Map,5);

  // Add  rows one-at-a-time
  // Need some vectors to help

  double Values[4];
  int Indices[4];
  int NumEntries;
  int left, right, lower, upper;
  double diag = 4.0;
  
  for( int i=0 ; i<NumMyElements; ++i ) {
    int NumEntries=0;
    get_neighbours(  MyGlobalElements[i], nx, ny, 
		     left, right, lower, upper);
    if( left != -1 ) {
	Indices[NumEntries] = left;
	Values[NumEntries] = -1.0;
	++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = -1.0;
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, 
				Values, Indices)==0);
    // Put in the diagonal entry
    assert(A.InsertGlobalValues(MyGlobalElements[i], 1, 
				&diag, MyGlobalElements+i)==0);
  }
  
  // Finish up
  assert(A.TransformToLocal()==0);

  // create x and b vectors
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);
  b.PutScalar(1.0);

  // ==================== AZTECOO INTERFACE ======================
  
  // create linear problem
  Epetra_LinearProblem Problem(&A,&x,&b);
  // create AztecOO instance
  AztecOO Solver(Problem);

  Solver.SetAztecOption( AZ_precond, AZ_Jacobi );
  
  Solver.Iterate(1000,1E-9);

  // ==================== END OF AZTECOO INTERFACE ================
  
  if( Comm.MyPID() == 0 ) {
    cout << "Solver performed " << Solver.NumIters() 
	 << "iterations.\n";
    cout << "Norm of the true residual = " << Solver.TrueResidual() << endl;
  }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void  get_neighbours( const int i, const int nx, const int ny,
		      int & left, int & right, 
		      int & lower, int & upper) 
{

  int ix, iy;
  ix = i%nx;
  iy = (i - ix)/nx;

  if( ix == 0 ) 
    left = -1;
  else 
    left = i-1;
  if( ix == nx-1 ) 
    right = -1;
  else
    right = i+1;
  if( iy == 0 ) 
    lower = -1;
  else
    lower = i-nx;
  if( iy == ny-1 ) 
    upper = -1;
  else
    upper = i+nx;

  return;

}
