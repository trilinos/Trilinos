
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

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

// ==================== //
// M A I N  D R I V E R //
// ==================== //
//
// This example will:
// 1.- Create an tridiagonal matrix;
// 2.- Call SymbolicFactorization();
// 3.- Change the numerical values of the matrix;
// 4.- Call NumericFactorization();
// 5.- Set the entries of the RHS;
// 6.- Call Solve().
//
// This example is intended to show the required data
// for each phase. Phase at point (2) requires matrix
// structure only. Phase at point (4) requires matrix
// structure (supposed unchanged from point (2)) and
// matrix data. Phase (6) requires RHS and solution vector.
//
// NOTE: this example can be run with one or more processes.
//
// NOTE2: More details about Amesos can be found in
// example example_AmesosFactory.cpp and the
// Amesos users' guide.

int main(int argc, char *argv[]) 
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumGlobalElements = 100; // global dimension of the problem.

  // ======================================================= //
  // B E G I N N I N G   O F   M A T R I X   C R E A T I O N //
  // ======================================================= //
 
  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Epetra_Map Map(NumGlobalElements, 0, Comm);

  // Create an empty EpetraCrsMatrix 
  Epetra_CrsMatrix A(Copy, Map, 0);

  // Create the structure of the matrix (tridiagonal)

  if (Comm.MyPID() == 0)
    cout << "Defining nonzero structure of the matrix..." << endl;
  int NumMyElements = Map.NumMyElements();

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1

  double *Values = new double[3];
  Values[0] = 0.0;
  Values[1] = 0.0;
  Values[2] = 0.0;
  int *Indices = new int[3];
  int NumEntries;
  /// global ID's of local ID's
  int* MyGlobalElements = Map.MyGlobalElements();

  // At this point we simply set the nonzero structure of A.
  // Actual values will be inserted later.
  for (int i=0; i<NumMyElements; i++) {

    if (MyGlobalElements[i] == 0) {
      Indices[0] = 0; Indices[1] = 1;
      Values[0]  = 2.0; Values[1]  = -1.0;
      NumEntries = 2;
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1) {
      Indices[0] = NumGlobalElements-1;
      Indices[1] = NumGlobalElements-2;
      Values[0]  = 2.0; Values[1]  = -1.0;
      NumEntries = 2;
    }
    else {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i]+1;
      Values[0] = -1.0; Values[1] = 2.0; Values[2] = -1.0;
      NumEntries = 3;
    }

    AMESOS_CHK_ERR(A.InsertGlobalValues(MyGlobalElements[i], 
					NumEntries, Values, Indices));
  }

  assert(A.TransformToLocal()==0);

  // Finish up. Now the matrix STRUCTURE is set. We cannot add
  // nonzero elements, but we can still change the
  // numerical values of all inserted elements (as we will
  // do later).

  // =========================================== //
  // E N D   O F   M A T R I X   C R E A T I O N //
  // =========================================== //

  // ===================================================== //
  // B E G I N N I N G   O F  T H E   AM E S O S   P A R T //
  // ===================================================== //

  Epetra_LinearProblem Problem;

  Problem.SetOperator(&A);

  // initialize Amesos solver. Here we solve with Amesos_Klu.

  Amesos_BaseSolver * Solver;
  Amesos Factory;
  Solver = Factory.Create("Amesos_Klu", Problem);

  // Factory.Create() returns 0 if the requested solver
  // is not available
 
  if (Solver == 0) 
    return(EXIT_FAILURE);

  // At this point we can perform the numeric factorization.
  // Note that the matrix contains 0's only.

  if (Comm.MyPID() == 0)
    cout << "Inserting values in the matrix..." << endl;

  Solver->SymbolicFactorization();

  // Now, we repopulate the matrix with entries corresponding
  // to a 1D Laplacian. LHS and RHS are still untouched.

  for (int i=0; i<NumMyElements; i++) {

    if (MyGlobalElements[i] == 0) {
      Indices[0] = 0; Indices[1] = 1;
      Values[0]  = 2.0; Values[1]  = -1.0;
      NumEntries = 2;
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1) {
      Indices[0] = NumGlobalElements-1;
      Indices[1] = NumGlobalElements-2;
      Values[0]  = 2.0; Values[1]  = -1.0;
      NumEntries = 2;
    }
    else {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i]+1;
      Values[0] = -1.0; Values[1] = 2.0; Values[2] = -1.0;
      NumEntries = 3;
    }

    AMESOS_CHK_ERR(A.ReplaceGlobalValues(MyGlobalElements[i], 
					NumEntries, Values, Indices));
  }

  // ... and we can compute the numeric factorization.
  Solver->NumericFactorization();

  // Finally, we set up the RHS vector (Random().
  Epetra_Vector b(Map);
  b.Random();
  Epetra_Vector x(Map);
  x.PutScalar(0.0);
  
  Problem.SetLHS(&x);
  Problem.SetRHS(&b);
  
  Solver->Solve();

  // =========================================== //
  // E N D   O F   T H E   A M E S O S   P A R T //
  // =========================================== //

  // ================== //
  // compute ||Ax - b|| //
  // ================== //

  double residual;

  Epetra_Vector Ax(b.Map());
  A.Multiply(false, x, Ax);
  Ax.Update(1.0, b, -1.0);
  Ax.Norm2(&residual);

  if( Comm.MyPID() == 0 ) {
    cout << "After AMESOS solution, ||b-Ax||_2 = " << residual << endl;
  }

  // delete Solver. MPI calls can occur.
  delete Solver;
  delete [] Values;
  delete [] Indices;
    
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (residual < 1e-5)
    return(EXIT_SUCCESS);
  else
    return(EXIT_FAILURE);

} // end of main()

