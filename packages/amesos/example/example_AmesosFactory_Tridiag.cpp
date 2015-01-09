
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Amesos_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"

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
// for each phase. Phase (2) requires the matrix structure only. 
// Phase (4) requires the matrix structure (supposed unchanged 
// from phase (2)) and the matrix data. Phase (6) requires the 
// RHS and solution vector.
//
// This example can be run with any number of processors.
//
// Author: Marzio Sala, SNL 9214
// Last modified: Apr-05.

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
 
  // Construct a Map that puts approximatively the same number of 
  // equations on each processor. `0' is the index base (that is,
  // numbering starts from 0.
  Epetra_Map Map(NumGlobalElements, 0, Comm);

  // Create an empty EpetraCrsMatrix 
  Epetra_CrsMatrix A(Copy, Map, 0);

  // Create the structure of the matrix (tridiagonal)
  int NumMyElements = Map.NumMyElements();

  // Add  rows one-at-a-time
  // Need some vectors to help

  double Values[3];
  // Right now, we put zeros only in the matrix.
  Values[0] = 0.0;
  Values[1] = 0.0;
  Values[2] = 0.0;
  int Indices[3];
  int NumEntries;
  /// global ID's of local ID's
  int* MyGlobalElements = Map.MyGlobalElements();

  // At this point we simply set the nonzero structure of A.
  // Actual values will be inserted later (now all zeros)
  for (int i = 0; i < NumMyElements; i++) 
  {
    if (MyGlobalElements[i] == 0) 
    {
      Indices[0] = 0; 
      Indices[1] = 1;
      NumEntries = 2;
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1) 
    {
      Indices[0] = NumGlobalElements-1;
      Indices[1] = NumGlobalElements-2;
      NumEntries = 2;
    }
    else 
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i]+1;
      NumEntries = 3;
    }

    AMESOS_CHK_ERR(A.InsertGlobalValues(MyGlobalElements[i], 
					NumEntries, Values, Indices));
  }

  // Finish up.
  A.FillComplete();

  // =========================================== //
  // E N D   O F   M A T R I X   C R E A T I O N //
  // =========================================== //

  // Now the matrix STRUCTURE is set. We cannot add
  // new nonzero elements, but we can still change the
  // numerical values of all inserted elements (as we will
  // do later).

  // ===================================================== //
  // B E G I N N I N G   O F  T H E   AM E S O S   P A R T //
  // ===================================================== //

  // For comments on the commands in this section, please
  // see file example_AmesosFactory.cpp.
  
  Epetra_LinearProblem Problem;

  Problem.SetOperator(&A);

  // Initializes Amesos solver. Here we solve with Amesos_Klu.

  Amesos_BaseSolver * Solver;
  Amesos Factory;
  Solver = Factory.Create("Amesos_Klu", Problem);

  // Factory.Create() returns 0 if the requested solver
  // is not available
 
  if (Solver == 0) {
    std::cerr << "Selected solver is not available" << std::endl;
    // return ok not to break the test harness
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  // At this point we can perform the numeric factorization.
  // Note that the matrix contains 0's only.

  Solver->SymbolicFactorization();

  // Now, we repopulate the matrix with entries corresponding
  // to a 1D Laplacian. LHS and RHS are still untouched.

  for (int i = 0; i < NumMyElements; i++) 
  {
    if (MyGlobalElements[i] == 0) 
    {
      Indices[0] = 0;   
      Indices[1] = 1;
      Values[0]  = 2.0; 
      Values[1]  = -1.0;
      NumEntries = 2;
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1) 
    {
      Indices[0] = NumGlobalElements - 1;
      Indices[1] = NumGlobalElements - 2;
      Values[0]  = 2.0; 
      Values[1]  = -1.0;
      NumEntries = 2;
    }
    else 
    {
      Indices[0] = MyGlobalElements[i] - 1;
      Indices[1] = MyGlobalElements[i];
      Indices[2] = MyGlobalElements[i] + 1;
      Values[0] = -1.0; 
      Values[1] = 2.0; 
      Values[2] = -1.0;
      NumEntries = 3;
    }

    AMESOS_CHK_ERR(A.ReplaceGlobalValues(MyGlobalElements[i], 
                                         NumEntries, Values, Indices));
  }

  // ... and we can compute the numeric factorization.
  Solver->NumericFactorization();

  // Finally, we set up the LHS and the RHS vector (Random().
  Epetra_Vector b(Map);
  b.Random();
  Epetra_Vector x(Map);
  x.PutScalar(0.0);
  
  Problem.SetLHS(&x);
  Problem.SetRHS(&b);
  
  Solver->Solve();

  // Print out the timing information and get it from the solver
  Solver->PrintTiming();
  
  Teuchos::ParameterList TimingsList;
  Solver->GetTiming( TimingsList );
  
  // You can find out how much time was spent in ...
  double sfact_time, nfact_time, solve_time;
  double mtx_conv_time, mtx_redist_time, vec_redist_time;

  // 1) The symbolic factorization
  //    (parameter doesn't always exist)
  sfact_time = TimingsList.get( "Total symbolic factorization time", 0.0 );

  // 2) The numeric factorization
  //    (always exists if NumericFactorization() is called)
  nfact_time = Teuchos::getParameter<double>( TimingsList, "Total numeric factorization time" );

  // 3) Solving the linear system
  //    (always exists if Solve() is called)
  solve_time = Teuchos::getParameter<double>( TimingsList, "Total solve time" );

  // 4) Converting the matrix to the accepted format for the solver
  //    (always exists if SymbolicFactorization() is called)
  mtx_conv_time = Teuchos::getParameter<double>( TimingsList, "Total solve time" );

  // 5) Redistributing the matrix for each solve to the accepted format for the solver
  mtx_redist_time = TimingsList.get( "Total matrix redistribution time", 0.0 );

  // 6) Redistributing the vector for each solve to the accepted format for the solver
  vec_redist_time = TimingsList.get( "Total vector redistribution time", 0.0 );

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

  if (!Comm.MyPID())
    std::cout << "After AMESOS solution, ||b-Ax||_2 = " << residual << std::endl;

  // delete Solver. Do this before MPI_Finalize()
  // as MPI calls can occur in the destructor.
  delete Solver;
    
  if (residual > 1e-5)
    return(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);

} // end of main()
