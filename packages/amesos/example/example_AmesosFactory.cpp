
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
// This example needs triutils to generate the linear system.
#ifdef HAVE_AMESOS_TRIUTILS
#include "Amesos.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
// following header file and namespace declaration
// are  required by this example to generate the linear system,
// not by Amesos itself.
#include "Trilinos_Util_CrsMatrixGallery.h"
using namespace Trilinos_Util;

// ==================== //
// M A I N  D R I V E R //
// ==================== //
//
// This example will:
// 1.- create a linear system, stored as
//     Epetra_LinearProblem. The matrix corresponds
//     to a 5pt Laplacian (2D on Cartesian grid).
//     The user can change the size of the problem 
//     by modifying variable NumGlobalRows.
// 2.- The linear system matrix, solution and rhs
//     are distributed among the available processors,
//     using a linear distribution. This is for 
//     simplicity only! Amesos can support any
//     Epetra_Map.
// 3.- Once the linear problem is created, we
//     create an Amesos Factory object.
// 4.- With the Factory, we create the required Amesos_BaseSolver
//     solver. Any supported (and compiled) Amesos
//     solver can be used. If the selected solver
//     is not available (that is, if Amesos has *not*
//     been configured with support for this solver),
//     the factory returns 0. Usually, Amesos_Klu
//     is always avaiable.
// 5.- At this point we can factorize the matrix,
//     and solve the linear system. Only three methods
//     should be used for an Amesos_BaseSolver object:
//     1.- NumericFactorization();
//     2.- SymbolicFactorization();
//     3.- Solve();
// 6.- We note that the header files of Amesos-supported
//     libraries are *not* required in this file. They are
//     actually needed to compile the Amesos library only.
//
// NOTE: this example can be run with one or more processors.

int main(int argc, char *argv[]) 
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumGlobalRows = 10000; // must be a square for the
                             // matrix generator.
  int NumVectors = 1;        // number of rhs's. Amesos
                             // supports single or
			     // multiple RHS.

  // initialize an Gallery object.
  // NOTE: only this example needs triutils, to 
  // build in an easy way the linear system matrix.
  // The user can easily change the matrix type;
  // consult the Trilinos tutorial on the triutils
  // chapter for more details.
  //
  // Amesos itself is INDEPENDENT from triutils.

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumGlobalRows);
  Gallery.Set("num_vectors", NumVectors);

  // get pointers to Gallery's objects. Matrix, LHS and
  // RHS are constructed by Gallery. The matrix is actually
  // stored as Epetra_CrsMatrix.
  // Problem will be used in the Amesos contruction.

  Epetra_RowMatrix* Matrix = Gallery.GetMatrix();
  Epetra_MultiVector* LHS  = Gallery.GetStartingSolution();
  Epetra_MultiVector* RHS  = Gallery.GetRHS();
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  // random RHS and zero LHS
  RHS->Random();
  LHS->PutScalar(0.0);

  // ===================================================== //
  // B E G I N N I N G   O F  T H E   AM E S O S   P A R T //
  // ===================================================== //

  // initialize Amesos solver. This is the base class for
  // Amesos. It is a pure virtual class (hence object of this
  // class cannot be allocated, only pointers or references).

  Amesos_BaseSolver * Solver;

  // initialize the Factory. Factory is a function class (a
  // class that contains methods only, no data). Factory
  // will be used to create Amesos_BaseSolver derived
  // objects.

  Amesos Factory;

  // solver can assume one of the following values:
  // - Amesos_Klu: for KLU solver
  // - Amesos_Superlu: for SuperLU
  // - Amesos_Superludist: for SuperLU_DIST 2.0 or later
  // - Amesos_MUMPS: for MUMPS 4.3.2 or later
  // 
  // Note that users can change solver simply changing
  // this parameter!

  string SolverType = "Amesos_Klu";
  Solver = Factory.Create((char*)SolverType.c_str(), *Problem);

  // Factory.Create() returns 0 if the requested solver
  // is not available
 
  if (Solver == 0) 
    return(EXIT_FAILURE);

  // Parameters for all Amesos solvers are set through
  // a call to SetParameters(List). List is a Teuchos
  // parameter list (Amesos requires Teuchos to compile).
  // In most cases, users can proceed without calling
  // SetParameters(). Please refer to the Amesos guide
  // for more details.
  //
  // Parameters in the list are set using 
  // List.set("parameter-name", ParameterValue);
  // In this example, we specify that we want more output.

  Teuchos::ParameterList List;
  List.set("PrintTiminig", true);
  List.set("PrintStatus", true);
  
  Solver->SetParameters(List);
  
  // Now we are ready to solve. Generally, users will
  // call SymbolicFactorization(), then NumericFactorization(),
  // and finally Solve(). Note that:
  // - the numerical values of the linear system matrix
  //   are *not* required before NumericFactorization();
  // - solution and rhs are *not* required before
  //   Solve().

  Solver->SymbolicFactorization();
  Solver->NumericFactorization();
  Solver->Solve();

  // =========================================== //
  // E N D   O F   T H E   A M E S O S   P A R T //
  // =========================================== //

  // At this point we can check that the residual is
  // small as it should be. NOTE: this check can be
  // performed inside Amesos as well. Use
  // List.set("ComputeTrueResidual",true) before
  // calling SetParameters(List).
  //
  double* residual;
  residual = new double[NumVectors];

  Gallery.ComputeResidual(residual);

  if( Comm.MyPID() == 0 ) {
    for (int i = 0 ; i < NumVectors ; ++i)
      cout << "After AMESOS solution, for vector " << i << ", ||b-Ax||_2 = " << residual[i] << endl;
  }

  // delete Solver. MPI calls can occur.
  delete Solver;
  delete [] residual;
    
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (residual[0] < 1e-5)
    return(EXIT_SUCCESS);
  else
    return(EXIT_FAILURE);

} // end of main()

#else

// Triutils is not available. Sorry, we have to give up.

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure AMESOS with --enable-triutils");
  puts("to run this example");
  
  return 0;
}

#endif

