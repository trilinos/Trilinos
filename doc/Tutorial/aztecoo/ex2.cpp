
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
//  Solve a linear system with Aztec, using Aztec as a preconditioner
// (recursive way)
//
// NOTE: this example implemenets minor modifications to one of the
// examples included in the AztecOO package. Please give a look  to file
// ${TRILINOS_HOME}/packages/aztecoo/examples/AztecOO_RecursiveCall/cxx_main.cpp
// for more details.
//
// (output reported at the end of the file)
//
// Marzio Sala, SNL, 9214, 20-Nov-2003

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  /* this example creates a tridiagonal matrix of type
   *
   *     |  2  -1            |
   *     | -1   2   -1       | 
   * A = |      ...  ... ... |
   *     |            -1  2  |
   */
  
  // set global dimension to 5, could be any number
  int NumGlobalElements = 5;
  
  // create a linear map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  // local number of rows
  int NumMyElements = Map.NumMyElements();
  // get update list
  int * MyGlobalElements = new int [NumMyElements];
  Map.MyGlobalElements( MyGlobalElements );

  // Create a Epetra_Matrix
  // create a CSR matrix

  Epetra_CrsMatrix A(Copy,Map,3);

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1

  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  int *Indices = new int[2];
  double two = 2.0;
  int NumEntries;

  for( int i=0 ; i<NumMyElements; ++i ) {
    if (MyGlobalElements[i]==0) {
	Indices[0] = 1;
	NumEntries = 1;
    } else if (MyGlobalElements[i] == NumGlobalElements-1) {
      Indices[0] = NumGlobalElements-2;
      NumEntries = 1;
    } else {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      NumEntries = 2;
    }
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
    // Put in the diagonal entry
    assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i)==0);
  }
  
  // Finish up
  assert(A.TransformToLocal()==0);

  // E N D   O F   M A T R I X   C O N S T R U C T I O N

  Epetra_Vector x(Map);
  Epetra_Vector b(Map);

  x.Random();

  // this is the linear problem to be solved: set the linear operator,
  // the solution and the right-hand side
  Epetra_LinearProblem A_Problem(&A, &x, &b);

  // and this is the AztecOO solver
  AztecOO A_Solver(A_Problem);

  // --- Here we define the precondioner ---
  // create P as a copy of A, in principle can be different
  Epetra_CrsMatrix P(A);

  // Here we create the linear problem which will be used as
  // preconditioner. This requires sereval steps.
  // (Note that all the P_ prefix indentify preconditioner'
  // objects)

  // 1. we create the linear system solve at each prec step
  Epetra_LinearProblem P_Problem;
  // and we assign the linear operator (in this case, the
  // matrix A itself)
  P_Problem.SetOperator(&P);

  // as we wish to use AztecOO to solve the prec step
  // (in a recursive way), we have to define an AztecOO
  // object. 
  AztecOO P_Solver(P_Problem);

  // now, we customize certain parameters
  P_Solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  P_Solver.SetAztecOption(AZ_output, AZ_none);
  P_Solver.SetAztecOption(AZ_solver, AZ_cg);
  // The last step is to create an AztecOO_Operator, so that
  // we can set the Aztec's preconditioner with.
  AztecOO_Operator P_Operator(&P_Solver, 10);

  // Here we set the user's defined preconditioners
  A_Solver.SetPrecOperator(&P_Operator);

  // --- up to here ---
  
  // Finally, we solve the linear system:
  
  int Niters=100;
  A_Solver.SetAztecOption(AZ_kspace, Niters);
  A_Solver.SetAztecOption(AZ_solver, AZ_GMRESR);

  A_Solver.Iterate(Niters, 1.0E-12);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0;
  
}

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:aztecoo]> mpirun -np 2 ./ex2

                *******************************************************
                ***** Preconditioned GMRESR solution
                ***** AztecOO Operator
                ***** No scaling
                *******************************************************

                iter:    0           residual = 1.000000e+00
                iter:    1           residual = 4.809435e-17


                Solution time: 0.002579 (sec.)
                total iterations: 1

*/		
