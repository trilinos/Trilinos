
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
// Solve a linear system with KLU
//
// NOTE: The Makefile may change, depending on how Amesos has been configured.
//       Please check linker flags

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Amesos.h"
#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // initialize the command line parser
  Trilinos_Util_CommandLineParser CLP(argc,argv);

  // initialize an Gallery object
  Trilinos_Util_CrsMatrixGallery Gallery("", Comm);

  // add default values
  if( CLP.Has("-problem_type") == false ) CLP.Add("-problem_type", "laplace_2d" ); 
  if( CLP.Has("-problem_size") == false ) CLP.Add("-problem_size", "100" ); 

  // initialize the gallery as specified in the command line
  Gallery.Set(CLP);

  // retrive pointers to matrix, StartingSolution, RHS, and linear problem
  // NOTE: StartingSolution and RHS are pointers to Gallery's internally stored
  // vectors. Using StartingSolution and RHS, we can verify the residual
  // after the solution of the linear system. However, users may define as well
  // their own vectors for solution and RHS. 
  Epetra_CrsMatrix * Matrix = Gallery.GetMatrix();
  Epetra_Vector * LHS = Gallery.GetStartingSolution();
  Epetra_Vector * RHS = Gallery.GetRHS();
  
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  // initialize Amesos solver  
  Amesos_BaseSolver * Solver;
  Amesos Amesos_Factory;

  // empty parameter list
  Teuchos::ParameterList List;
  
  // use KLU (other choices are valid as well)
  int choice = 0;

  switch( choice ) {
  case 0:
    Solver = A_Factory.Create("Amesos_Klu", *Problem);
    break;

  case 1:
    Solver = A_Factory.Create("Amesos_Umfpack", *Problem);
    break;

  }
  
  // start solving
  Solver->SymbolicFactorization();
  Solver->NumericFactorization();
  Solver->Solve();

  // verify that residual is really small  
  double residual, diff;

  Gallery.ComputeResidual(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(diff);

  if( Comm.MyPID() == 0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

  // delete Solver
  delete Solver;
    
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}
