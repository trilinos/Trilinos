
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
// Solve a linear system with AztecOO. Linear system is created using MatrixGallery
//
// NOTE: if AztecOO has been configured with --enable-aztecoo-azlu, then you need to
// specify location and name of Y12M library in linking phase

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "AztecOO.h"
#include "Epetra_CrsMatrix.h"

#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"


int main(int argc, char *argv[])
{
    
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
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

  // retrive pointers to matrix and linear problem. LHS and RHS are
  // created an initialized by Gallery 
  Epetra_CrsMatrix * Matrix = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem= Gallery.GetLinearProblem();

  // initialize the AztecOO solve object, based on current linear problem
  AztecOO solver(*Problem);

  // here set some AztecOO options:
  // - symmetric problem;
  // - domain decomposition preconditioner
  // - ICC factorization on each subdomain
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_overlap,0);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_icc);

  // solve the linear system
  solver.Iterate(1550, 1e-12);

  // AztecOO defined a certain number of output parameters, and store them
  // in a double vector called status. 
  double status[AZ_STATUS_SIZE];
  solver.GetAllAztecStatus(status);
  
  // verify that linear system has been solved as required
  double residual, diff;

  Gallery.ComputeResidual(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }
  
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return 0 ;

}
