
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
// 
// Use of ML as a black-box smoothed aggregation preconditioner
//
// NOTE: The linker line in the Makefile can change considerably, depending on how
// you configures ML and Amesos (if ML was configured with Amesos support). You may need to add
// -lteuchos -lamesos -lifpack -lparmetis-3.1 -ly12m -lumfpack -lamd -lmetis-4.0
// to the Makefile

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"

// includes required by ML
#include "ml_epetra_preconditioner.h"

#include "Trilinos_Util_CommandLineParser.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Teuchos;

#include <iostream>

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  Epetra_Time Time(Comm);

  // initialize the command line parser
  Trilinos_Util_CommandLineParser CLP(argc,argv);

  // initialize an Gallery object
  Trilinos_Util_CrsMatrixGallery Gallery("", Comm);

  // add default values
  if( CLP.Has("-problem_type") == false ) CLP.Add("-problem_type", "laplace_2d" ); 
  if( CLP.Has("-problem_size") == false ) CLP.Add("-problem_size", "10000" ); 

  // initialize the gallery as specified in the command line
  Gallery.Set(CLP);

  // retrive pointers to matrix and linear problem
  Epetra_RowMatrix * A = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();
  
  // Construct a solver object for this problem
  AztecOO solver(*Problem);
  
  // create the preconditioner object and compute hierarchy
  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, CLP, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
 
  double rthresh = 1.4;
  solver.SetAztecParam(AZ_rthresh, rthresh);
  solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
  solver.SetAztecOption(AZ_output, 32);
  double athresh = 10.0;
  solver.SetAztecParam(AZ_athresh, athresh);
  solver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);

  int Niters = 500;
  solver.SetAztecOption(AZ_kspace, 160);
   
  solver.Iterate(Niters, 1e-12);

  // print out some information about the preconditioner
  if( Comm.MyPID() == 0 ) cout << MLPrec->GetOutputList();
  
  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
  
}
