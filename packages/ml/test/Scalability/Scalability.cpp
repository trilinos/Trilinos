
//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif
#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Trilinos_Util;

// Simple test to check basic scalability issues.
//
// This problem creates a matrix using the matrix gallery, using
// a linear partition among the processors. (Note that a linear partition is
// not an optimal partition for domain decomposition methods for some 2D and
// 3D problems). The ML preconditioner is created, and the linear problem is
// solved using AztecOO (CG).
//
// The code prints out the time required to build the ML preconditioner, to
// apply it once (this number can be negative for small problems), and the
// AztecOO time. 
//
// Users can specify the matrix name using --matrix (refer to the Trilnos
// tutorial in the Triutils chapters for more details), and the global size
// using --size.
// 
// A simple scalability test can read as follows:
//
// $ mpirun -np 1 ./Scalability.exe --matrix=laplace_2d --size=90000
// $ mpirun -np 2 ./Scalability.exe --matrix=laplace_2d --size=90000
// $ mpirun -np 4 ./Scalability.exe --matrix=laplace_2d --size=90000
// $ mpirun -np 8 ./Scalability.exe --matrix=laplace_2d --size=90000
//
// If you have a parallel machine, then the time per iteration (TPI) should
// decrease linearly. If you have a single-processor machine, then the
// TPI should remain roughly constant.
//
// Here there are some results obtained on a single-processor LINUX
// (Intel(R) Pentium(R) M processor 1700MHz, GNU compilers, LAM/MPI)
// (TPI should remain constant)
// 
// Processors      TPI
// ----------      ---
//  1               0.29618
//  2               0.301127
//  4               0.305803
//  8               0.315509
//
// Here are instead some results obtained on a LINUX cluster
// (TPI * NumProcs should remain constant)
//
// Processors      TPI
// ----------      ---
//  1               0.96114
//  2               0.48463
//  4               0.25713 
//  8               
//  
//  
// Author: Marzio Sala, SNL 9214
// Last modified: Nov-04
// 
// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Epetra_Time Time(Comm);

  double ConstructionTime, ApplicationTime, SolutionTime;

  // process the command line
  Teuchos::CommandLineProcessor CLP;
  // matrix name
  string MatrixName = "laplace_2d";
  int GlobalSize = 90000;
  CLP.setOption("matrix", &MatrixName, "Matrix name for Gallery (ex: laplace_2d)");
  CLP.setOption("size", &GlobalSize, "Size of the problem. Note: this may need to be a square/cube depending on specified matrix name");

  CLP.throwExceptions(false);
  CLP.parse(argc,argv);

  CrsMatrixGallery Gallery(MatrixName.c_str(), Comm);
  Gallery.Set("problem_size", GlobalSize);
  
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  AztecOO solver(*Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation
  ML_Epetra::SetDefaults("SA",MLList);
  
  MLList.set("output", 0);

  // create the preconditioning object. 
  Time.ResetStartTime();
  ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);
  ConstructionTime = Time.ElapsedTime();

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // apply the preconditioner once, track required time
  {
    Epetra_MultiVector FakeX(*(Problem->GetLHS()));
    Epetra_MultiVector FakeY(*(Problem->GetRHS()));

    Time.ResetStartTime();
    MLPrec->Apply(FakeX, FakeY);
    ApplicationTime = Time.ElapsedTime();
  }
  
  // =========================== end of ML part =============================
  
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 32);

  Time.ResetStartTime();
  solver.Iterate(500, 1e-12);
  SolutionTime = Time.ElapsedTime();

  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);
  
  if (Comm.MyPID() == 0) {
    string Prefix = "[ML test] ";
    cout << Prefix << "# of processors  = " << Comm.NumProc() << endl;
    cout << Prefix << "# of global rows = " << GlobalSize << endl;
    cout << Prefix << "ML construction time    = " << ConstructionTime << " (s)" << endl;
    cout << Prefix << "one ML application time = " << ApplicationTime << " (s)" << endl;
    cout << Prefix << "AztecOO solution time   = " << SolutionTime << " (s)" << endl;
    cout << Prefix << "AztecOO iterations      = " << solver.NumIters() << endl;
    cout << endl;
    cout << Prefix << "time per iteration      = " 
         << SolutionTime / solver.NumIters() << " (s)" << endl;
    cout << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(0);
  
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");
  
  return(0);
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
