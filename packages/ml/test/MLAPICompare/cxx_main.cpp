
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

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

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
#include "AztecOO.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
// includes required by ML

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiLevelSA.h"
#include "MLAPI_MATLABStream.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_EpetraBaseOperator.h"

using namespace Teuchos;
using namespace Trilinos_Util;
using namespace ML_Epetra;
using namespace MLAPI;

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

  int ProblemSize = 10000;
  if (argc == 2) {
    ProblemSize = atoi(argv[1]);
  }
    
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", ProblemSize);
  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  Epetra_MultiVector* LHS = Gallery.GetStartingSolution();
  Epetra_MultiVector* RHS = Gallery.GetRHS();

  AztecOO solver(*Problem);

  // =========================== parameters =================================
  
  double DampingFactor = 1.333;
  int    Output = 16;

  int MLPIters, MLAPIIters;
  double MLPResidual, MLAPIResidual;
  double MLPConstructionTime, MLAPIConstructionTime;
  double MLPSolveTime, MLAPISolveTime;
  Epetra_Time Time(Comm);
    
  // =========================== begin of ML part ===========================
  
  ParameterList MLList;
  ParameterList IFPACKList;
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("max levels",10);
  MLList.set("increasing or decreasing","increasing");
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("aggregation: damping factor", DampingFactor); 
  MLList.set("coarse: max size",32);
  MLList.set("smoother: pre or post", "both");
  MLList.set("coarse: type","Amesos-KLU");
  MLList.set("smoother: type","IFPACK");
  MLList.set("smoother: ifpack type","point relaxation stand-alone");
  MLList.set("smoother: sweeps",1);
  MLList.set("smoother: damping factor",0.67);
  IFPACKList.set("relaxation: sweeps", 1);
  IFPACKList.set("relaxation: damping factor", 0.67);
  IFPACKList.set("relaxation: type", "symmetric Gauss-Seidel");
  IFPACKList.set("relaxation: zero starting solution", false);
  MLList.set("smoother: ifpack list", IFPACKList);
  
  // ======================================================== //
  // build the MultiLevelPreconditioner first, and track down //
  // the number of iterations, the residual and the CPU-time  //
  // ======================================================== //
  
  Time.ResetStartTime();

  MultiLevelPreconditioner* MLPPrec;
  MLPPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);
  MLPConstructionTime = Time.ElapsedTime();
  Time.ResetStartTime();

  LHS->PutScalar(0.0);
  RHS->PutScalar(1.0);

  solver.SetPrecOperator(MLPPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 16);
  solver.Iterate(1550, 1e-5);

  MLPIters = solver.NumIters();
  MLPResidual = solver.TrueResidual();

  delete MLPPrec;
  MLPSolveTime = Time.ElapsedTime();

  // ======================================================= //
  // now we aim to define the same preconditioner, using the //
  // MLAPI interface. Iterations and residual should be the  //
  // same, then we compare the CPU time as well.             //
  // ======================================================= //
  
  MLList.set("smoother: type", "symmetric Gauss-Seidel");

  Init();
  int size = A->NumMyRows();
  Space FineSpace(-1,size);

  Operator AA(FineSpace,FineSpace,A,false);

  Time.ResetStartTime();

  MultiLevelSA* Cycle = new MultiLevelSA(AA,MLList);
  Epetra_Operator* MLAPIPrec = 
    new EpetraBaseOperator(A->RowMatrixRowMap(),*Cycle);
  MLAPIConstructionTime = Time.ElapsedTime();
  Time.ResetStartTime();

  LHS->PutScalar(0.0);
  RHS->PutScalar(1.0);

  solver.SetPrecOperator(MLAPIPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 16);
  solver.Iterate(1550, 1e-5);

  MLAPIIters = solver.NumIters();
  MLAPIResidual = solver.TrueResidual();

  delete MLAPIPrec;
  MLAPISolveTime = Time.ElapsedTime();

  // ================ //
  // output some crap //
  // ================ //

  if( Comm.MyPID()==0 ) {
    cout << "Iterations: " << MLPIters    << " vs. " << MLAPIIters    << endl;
    cout << "Residual  : " << MLPResidual << " vs. " << MLAPIResidual << endl;
    cout << "CPU-time (const) : " << MLPConstructionTime     << " vs. " 
         << MLAPIConstructionTime     << endl;
    cout << "CPU-time (solve) : " << MLPSolveTime     << " vs. " 
         << MLAPISolveTime     << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
  
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
