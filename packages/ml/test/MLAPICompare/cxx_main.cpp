
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
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  Epetra_Time Time(Comm);

  int ProblemSize = 100;
  if (argc > 1) {
    ProblemSize = atoi(argv[1]);
  }
    
  int NumPDEEqns = 4;

  VbrMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", ProblemSize);
  Epetra_RowMatrix* A = Gallery.GetVbrMatrix(NumPDEEqns);
  Epetra_LinearProblem* Problem = Gallery.GetVbrLinearProblem();
  Epetra_MultiVector* LHS = Gallery.GetVbrStartingSolution();
  Epetra_MultiVector* RHS = Gallery.GetVbrRHS();

  AztecOO solver(*Problem);

  if (argc == 3) {

    int rep = 100;
    double res_ML, res_MLAPI;

    // want to compute r = b - A * x
    Epetra_Vector x_ML(A->OperatorDomainMap());
    Epetra_Vector b_ML(A->OperatorRangeMap());
    Epetra_Vector r_ML(A->OperatorRangeMap());

    // test the Epetra normal way
    Time.ResetStartTime();
    for (int i = 0 ; i < rep ; ++i) {

      x_ML.PutScalar(1.0);
      b_ML.PutScalar(1.0);
      A->Multiply(false,x_ML,r_ML);
      r_ML.Update(1.0, b_ML, -1.0);

    }
    res_ML = Time.ElapsedTime();
    cout << "time - ML    = " << res_ML / rep << endl;

    // now the MAPI way
    Init();
    Space S(-1,A->NumMyRows());
    Operator A_MLAPI(S,S,A,false);
    MultiVector x_MLAPI(S), b_MLAPI(S), r_MLAPI(S);

    Time.ResetStartTime();
    for (int i = 0 ; i < rep ; ++i) {

      x_MLAPI = 1.0;
      b_MLAPI = 1.0;

#if 1
      A_MLAPI.Apply(x_MLAPI, r_MLAPI);
      r_MLAPI.Update(1.0, b_MLAPI, -1.0);
#else
      r_MLAPI = b_MLAPI - A_MLAPI * x_MLAPI;
#endif

    }
    res_MLAPI = Time.ElapsedTime();
    cout << "time - MLAPI = " << res_MLAPI / rep << endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    exit(0);
  }

  // =========================== parameters =================================

  double DampingFactor = 1.333;
  int    OutputLevel = 16;

  int MLPIters, MLAPIIters;
  double MLPResidual, MLAPIResidual;
  double MLPConstructionTime, MLAPIConstructionTime;
  double MLPSolveTime, MLAPISolveTime;

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
  MLList.set("PDE equations", NumPDEEqns);
  
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
  solver.SetAztecOption(AZ_output, OutputLevel);
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
  solver.SetAztecOption(AZ_output, OutputLevel);
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
