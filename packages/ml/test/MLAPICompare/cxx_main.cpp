
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
#include "ml_common.h"

#if defined(HAVE_ML_MLAPI) && defined(HAVE_ML_GALERI)

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
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
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
using namespace Galeri;
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

  int ProblemSize = 50;
  if (argc > 1) {
    ProblemSize = atoi(argv[1]);
  }
    
  Init();
  
  ParameterList GalerList;
  GalerList.set("nx", ProblemSize);
  GalerList.set("ny", ProblemSize * Comm.NumProc());
  GalerList.set("mx", 1);
  GalerList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GalerList);
  Epetra_CrsMatrix* A_Epetra = CreateCrsMatrix("Laplace2D", Map, GalerList);
  Epetra_MultiVector X_Epetra(*Map, 1);
  Epetra_MultiVector B_Epetra(*Map, 1);
  Epetra_MultiVector R_Epetra(X_Epetra);
  R_Epetra.PutScalar(0.0);

  Epetra_LinearProblem Problem(A_Epetra, &X_Epetra, &B_Epetra);

  AztecOO solver(Problem);

  Space S(A_Epetra->OperatorDomainMap());
  Operator A_MLAPI(S, S, A_Epetra, false);

  int rep = 3;
  int out_rep = 2;

  double time_Epetra = 100000.0;
  double time_MLAPI  = 100000.0;
  double norm_Epetra, norm_MLAPI;

  // ============== //
  // CHECK RESIDUAL //
  // ============== //
  
  for (int ii = 0 ; ii < out_rep ; ++ii)
  {
    // =========== //
    // EPETRA PART //
    // =========== //

    X_Epetra.PutScalar(1.0);
    B_Epetra.PutScalar(1.0);
    R_Epetra.PutScalar(0.0);

    Time.ResetStartTime();
    for (int i = 0 ; i < rep ; ++i) 
    {
      A_Epetra->Multiply(false, X_Epetra, R_Epetra);
      R_Epetra.Update(1.0, B_Epetra, -1.0);

    }
    time_Epetra = EPETRA_MIN(time_Epetra, Time.ElapsedTime());
    R_Epetra.Norm2(&norm_Epetra);

    // ========== //
    // MLAPI PART //
    // ========== //

    MultiVector X_MLAPI(S); X_MLAPI = 1.0;
    MultiVector B_MLAPI(S); B_MLAPI = 1.0;
    MultiVector R_MLAPI(S); R_MLAPI = 0.0;

    Time.ResetStartTime();
    for (int i = 0 ; i < rep ; ++i) 
    {
      R_MLAPI = B_MLAPI - A_MLAPI * X_MLAPI;

    }
    time_MLAPI = EPETRA_MIN(time_MLAPI, Time.ElapsedTime());
    norm_MLAPI = R_MLAPI.Norm2();

    if (abs(norm_Epetra - norm_MLAPI) > 1e-8)
      exit(EXIT_FAILURE);
  }

  time_Epetra /= out_rep;
  time_MLAPI /= out_rep;

  if (GetMyPID() == 0)
  {
    cout << "Checking residual..." << endl;
    cout << "Epetra = " << time_Epetra << endl;
    cout << "MLAPI  = " << time_MLAPI << endl;
    cout << "diff = " << (time_MLAPI - time_Epetra ) / time_Epetra * 100 << endl;
    cout << endl;
  }


  // =================== //
  // TEST PRECONDITIONER //
  // =================== //
  
  // =========================== parameters =================================

  double DampingFactor = 1.333;
  int    OutputLevel = 16;

  int MLPIters, MLAPIIters;
  double MLPResidual, MLAPIResidual;
  double MLPConstructionTime, MLAPIConstructionTime;
  double MLPSolveTime, MLAPISolveTime;
  int MaxIters = 1550;

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
  //IFPACKList.set("relaxation: zero starting solution", false);
  MLList.set("smoother: ifpack list", IFPACKList);
  MLList.set("eigen-analysis: type", "Anorm");
  
  MLList.set("output", 0);
  MLList.set("energy minimization: enable", false);
  MLList.set("energy minimization: type", 16);

  // ======================================================== //
  // build the MultiLevelPreconditioner first, and track down //
  // the number of iterations, the residual and the CPU-time  //
  // ======================================================== //
  
  Time.ResetStartTime();

  MultiLevelPreconditioner* MLPPrec;
  MLPPrec = new ML_Epetra::MultiLevelPreconditioner(*A_Epetra, MLList, true);
  MLPConstructionTime = Time.ElapsedTime();
  Time.ResetStartTime();

  X_Epetra.PutScalar(0.0);
  B_Epetra.PutScalar(1.0);

  solver.SetPrecOperator(MLPPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, OutputLevel);
  solver.Iterate(MaxIters, 1e-5);

  MLPIters = solver.NumIters();
  MLPResidual = solver.TrueResidual();

  delete MLPPrec;
  MLPSolveTime = Time.ElapsedTime();

  // ======================================================= //
  // now we aim to define the same preconditioner, using the //
  // MLAPI interface. Iterations and residual should be the  //
  // same, then we compare the CPU time as well.             //
  // ======================================================= //
  
  SetPrintLevel(0);
  MLList.set("smoother: type", "symmetric Gauss-Seidel");
  MLList.set("output level", 0);

  Time.ResetStartTime();

  MultiLevelSA* Cycle = new MultiLevelSA(A_MLAPI,MLList);
  MLAPIConstructionTime = Time.ElapsedTime();

  Epetra_Operator* MLAPIPrec = 
    new EpetraBaseOperator(A_Epetra->RowMatrixRowMap(),*Cycle);
  Time.ResetStartTime();

  X_Epetra.PutScalar(0.0);
  B_Epetra.PutScalar(1.0);

  solver.SetPrecOperator(MLAPIPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, OutputLevel);
  solver.Iterate(MaxIters, 1e-5);

  MLAPIIters = solver.NumIters();
  MLAPIResidual = solver.TrueResidual();

  delete MLAPIPrec;
  MLAPISolveTime = Time.ElapsedTime();

  if (Comm.MyPID() == 0) 
  {
    cout << "Iterations: " << MLPIters    << " vs. " << MLAPIIters    << endl;
    cout << "Residual  : " << MLPResidual << " vs. " << MLAPIResidual << endl;
    cout << "CPU-time (const) : " << MLPConstructionTime     << " vs. " 
         << MLAPIConstructionTime    
         << ", diff = " << (MLAPIConstructionTime - MLPConstructionTime) / MLPConstructionTime  * 100 << endl;
    cout << "CPU-time (solve) : " << MLPSolveTime     << " vs. " 
         << MLAPISolveTime    
         << ", diff = " << (MLAPISolveTime - MLPSolveTime) / MLPSolveTime * 100 << endl;
    double MLPTotal = MLPConstructionTime + MLPSolveTime;
    double MLAPITotal = MLAPIConstructionTime + MLAPISolveTime;
    cout << "CPU-time (total) : " << MLPTotal     << " vs. " 
         << MLAPITotal
         << ", diff = " << (MLAPITotal - MLPTotal) / MLPTotal * 100 << endl;
  }

  delete A_Epetra;
  delete Map;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-galeri --enable-amesos --enable-ifpack");
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}
#endif
