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

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_MatrixFreePreconditioner.h"
#include "ml_epetra_utils.h"

using namespace ML_Epetra;
using namespace Teuchos;
using namespace Galeri;

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

  int nx;
  if (argc > 1)
    nx = (int) strtol(argv[1],NULL,10);
  else
    nx = 8;
  int ny = nx * Comm.NumProc(); // each subdomain is a square

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());
  GaleriList.set("conv", 1e0); // for Recirc2D
  GaleriList.set("diff", 1e0); // for Recirc2D

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* CrsA = CreateCrsMatrix("Laplace2D", Map, GaleriList);
  int NumPDEEqns = 2;
  int NullSpaceDim = 3;
  Epetra_VbrMatrix* VbrA = CreateVbrMatrix(CrsA, NumPDEEqns);
  delete CrsA;

  Epetra_MultiVector OrigNullSpace(VbrA->Map(), NullSpaceDim);
  OrigNullSpace.Random();
  Epetra_MultiVector NullSpace(OrigNullSpace);
  Teuchos::ParameterList MLList;
  MLList.set("ML output", 0);
  MLList.set("prec: type", "hybrid");
  MLList.set("smoother: type", "Chebyshev");
  MLList.set("low memory", true);
  MLList.sublist("ML list").set("ML output", 0);
  MLList.sublist("ML list").set("max levels", 10);
  MLList.sublist("ML list").set("smoother: type", "Aztec");
  MLList.sublist("ML list").set("aggregation: damping factor", 0.0);
  MLList.sublist("ML list").set("smoother: pre or post", "both");
  MLList.sublist("ML list").set("coarse: max size", 32);

  Epetra_Vector PointDiagonal(VbrA->Map());
  VbrA->ExtractDiagonalCopy(PointDiagonal);

  // compute the preconditioner using the matrix-free approach
  MatrixFreePreconditioner* MFP = new
    MatrixFreePreconditioner(*VbrA, VbrA->Graph(), NullSpace,
                             PointDiagonal, MLList);

  NullSpace = OrigNullSpace;

  assert (MFP->IsComputed() == true);

  /* check to verify that all the objects are SPD
  int NumChecks = 10;
  MFP->CheckSPD(*VbrA, true, NumChecks);
  MFP->CheckSPD(MFP->C(), true, NumChecks);
  MFP->CheckSPD(MFP->MLP(), false, NumChecks);
  */

  // =========== //
  // CHECK START //
  // =========== //
  
  // wrap VbrA as an ML_Operator
  ML_Operator* VbrA_ML = ML_Operator_Create(MFP->Comm_ML());
  ML_Operator_WrapEpetraMatrix(VbrA, VbrA_ML);
  // then build P, R and C using ML, based on VbrA and default null space
  ML_Aggregate* Aggregates_ML;
  ML_Operator* P_ML,* R_ML,* C_ML;
  MFP->Coarsen(VbrA_ML, &Aggregates_ML, &P_ML, &R_ML, &C_ML, NumPDEEqns,
               NullSpace.NumVectors(), NullSpace.Values());
  ML_Aggregate_Destroy(&Aggregates_ML);

  // ========================================== //
  // CHECK 1: Non-smoothed Restriction Operator //
  // ========================================== //
  
  Epetra_CrsMatrix* R;
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(R_ML, R));

  assert (R->OperatorDomainMap().NumGlobalElements() == 
          MFP->R().OperatorDomainMap().NumGlobalElements());
  assert (R->OperatorRangeMap().NumGlobalElements() == 
          MFP->R().OperatorRangeMap().NumGlobalElements());

  //cout << *R;
  //cout << MFP->R();

  Epetra_Vector x(R->OperatorDomainMap());
  Epetra_Vector y(R->OperatorRangeMap());
  Epetra_Vector z(R->OperatorRangeMap());

  for (int t = 0; t < 2; ++t)
  {
    x.Random();
    R->Apply(x, y);
    MFP->R().Apply(x, z);
    y.Update(1.0, z, -1.0);
    double norm;
    y.Norm2(&norm);

    if (Comm.MyPID() == 0)
      cout << "test " << t << ", ||(R_ML - R_MFP) * y||_2 = " << norm << endl;

    if (norm > 1e-5) exit(EXIT_FAILURE);
  }

  // =========================================== //
  // CHECK 2: Coarse-level operator.             //
  // aggregate and build C using ML on VbrA, and //
  // check that the coarse operator is the same  //
  // =========================================== //

  Epetra_CrsMatrix* C;
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(C_ML, C));

  assert (C->OperatorRangeMap().SameAs(MFP->C().OperatorRangeMap()));
  assert (C->OperatorDomainMap().SameAs(MFP->C().OperatorDomainMap()));

  Epetra_Vector xx(C->OperatorDomainMap());
  Epetra_Vector yy(C->OperatorRangeMap());
  Epetra_Vector zz(C->OperatorRangeMap());

  for (int t = 0; t < 2; ++t)
  {
    xx.Random();
    C->Apply(xx, yy);
    MFP->C().Apply(xx, zz);
    yy.Update(1.0, zz, -1.0);
    double norm;
    yy.Norm2(&norm);

    if (Comm.MyPID() == 0)
      cout << "test " << t << ", ||(C_ML - C_MFP) * y||_2 = " << norm << endl;
    
    if (norm > 1e-5) exit(EXIT_FAILURE);
  }

  delete C;
  ML_Operator_Destroy(&P_ML);
  ML_Operator_Destroy(&R_ML);
  ML_Operator_Destroy(&C_ML);
  ML_Operator_Destroy(&VbrA_ML);

  // ================= //
  // CHECK 3: solution //
  // ================= //

  Epetra_Vector LHS(VbrA->OperatorDomainMap()); 
  Epetra_Vector RHS(VbrA->OperatorRangeMap()); 
  Epetra_LinearProblem Problem(VbrA, &LHS, &RHS);
  AztecOO solver(Problem);

  LHS.PutScalar(1.0);
  RHS.PutScalar(0.0);

  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.SetPrecOperator(MFP);
  solver.Iterate(500, 1e-5);

  int MFPIters = solver.NumIters();

  delete MFP;

  // compare with the same stuff without low memory. All 
  // parameters are exactly as before, except "low memory".

  NullSpace = OrigNullSpace;

  MLList.set("low memory", false);
  MatrixFreePreconditioner* MFP_noLowMemory = new
    MatrixFreePreconditioner(*VbrA, VbrA->Graph(), NullSpace,
                             PointDiagonal, MLList);

  assert (MFP_noLowMemory->IsComputed() == true);

  LHS.PutScalar(1.0);
  RHS.PutScalar(0.0);

  solver.SetPrecOperator(MFP_noLowMemory);
  solver.Iterate(500, 1e-5);

  int MFP2Iters = solver.NumIters();

  assert (MFP2Iters == MFPIters);

  delete MFP_noLowMemory;

  LHS.PutScalar(1.0);
  RHS.PutScalar(0.0);

  Teuchos::ParameterList MLList2;
  MLList2.set("PDE equations", NumPDEEqns);
  MLList2.set("ML output", 0);
  MLList2.set("max levels", 10);
  MLList2.set("coarse: max size", 32);
  MLList2.set("smoother: type", "Chebyshev");
  MLList2.set("smoother: sweeps", 3);
  MLList2.set("smoother: type (level 1)", "Aztec");
  //MLList2.set("prec type", "two-level-additive");
  MLList2.set("aggregation: damping factor", 0.0);
  MLList2.set("smoother: pre or post", "both");
  MLList2.set("null space: type", "pre-computed");
  MLList2.set("null space: dimension", OrigNullSpace.NumVectors());
  MLList2.set("null space: vectors", OrigNullSpace.Values());

  MultiLevelPreconditioner* MLP = new MultiLevelPreconditioner(*VbrA, MLList2);

  solver.SetPrecOperator(MLP);
  solver.Iterate(500, 1e-5);

  int MLPIters = solver.NumIters();

  delete MLP;

  assert (abs(MLPIters < MFPIters) < 2);

  delete VbrA;
  delete Map;

  if (Comm.MyPID() == 0)
    cout << "TEST PASSED" << endl;

#ifdef HAVE_MPI
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
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with:");
  puts("--enable-epetra");
  puts("--enable-epetraext");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");
  puts("--enable-galeri");
  puts("--enable-ifpack");
  puts("--enable-mpi");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  exit(EXIT_SUCCESS);
}

#endif //#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)
