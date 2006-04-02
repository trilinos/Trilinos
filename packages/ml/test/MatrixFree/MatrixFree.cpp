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
    nx = 4;
  int ny = nx * Comm.NumProc(); // each subdomain is a square

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* CrsA = CreateCrsMatrix("Laplace2D", Map, GaleriList);
  int NumPDEEqns = 2;
  int NullSpaceDim = 3;
  Epetra_VbrMatrix* VbrA = CreateVbrMatrix(CrsA, NumPDEEqns);

  Epetra_MultiVector OrigNullSpace(VbrA->Map(), NullSpaceDim);
  //OrigNullSpace.PutScalar(1.0);
  OrigNullSpace.Random();
  Epetra_MultiVector NullSpace(OrigNullSpace);
  Teuchos::ParameterList MLList;
  MLList.sublist("ML list").set("max levels", 1);
  //MLList.sublist("ML list").set("max levels", 2);
  //MLList.sublist("ML list").set("coarse: max size", 1);
  //MLList.sublist("ML list").set("cycle applications", 10);

  // compute the preconditioner using the matrix-free approach
  MatrixFreePreconditioner* MFP = new
    MatrixFreePreconditioner(*VbrA, VbrA->Graph(), MLList, NullSpace);

  NullSpace = OrigNullSpace;

  assert (MFP->IsComputed() == true);

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

  for (int t = 0; t < 3; ++t)
  {
    x.Random();
    R->Apply(x, y);
    MFP->R().Apply(x, z);
    y.Update(1.0, z, -1.0);
    double norm;
    y.Norm2(&norm);

    if (Comm.MyPID() == 0)
      cout << "test " << t << ", ||(R_ML - R_MFP) * y||_2 = " << norm << endl;

    if (norm > 1e-10) exit(EXIT_FAILURE);
  }

  // =========================================== //
  // CHECK 2: Coarse-level operator.             //
  // aggregate and build C using ML on VbrA, and //
  // check that the coarse operator is the same  //
  // =========================================== //

  Epetra_CrsMatrix* C;
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(C_ML, C));
  
  // cout << *C;
  // cout << MFP->C();

  assert (C->OperatorRangeMap().SameAs(MFP->C().OperatorRangeMap()));
  assert (C->OperatorDomainMap().SameAs(MFP->C().OperatorDomainMap()));

  Epetra_Vector xx(C->OperatorDomainMap());
  Epetra_Vector yy(C->OperatorRangeMap());
  Epetra_Vector zz(C->OperatorRangeMap());

  for (int t = 0; t < 3; ++t)
  {
    xx.Random();
    C->Apply(xx, yy);
    MFP->C().Apply(xx, zz);
    yy.Update(1.0, zz, -1.0);
    double norm;
    yy.Norm2(&norm);

    if (Comm.MyPID() == 0)
      cout << "test " << t << ", ||(C_ML - C_MFP) * y||_2 = " << norm << endl;
    
    if (norm > 1e-10) exit(EXIT_FAILURE);
  }

  // ================= //
  // CHECK 3: solution //
  // ================= //

  Teuchos::ParameterList MLList2;
  MLList2.set("PDE equations", NumPDEEqns);
  MLList2.set("max levels", 2);
  MLList2.set("coarse: max size", 4);
  MLList2.set("smoother: type", "Jacobi");
  MLList2.set("prec type", "two-level-additive");
  MLList2.set("aggregation: damping factor", 0.0);
  MLList2.set("smoother: pre or post", "post");
  MLList2.set("null space: type", "pre-computed");
  MLList2.set("null space: dimension", OrigNullSpace.NumVectors());
  MLList2.set("null space: vectors", OrigNullSpace.Values());

  MultiLevelPreconditioner* MLP = new MultiLevelPreconditioner(*VbrA, MLList2);

  Epetra_Vector LHS(VbrA->OperatorDomainMap()); 
  Epetra_Vector RHS(VbrA->OperatorRangeMap()); 
  Epetra_LinearProblem Problem(VbrA, &LHS, &RHS);
  AztecOO solver(Problem);

  LHS.PutScalar(1.0);
  RHS.PutScalar(0.0);

  solver.SetPrecOperator(MLP);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-5);

  int MLPIters = solver.NumIters();
  double MLPResidual = solver.TrueResidual();

  delete MLP;

  LHS.PutScalar(1.0);
  RHS.PutScalar(0.0);

  solver.SetPrecOperator(MFP);
  solver.Iterate(500, 1e-5);

  int MFPIters = solver.NumIters();
  double MFPResidual = solver.TrueResidual();

  //assert (MLPIters == MFPIters);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
