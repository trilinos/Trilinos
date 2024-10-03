// @HEADER
// *****************************************************************************
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack.h"
#include "Ifpack_Hypre.h"

#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"

#include "Epetra_MultiVector.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_InvOperator.h"

#include "BelosEpetraAdapter.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "hypre_Helpers.hpp"

#include <string>
#include <stdio.h>
#include <map>


TEUCHOS_UNIT_TEST(Belos_Hypre, Laplace2D){
  const double tol = 1E-7;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  typedef Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>  LinearProblem;

  //
  // Create Laplace2D
  //
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm();
#endif
  Teuchos::ParameterList GaleriList;
  int nx = 10 * Comm.NumProc();
  int ny = 10 * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  Epetra_Map Map(nx*ny,0,Comm);
  RCP<Epetra_CrsMatrix> Crs_Matrix   = rcp(Galeri::CreateCrsMatrix("Laplace2D", &Map, GaleriList));
  int NumProc = Crs_Matrix->Comm().NumProc();

  //
  // Create the hypre preconditioner
  //
  RCP<Ifpack_Hypre> preconditioner = rcp(new Ifpack_Hypre(Crs_Matrix.get()));
  TEST_EQUALITY(preconditioner->Initialize(),0);
  TEST_EQUALITY(preconditioner->SetParameter(Preconditioner, ParaSails),0); // Use a Euclid Preconditioner (but not really used)
  TEST_EQUALITY(preconditioner->SetParameter(Preconditioner),0); // Solve the problem
  TEST_EQUALITY(preconditioner->Compute(),0);

  //
  // Create the solution vector and rhs
  //
  int numVec = 1;
  RCP<Epetra_MultiVector> X = rcp(new Epetra_MultiVector(Crs_Matrix->OperatorDomainMap(), numVec));
  RCP<Epetra_MultiVector> KnownX = rcp(new Epetra_MultiVector(Crs_Matrix->OperatorDomainMap(), numVec));
  KnownX->Random();
  RCP<Epetra_MultiVector> B = rcp(new Epetra_MultiVector(Crs_Matrix->OperatorRangeMap(), numVec));
  Crs_Matrix->Apply(*KnownX, *B);

  //
  // Test the EpetraExt wrapper
  // amk November 24, 2015: Should we deprecate this?
  //
//  RCP<ParameterList> pl = rcp(new ParameterList());
//  TEST_EQUALITY(X->PutScalar(0.0),0);
//  HYPRE_IJMatrix hypre_mat = preconditioner->HypreMatrix();
//  RCP<EpetraExt_HypreIJMatrix> Hyp_Matrix = rcp(new EpetraExt_HypreIJMatrix(hypre_mat));
//  TEST_EQUALITY(Hyp_Matrix->SetParameter(Preconditioner, ParaSails),0);
//  TEST_EQUALITY(Hyp_Matrix->SetParameter(Preconditioner),0);
//  TEST_EQUALITY(EquivalentMatrices(*Hyp_Matrix, *Crs_Matrix, tol), true);
//  RCP<LinearProblem> problem1 = rcp(new LinearProblem(Crs_Matrix,X,B));
//  problem1->setLeftPrec(Hyp_Matrix);
//  TEST_EQUALITY(problem1->setProblem(),true);
//  Belos::PseudoBlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator> solMgr1(problem1,pl);
//  Belos::ReturnType rv1 = solMgr1.solve(); // TEST_EQUALITY(solMgr2.solve(),Belos::Converged);
//  TEST_EQUALITY(rv1,Belos::Converged);
//  TEST_EQUALITY(EquivalentVectors(*X, *KnownX, tol*10*pow(10.0,NumProc)), true);

  //
  // Test the Ifpack hypre interface
  //
  RCP<ParameterList> pl2 = rcp(new ParameterList());
  RCP<Epetra_Operator> invOp = rcp(new Epetra_InvOperator(preconditioner.get()));
  TEST_EQUALITY(X->PutScalar(0.0),0);
  RCP<LinearProblem> problem2 = rcp(new LinearProblem(Crs_Matrix,X,B));
  problem2->setLeftPrec(invOp);
  TEST_EQUALITY(problem2->setProblem(),true);
  Belos::PseudoBlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator> solMgr2(problem2,pl2);
  Belos::ReturnType rv2 = solMgr2.solve(); // TEST_EQUALITY(solMgr2.solve(),Belos::Converged);
  TEST_EQUALITY(rv2,Belos::Converged);
  TEST_EQUALITY(EquivalentVectors(*X, *KnownX, tol*10*pow(10.0,NumProc)), true);
}
