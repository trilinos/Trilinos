// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teuchos includes /*@ \label{lned:being-includes} @*/
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Epetra includes
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StridedEpetraOperator.hpp"
#include "Teko_InverseFactoryOperator.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

#include <iostream> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char* argv[]) {
  // calls MPI_Init and MPI_Finalize
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // read in parameter list
  Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::getParametersFromXmlFile("simple_example.xml");

  // build global communicator
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Read in the matrix, store pointer as an RCP
  Epetra_CrsMatrix* ptrA = 0;
  EpetraExt::MatrixMarketFileToCrsMatrix("../data/nsjac.mm", Comm, ptrA);
  RCP<Epetra_CrsMatrix> A = rcp(ptrA);

  // read in the RHS vector
  Epetra_Vector* ptrb = 0;
  EpetraExt::MatrixMarketFileToVector("../data/nsrhs_test.mm", A->OperatorRangeMap(), ptrb);
  RCP<Epetra_Vector> b = rcp(ptrb);

  // allocate vectors
  RCP<Epetra_Vector> x = rcp(new Epetra_Vector(A->OperatorDomainMap()));
  x->PutScalar(0.0);

  // Break apart the strided linear system
  /////////////////////////////////////////////////////////

  // Block the linear system using a strided epetra operator
  std::vector<int> vec(2);
  vec[0] = 2;
  vec[1] = 1; /*@ \label{lned:define-strided} @*/
  Teuchos::RCP<const Epetra_Operator> strided_A =
      Teuchos::rcp(new Teko::Epetra::StridedEpetraOperator(vec, A));

  // Build the preconditioner /*@ \label{lned:construct-prec} @*/
  /////////////////////////////////////////////////////////

  // build an InverseLibrary and inverse factory
  RCP<Teko::InverseLibrary> invLib  = Teko::InverseLibrary::buildFromParameterList(*paramList);
  RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory("SIMPLE");

  // Create the initial preconditioner, and build it from strided_A
  RCP<Teko::Epetra::InverseFactoryOperator> prec_A =
      rcp(new Teko::Epetra::InverseFactoryOperator(inverse));
  prec_A->initInverse();
  prec_A->rebuildInverseOperator(strided_A.getConst());

  // Build and solve the linear system
  /////////////////////////////////////////////////////////

  // Setup the linear solve: notice A is used directly
  Epetra_LinearProblem problem(&*A, &*x, &*b); /*@ \label{lned:aztec-solve} @*/

  // build the solver
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_precond, AZ_none);
  solver.SetAztecOption(AZ_kspace, 1000);
  solver.SetAztecOption(AZ_output, 10);
  solver.SetPrecOperator(&*prec_A);

  // solve the linear system
  solver.Iterate(1000, 1e-5);

  return 0;
}
