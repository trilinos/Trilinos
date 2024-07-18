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
#include "Teko_StratimikosFactory.hpp"

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_VectorBase.hpp"

#include <iostream> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char* argv[]) {
  // calls MPI_Init and MPI_Finalize
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // read in parameter list
  Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::getParametersFromXmlFile("strat_example.xml");

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
  RCP<const Epetra_Vector> b = rcp(ptrb);

  // allocate vectors
  RCP<Epetra_Vector> x = rcp(new Epetra_Vector(A->OperatorDomainMap()));
  x->PutScalar(0.0);

  // Build Thyra linear algebra objects
  RCP<const Thyra::LinearOpBase<double> > th_A = Thyra::epetraLinearOp(A);
  RCP<const Thyra::VectorBase<double> > th_b   = Thyra::create_Vector(b, th_A->range());
  RCP<Thyra::VectorBase<double> > th_x         = Thyra::create_Vector(x, th_A->domain());

  // Build stratimikos solver
  /////////////////////////////////////////////////////////

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

  Teko::addTekoToStratimikosBuilder(linearSolverBuilder);
  linearSolverBuilder.setParameterList(paramList);

  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
      Thyra::createLinearSolveStrategy(linearSolverBuilder);

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > th_invA =
      Thyra::linearOpWithSolve(*lowsFactory, th_A);

  Thyra::assign(th_x.ptr(), 0.0);
  Thyra::SolveStatus<double> status =
      Thyra::solve<double>(*th_invA, Thyra::NOTRANS, *th_b, th_x.ptr());

  return 0;
}
