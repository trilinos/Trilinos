// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teuchos includes /*@ \label{lned:being-includes} @*/
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Tpetra includes
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Teko-Package includes
#include "Teko_StratimikosFactory.hpp"
#include "Teko_ConfigDefs.hpp"

// Stratimikos includes
#include "Stratimikos_MueLuHelpers.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"

#include <iostream> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

namespace{
Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> create_linear_solver_builder()
{
  Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> strat =
      Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder("",
          "",
          "",
          "linear-solver-params-file",
          "extra-linear-solver-params",
          "linear-solver-params-used-file"));
  Stratimikos::enableMueLu<double, Teko::LO, Teko::GO, Teko::NT>(*strat);
  Teko::addToStratimikosBuilder(strat);

  return strat;
}
}

int main(int argc, char* argv[]) {
  // calls MPI_Init and MPI_Finalize
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  using ST = double;
  using LO = Teko::LO;
  using GO = Teko::GO;
  using NT = Teko::NT;

  using crs_matrix_type = Tpetra::CrsMatrix<ST, LO, GO, NT>;
  using vector_type = Tpetra::Vector<ST, LO, GO, NT>;
  using map_type = Tpetra::Map<LO, GO, NT>;

  // build global communicator
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

  // read in parameter list
  Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::getParametersFromXmlFile("strat_example.xml");

  // Read in the matrix, store pointer as an RCP
  RCP<crs_matrix_type> A =
      Tpetra::MatrixMarket::Reader<crs_matrix_type>::readSparseFile("data/nsjac.mm", Comm);

  // read in the RHS vector
  RCP<const map_type> rangeMap = A->getRangeMap();
  RCP<vector_type> b =
      Tpetra::MatrixMarket::Reader<crs_matrix_type>::readVectorFile(
          "data/nsrhs_test.mm", Comm, rangeMap, false, false);

  // allocate vectors
  RCP<vector_type> x = rcp(new vector_type(A->getDomainMap()));
  x->putScalar(0.0);

  // Build Thyra linear algebra objects
  RCP<const Thyra::VectorSpaceBase<ST> > thyraRangeSpace =
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A->getRangeMap());

  RCP<const Thyra::VectorSpaceBase<ST> > thyraDomainSpace =
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A->getDomainMap());

  RCP<const Thyra::LinearOpBase<ST> > th_A =
      Thyra::tpetraLinearOp<ST, LO, GO, NT>(thyraRangeSpace, thyraDomainSpace, A);

  RCP<const Thyra::VectorBase<ST> > th_b =
      Thyra::createConstVector<ST, LO, GO, NT>(b, thyraRangeSpace);

  RCP<Thyra::VectorBase<ST> > th_x =
      Thyra::createVector<ST, LO, GO, NT>(x, thyraDomainSpace);

  // Build stratimikos solver
  /////////////////////////////////////////////////////////

  auto linearSolverBuilder = create_linear_solver_builder();

  Teko::addTekoToStratimikosBuilder(*linearSolverBuilder);
  linearSolverBuilder->setParameterList(paramList);

  RCP<Thyra::LinearOpWithSolveFactoryBase<ST> > lowsFactory =
      Thyra::createLinearSolveStrategy(*linearSolverBuilder);

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<ST> > th_invA =
      Thyra::linearOpWithSolve(*lowsFactory, th_A);

  Thyra::assign(th_x.ptr(), 0.0);
  Thyra::SolveStatus<ST> status =
      Thyra::solve<ST>(*th_invA, Thyra::NOTRANS, *th_b, th_x.ptr());

  if (Comm->getRank() == 0) {
    std::cout << "Solve status: " << status.message << std::endl;
  }

  return 0;
}