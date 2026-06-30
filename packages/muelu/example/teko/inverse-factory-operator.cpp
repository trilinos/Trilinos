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
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StridedTpetraOperator.hpp"
#include "Teko_TpetraInverseFactoryOperator.hpp"
#include "Teko_ConfigDefs.hpp"

// Belos includes
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"

// Stratimikos includes
#include "Stratimikos_MueLuHelpers.hpp"

#include <iostream>
#include <vector> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

namespace {
Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> create_linear_solver_builder() {
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
}  // namespace

int main(int argc, char* argv[]) {
  // calls MPI_Init and MPI_Finalize
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  using ST = double;
  using LO = Teko::LO;
  using GO = Teko::GO;
  using NT = Teko::NT;

  using crs_matrix_type = Tpetra::CrsMatrix<ST, LO, GO, NT>;
  using vector_type     = Tpetra::Vector<ST, LO, GO, NT>;
  using mv_type         = Tpetra::MultiVector<ST, LO, GO, NT>;
  using operator_type   = Tpetra::Operator<ST, LO, GO, NT>;
  using map_type        = Tpetra::Map<LO, GO, NT>;

  // read in parameter list
  Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::getParametersFromXmlFile("simple_example.xml");

  // build global communicator
  RCP<const Teuchos::Comm<int> > Comm = Tpetra::getDefaultComm();

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

  // Break apart the strided linear system
  /////////////////////////////////////////////////////////

  // Block the linear system using a strided tpetra operator
  std::vector<int> vec(2);
  vec[0] = 2;
  vec[1] = 1; /*@ \label{lned:define-strided} @*/

  RCP<Teko::TpetraHelpers::StridedTpetraOperator> stridedAConcrete =
      Teuchos::rcp(new Teko::TpetraHelpers::StridedTpetraOperator(vec, A));

  // Build the preconditioner /*@ \label{lned:construct-prec} @*/
  /////////////////////////////////////////////////////////

  // build an InverseLibrary and inverse factory
  auto strat                        = create_linear_solver_builder();
  RCP<Teko::InverseLibrary> invLib  = Teko::InverseLibrary::buildFromParameterList(*paramList, strat);
  RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory("SIMPLE");

  // Create the initial preconditioner, and build it from strided_A
  RCP<Teko::TpetraHelpers::InverseFactoryOperator> prec_A =
      rcp(new Teko::TpetraHelpers::InverseFactoryOperator(inverse));
  prec_A->initInverse();

  RCP<operator_type> strided_A =
      Teuchos::rcp_dynamic_cast<operator_type>(stridedAConcrete);
  prec_A->rebuildInverseOperator(strided_A);

  // Build and solve the linear system
  /////////////////////////////////////////////////////////

  // Setup the linear solve: notice A is used directly
  RCP<Belos::LinearProblem<ST, mv_type, operator_type> > problem =
      rcp(new Belos::LinearProblem<ST, mv_type, operator_type>(A, x, b)); /*@ \label{lned:aztec-solve} @*/

  // build the solver
  RCP<operator_type> precOp = Teuchos::rcp_dynamic_cast<operator_type>(prec_A);
  problem->setRightPrec(precOp);

  const bool set = problem->setProblem();
  if (!set) {
    if (Comm->getRank() == 0) {
      std::cerr << "Belos::LinearProblem failed to set up." << std::endl;
    }
    return -1;
  }

  Belos::SolverFactory<ST, mv_type, operator_type> factory;
  RCP<Teuchos::ParameterList> belosList = rcp(new Teuchos::ParameterList());
  belosList->set("Maximum Iterations", 1000);
  belosList->set("Convergence Tolerance", 1e-5);
  belosList->set("Num Blocks", 1000);
  belosList->set("Verbosity",
                 Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::FinalSummary);

  RCP<Belos::SolverManager<ST, mv_type, operator_type> > solver =
      factory.create("GMRES", belosList);

  solver->setProblem(problem);

  // solve the linear system
  Belos::ReturnType result = solver->solve();

  if (Comm->getRank() == 0) {
    std::cout << "Belos GMRES result: "
              << (result == Belos::Converged ? "CONVERGED" : "NOT CONVERGED")
              << std::endl;
  }

  return (result == Belos::Converged) ? 0 : 1;
}
