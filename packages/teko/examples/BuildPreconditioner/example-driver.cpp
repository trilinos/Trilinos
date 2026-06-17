// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teuchos includes
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

// Tpetra includes
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StridedTpetraOperator.hpp"
#include "Teko_TpetraBlockPreconditioner.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#include "Teko_ConfigDefs.hpp"

// Belos includes
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"

#include <iostream>
#include <vector>

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  using ST = double;
  using LO = Teko::LO;
  using GO = Teko::GO;
  using NT = Teko::NT;

  using map_t = Tpetra::Map<LO, GO, NT>;
  using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;
  using vec_t = Tpetra::Vector<ST, LO, GO, NT>;
  using mv_t  = Tpetra::MultiVector<ST, LO, GO, NT>;
  using op_t  = Tpetra::Operator<ST, LO, GO, NT>;

  auto comm = Tpetra::getDefaultComm();

  // Read in the matrix
  RCP<crs_t> A = Tpetra::MatrixMarket::Reader<crs_t>::readSparseFile("../data/nsjac.mm", comm);

  TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), std::runtime_error,
                             "Failed to read matrix from ../data/nsjac.mm");

  // Read in the RHS vector
  RCP<const map_t> rangeMap  = A->getRangeMap();
  RCP<const map_t> domainMap = A->getDomainMap();

  RCP<vec_t> b = Tpetra::MatrixMarket::Reader<crs_t>::readVectorFile("../data/nsrhs_test.mm", comm,
                                                                     rangeMap, false, false);

  TEUCHOS_TEST_FOR_EXCEPTION(b.is_null(), std::runtime_error,
                             "Failed to read RHS vector from ../data/nsrhs_test.mm");

  // Allocate solution vector
  RCP<vec_t> x = rcp(new vec_t(domainMap));
  x->putScalar(0.0);

  // Break apart the strided linear system
  /////////////////////////////////////////////////////////

  std::vector<int> vec(2);
  vec[0] = 2;
  vec[1] = 1;

  Teuchos::RCP<Teko::TpetraHelpers::StridedTpetraOperator> sA =
      Teuchos::rcp(new Teko::TpetraHelpers::StridedTpetraOperator(vec, A));

  // Build the preconditioner
  /////////////////////////////////////////////////////////

  RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

  RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory("Amesos2");

  RCP<Teko::NS::LSCStrategy> strategy = rcp(new Teko::NS::InvLSCStrategy(inverse, true));

  RCP<Teko::BlockPreconditionerFactory> precFact =
      rcp(new Teko::NS::LSCPreconditionerFactory(strategy));

  Teko::TpetraHelpers::TpetraBlockPreconditioner prec(precFact);
  prec.buildPreconditioner(sA);

  // Build and solve the linear system
  /////////////////////////////////////////////////////////

  RCP<mv_t> X = x;
  RCP<mv_t> B = b;

  using problem_t        = Belos::LinearProblem<ST, mv_t, op_t>;
  RCP<problem_t> problem = rcp(new problem_t(A, X, B));

  RCP<op_t> precOp = Teuchos::rcp(&prec, false);
  problem->setRightPrec(precOp);

  const bool set = problem->setProblem();
  TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
                             "Belos::LinearProblem::setProblem() failed.");

  using solver_t = Belos::BlockGmresSolMgr<ST, mv_t, op_t>;

  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations", 1000);
  belosList.set("Convergence Tolerance", 1e-5);
  belosList.set("Num Blocks", 1000);
  belosList.set("Verbosity",
                Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::FinalSummary);
  belosList.set("Output Frequency", 10);

  solver_t solver(problem, rcpFromRef(belosList));
  Belos::ReturnType result = solver.solve();

  TEUCHOS_TEST_FOR_EXCEPTION(result != Belos::Converged, std::runtime_error,
                             "Belos solver failed to converge.");

  return 0;
}
