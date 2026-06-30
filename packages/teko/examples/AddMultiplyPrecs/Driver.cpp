// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommHelpers.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

// Tpetra includes
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"
#include "Teko_GaussSeidelPreconditionerFactory.hpp"
#include "Teko_BlockInvDiagonalStrategy.hpp"
#include "Teko_StridedTpetraOperator.hpp"
#include "Teko_TpetraBlockPreconditioner.hpp"
#include "Teko_AddPreconditionerFactory.hpp"
#include "Teko_MultPreconditionerFactory.hpp"
#include "Teko_ConfigDefs.hpp"

// Belos includes
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"

#include <iostream>
#include <vector>

using Teuchos::FancyOStream;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;

void run_driver() {
  RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

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

  RCP<FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(0);

  const int numProc = comm->getSize();
  const int myPID   = comm->getRank();

  *out << "Approaching Barrier: proc = " << numProc << ", pid = " << myPID << std::endl;
  Teuchos::barrier(*comm);

  RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("solverparams.xml");

  *fos << "Reading matrix market files" << std::endl;

  RCP<crs_t> A = Tpetra::MatrixMarket::Reader<crs_t>::readSparseFile("./modified.mm", comm);

  RCP<const map_t> rangeMap  = A->getRangeMap();
  RCP<const map_t> domainMap = A->getDomainMap();

  RCP<vec_t> b = Tpetra::MatrixMarket::Reader<crs_t>::readVectorFile("./rhs_test.mm", comm,
                                                                     rangeMap, false, false);

  RCP<vec_t> x = Tpetra::MatrixMarket::Reader<crs_t>::readVectorFile("./lhs_test.mm", comm,
                                                                     domainMap, false, false);

  *fos << "Building strided operator" << std::endl;

  std::vector<int> vars(2);
  vars[0] = 2;
  vars[1] = 1;

  Teuchos::RCP<Teko::TpetraHelpers::StridedTpetraOperator> sA =
      Teuchos::rcp(new Teko::TpetraHelpers::StridedTpetraOperator(vars, A));

  RCP<Teko::InverseFactory> inverse = Teko::invFactoryFromParamList(*paramList, "Amesos2");

  RCP<Teko::BlockInvDiagonalStrategy> strategy = rcp(new Teko::InvFactoryDiagStrategy(inverse));

  RCP<Teko::BlockPreconditionerFactory> GSFactory =
      rcp(new Teko::GaussSeidelPreconditionerFactory(Teko::GS_UseLowerTriangle, strategy));

  RCP<Teko::BlockPreconditionerFactory> JacobiFactory =
      rcp(new Teko::JacobiPreconditionerFactory(strategy));

#ifdef ADD_PREC
  RCP<Teko::BlockPreconditionerFactory> MasterFactory =
      rcp(new Teko::AddPreconditionerFactory(GSFactory, JacobiFactory));
#else
  RCP<Teko::BlockPreconditionerFactory> MasterFactory =
      rcp(new Teko::MultPreconditionerFactory(GSFactory, JacobiFactory));
#endif

  Teko::TpetraHelpers::TpetraBlockPreconditioner MyPreconditioner(MasterFactory);
  MyPreconditioner.buildPreconditioner(sA);

  RCP<mv_t> X = x;
  RCP<mv_t> B = rcp(new mv_t(b->getMap(), 1));

  using problem_t        = Belos::LinearProblem<ST, mv_t, op_t>;
  RCP<problem_t> problem = rcp(new problem_t(A, X, B));

  RCP<op_t> precOp = Teuchos::rcp(&MyPreconditioner, false);
  problem->setRightPrec(precOp);

  const bool set = problem->setProblem();
  TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
                             "Belos::LinearProblem::setProblem() failed.");

  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations", 100);
  belosList.set("Convergence Tolerance", 1e-5);
  belosList.set("Num Blocks", 50);
  belosList.set("Verbosity",
                Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::FinalSummary);
  belosList.set("Output Frequency", 10);

  using solver_t = Belos::BlockGmresSolMgr<ST, mv_t, op_t>;
  solver_t solver(problem, rcpFromRef(belosList));

  *fos << "Starting Belos solve" << std::endl;

  Belos::ReturnType result = solver.solve();

  TEUCHOS_TEST_FOR_EXCEPTION(result != Belos::Converged, std::runtime_error,
                             "Belos solver failed to converge.");

  *fos << "Solve converged" << std::endl;
}

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  { run_driver(); }
  return 0;
}
