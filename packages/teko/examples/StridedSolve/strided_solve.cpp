// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <sys/types.h>
#include <unistd.h>

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
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
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_StridedTpetraOperator.hpp"
#include "Teko_TpetraBlockPreconditioner.hpp"
#include "Teko_ConfigDefs.hpp"

// Belos includes
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

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

  std::string solveName = "Amesos2";
  if (argc > 1) solveName = argv[1];

  RCP<FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(0);

  *fos << "Using \"" << solveName << "\" for approximate solve" << std::endl;

  const int numProc = comm->getSize();
  const int myPID   = comm->getRank();

  std::cout << "MPI_PID = " << myPID << ", UNIX_PID = " << getpid() << std::endl;

  *out << "Approaching Barrier: proc = " << numProc << ", pid = " << myPID << std::endl;
  Teuchos::barrier(*comm);

  RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile("solverparams.xml");

  *fos << "Reading matrix market file" << std::endl;

  RCP<crs_t> A = Tpetra::MatrixMarket::Reader<crs_t>::readSparseFile("../data/nsjac.mm", comm);

  RCP<const map_t> rangeMap  = A->getRangeMap();
  RCP<const map_t> domainMap = A->getDomainMap();

  RCP<vec_t> b = Tpetra::MatrixMarket::Reader<crs_t>::readVectorFile("../data/nsrhs_test.mm", comm,
                                                                     rangeMap, false, false);

  RCP<vec_t> x = Tpetra::MatrixMarket::Reader<crs_t>::readVectorFile("../data/nslhs_test.mm", comm,
                                                                     domainMap, false, false);

  *fos << "Building strided operator" << std::endl;

  std::vector<int> vars(2);
  vars[0] = 2;
  vars[1] = 1;

  Teuchos::RCP<Teko::TpetraHelpers::StridedTpetraOperator> sA =
      Teuchos::rcp(new Teko::TpetraHelpers::StridedTpetraOperator(vars, A));

  double alpha                      = 0.9;
  RCP<Teko::InverseFactory> inverse = Teko::invFactoryFromParamList(*paramList, solveName);

  RCP<Teko::BlockPreconditionerFactory> precFact =
      rcp(new Teko::NS::SIMPLEPreconditionerFactory(inverse, alpha));

  *fos << "Preconditioner factory built" << std::endl;

  Teko::TpetraHelpers::TpetraBlockPreconditioner prec(precFact);
  prec.buildPreconditioner(sA);

  *fos << "Preconditioner built" << std::endl;

  RCP<mv_t> X = x;
  RCP<mv_t> B = b;

  using problem_t        = Belos::LinearProblem<ST, mv_t, op_t>;
  RCP<problem_t> problem = rcp(new problem_t(A, X, B));

  RCP<op_t> precOp = Teuchos::rcp(&prec, false);
  problem->setRightPrec(precOp);

  const bool set = problem->setProblem();
  TEUCHOS_TEST_FOR_EXCEPTION(!set, std::runtime_error,
                             "Belos::LinearProblem::setProblem() failed.");

  *fos << "Setting solver parameters" << std::endl;

  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations", 1000);
  belosList.set("Convergence Tolerance", 1e-5);
  belosList.set("Num Blocks", 50);
  belosList.set("Verbosity",
                Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::FinalSummary);
  belosList.set("Output Frequency", 10);

  using solver_t = Belos::BlockGmresSolMgr<ST, mv_t, op_t>;
  solver_t solver(problem, rcpFromRef(belosList));

  *fos << "Solving" << std::endl;

  Belos::ReturnType result = solver.solve();

  TEUCHOS_TEST_FOR_EXCEPTION(result != Belos::Converged, std::runtime_error,
                             "Belos solver failed to converge.");

  *fos << "Solve converged" << std::endl;

  return 0;
}
