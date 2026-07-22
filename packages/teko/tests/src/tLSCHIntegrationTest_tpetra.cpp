// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tLSCHIntegrationTest_tpetra.hpp"

#include <string>

// Teko-Package includes
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_TpetraBlockPreconditioner.hpp"

#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

// Belos includes
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Tpetra includes
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_TpetraBlockPreconditioner.hpp"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;

void tLSCHIntegrationTest_tpetra::initializeTest() { tolerance_ = 1.0e-7; }

int tLSCHIntegrationTest_tpetra::runTest(int verbosity, std::ostream& stdstrm,
                                         std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tLSCHIntegrationTest_tpetra";

  status = test_hScaling(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"hScaling\" ... PASSED", "   \"hScaling\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tLSCHIntegrationTest_tpetra...PASSED",
                         "tLSCHIntegrationTest_tpetra...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tLSCHIntegrationTest_tpetra...FAILED");
  }

  return failcount;
}

bool tLSCHIntegrationTest_tpetra::test_hScaling(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  RCP<const Teuchos::Comm<int> > comm = GetComm_tpetra();

  // build some operators
  Teko::LinearOp F = Teko::Test::build2x2(comm, 1, 2, 2, 1);
  Teko::LinearOp G = Teko::Test::build2x2(comm, 1, -1, -3, 1);
  Teko::LinearOp D = Teko::Test::build2x2(comm, 1, -3, -1, 1);

  ST diag[2];

  diag[0]          = 1.0 / 3.0;
  diag[1]          = 1.0 / 2.0;
  Teko::LinearOp M = Teko::Test::DiagMatrix_tpetra(2, diag, "M");

  diag[0]          = 5.0;
  diag[1]          = 9.0;
  Teko::LinearOp H = Teko::Test::DiagMatrix_tpetra(2, diag, "H");

  Teko::LinearOp A = Thyra::block2x2<ST>(F, G, D, Teuchos::null);

  Teko::LinearOp exact;
  {
    // build some operators
    Teko::LinearOp D0 = Teko::Test::build2x2(comm, -1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0);
    Teko::LinearOp D1 = Teko::Test::build2x2(comm, -1.5, -3.0, -3.0, -5.5);
    Teko::LinearOp U  = Teko::Test::build2x2(comm, -0.5, -1.5, -0.5, -0.5);

    exact = Thyra::block2x2<ST>(D0, U, Teuchos::null, D1);
  }

  RCP<Teko::InverseLibrary> invLib       = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> invFact      = invLib->getInverseFactory("Ifpack2");
  RCP<Teko::NS::InvLSCStrategy> strategy = rcp(new Teko::NS::InvLSCStrategy(invFact, M));
  strategy->setHScaling(Teko::getDiagonal(H));
  strategy->setUseFullLDU(false);

  RCP<Teko::BlockPreconditionerFactory> precFact =
      rcp(new Teko::NS::LSCPreconditionerFactory(strategy));
  RCP<Teko::BlockPreconditionerState> bps =
      Teuchos::rcp_dynamic_cast<Teko::BlockPreconditionerState>(
          precFact->buildPreconditionerState());
  Teko::LinearOp prec = precFact->buildPreconditionerOperator(A, *bps);

  Teko::BlockedLinearOp bA = Teko::toBlockedLinearOp(A);
  std::stringstream ss;
  ss << "invF = " << Teuchos::describe(*strategy->getInvF(bA, *bps), Teuchos::VERB_EXTREME)
     << std::endl;
  ss << "invBQBt = " << Teuchos::describe(*strategy->getInvBQBt(bA, *bps), Teuchos::VERB_EXTREME)
     << std::endl;
  ss << "invF = " << Teuchos::describe(*strategy->getInvBHBt(bA, *bps), Teuchos::VERB_EXTREME)
     << std::endl;
  ss << "invMass = " << Teuchos::describe(*strategy->getInvMass(bA, *bps), Teuchos::VERB_EXTREME)
     << std::endl;
  ss << "HScaling = " << Teuchos::describe(*strategy->getHScaling(bA, *bps), Teuchos::VERB_EXTREME)
     << std::endl;
  ss << "prec = " << Teuchos::describe(*prec, Teuchos::VERB_EXTREME) << std::endl;

  // construct a couple of vectors
  const RCP<Tpetra::Map<LO, GO, NT> > map =
      rcp(new Tpetra::Map<LO, GO, NT>(2, 0, GetComm_tpetra()));
  Tpetra::Vector<ST, LO, GO, NT> ea(map), eb(map);
  const RCP<const Thyra::MultiVectorBase<ST> > x = BlockVector(ea, eb, prec->domain());
  const RCP<Thyra::MultiVectorBase<ST> > y       = Thyra::createMembers(prec->range(), 1);
  ea.replaceGlobalValue(0, 1.0);
  ea.replaceGlobalValue(1, 0.0);
  eb.replaceGlobalValue(0, 0.0);
  eb.replaceGlobalValue(1, 0.0);
  Thyra::apply(*prec, Thyra::NOTRANS, *x, y.ptr());
  ss << "prec = " << Teuchos::describe(*y, Teuchos::VERB_EXTREME) << std::endl;
  ea.replaceGlobalValue(0, 0.0);
  ea.replaceGlobalValue(1, 1.0);
  eb.replaceGlobalValue(0, 0.0);
  eb.replaceGlobalValue(1, 0.0);
  Thyra::apply(*prec, Thyra::NOTRANS, *x, y.ptr());
  ss << "prec = " << Teuchos::describe(*y, Teuchos::VERB_EXTREME) << std::endl;
  ea.replaceGlobalValue(0, 0.0);
  ea.replaceGlobalValue(1, 0.0);
  eb.replaceGlobalValue(0, 1.0);
  eb.replaceGlobalValue(1, 0.0);
  Thyra::apply(*prec, Thyra::NOTRANS, *x, y.ptr());
  ss << "prec = " << Teuchos::describe(*y, Teuchos::VERB_EXTREME) << std::endl;
  ea.replaceGlobalValue(0, 0.0);
  ea.replaceGlobalValue(1, 0.0);
  eb.replaceGlobalValue(0, 0.0);
  eb.replaceGlobalValue(1, 1.0);
  Thyra::apply(*prec, Thyra::NOTRANS, *x, y.ptr());
  ss << "prec = " << Teuchos::describe(*y, Teuchos::VERB_EXTREME) << std::endl;

  Thyra::LinearOpTester<ST> tester;
  tester.show_all_tests(true);
  Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss), "|||");
  const bool result = tester.compare(*prec, *exact, Teuchos::ptrFromRef(fos));
  TEST_ASSERT(result, std::endl
                          << "   tLSCHIntegration::test_hScaling "
                          << ": Comparing preconditioner to exactly computed version");

  const bool result2 =
      tester.compare(*H, *strategy->getHScaling(bA, *bps), Teuchos::ptrFromRef(fos));
  TEST_ASSERT(result2, std::endl
                           << "   tLSCHIntegration::test_hScaling "
                           << ": Comparing scaling of H operator");

  const bool result3 = tester.compare(*build2x2(comm, 0.208333333333333, 0.375, 0.375, 0.875),
                                      *strategy->getInvBQBt(bA, *bps), Teuchos::ptrFromRef(fos));
  TEST_ASSERT(result3, std::endl
                           << "   tLSCHIntegration::test_hScaling "
                           << ": Comparing inv(BQBt) operator");

  const bool result4 = tester.compare(
      *build2x2(comm, 0.077777777777778, 0.177777777777778, 0.177777777777778, 0.477777777777778),
      *strategy->getInvBHBt(bA, *bps), Teuchos::ptrFromRef(fos));
  TEST_ASSERT(result4, std::endl
                           << "   tLSCHIntegration::test_hScaling "
                           << ": Comparing inv(BHBt) operator");

  if (not allPassed || verbosity >= 10) os << ss.str();

  return allPassed;
}

}  // namespace Test
}  // end namespace Teko
