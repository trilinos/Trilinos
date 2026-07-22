// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tPCDStrategy.hpp"

#include <string>

// Epetra includes
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

// Test-rig
#include "Test_Utils.hpp"

// Teko includes
#include "Teko_PCDStrategy.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StaticRequestCallback.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

RCP<Thyra::LinearOpBase<double> > buildExampleOp(int type, const Epetra_Comm& comm, double scale) {
  Epetra_Map map(3, 0, comm);
  RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy, map, 3));

  double values[3];
  int indices[3] = {0, 1, 2};

  switch (type) {
    case 3:
      values[0] = 9.0;
      values[1] = 8.0;
      values[2] = 7.0;
      mat->InsertGlobalValues(0, 3, values, indices);

      values[0] = 1.0;
      values[1] = 2.0;
      values[2] = 3.0;
      mat->InsertGlobalValues(1, 3, values, indices);

      values[0] = 4.0;
      values[1] = 5.0;
      values[2] = 6.0;
      mat->InsertGlobalValues(2, 3, values, indices);
      break;
    case 2:
      values[0] = 5.0;
      values[1] = 0.0;
      values[2] = 6.0;
      mat->InsertGlobalValues(0, 3, values, indices);

      values[0] = 2.0;
      values[1] = 6.0;
      values[2] = 0.0;
      mat->InsertGlobalValues(1, 3, values, indices);

      values[0] = 0.0;
      values[1] = 3.0;
      values[2] = 7.0;
      mat->InsertGlobalValues(2, 3, values, indices);
      break;
    case 1:
      values[0] = 1.0;
      values[1] = 0.0;
      values[2] = 2.0;
      mat->InsertGlobalValues(0, 3, values, indices);

      values[0] = 2.0;
      values[1] = 5.0;
      values[2] = 0.0;
      mat->InsertGlobalValues(1, 3, values, indices);

      values[0] = 0.0;
      values[1] = 3.0;
      values[2] = 1.0;
      mat->InsertGlobalValues(2, 3, values, indices);
    default: break;
  }

  mat->FillComplete();
  mat->Scale(scale);

  return Thyra::nonconstEpetraLinearOp(mat);
}

void tPCDStrategy::initializeTest() { tolerance_ = 1.0e-15; }

int tPCDStrategy::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                          int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tPCDStrategy";

  status = test_PCDStrategy(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"PCDStrategy\" ... PASSED", "   \"PCDStrategy\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tPCDStrategy...PASSED", "tPCDStrategy...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tPCDStrategy...FAILED");
  }

  return failcount;
}

bool tPCDStrategy::test_PCDStrategy(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  Thyra::LinearOpTester<double> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(tolerance_);

  Teko::LinearOp M1 = buildExampleOp(1, *GetComm(), 1.0);
  Teko::LinearOp M2 = buildExampleOp(2, *GetComm(), 1.0);
  Teko::LinearOp A  = Thyra::nonconstBlock2x2<double>(
      buildExampleOp(2, *GetComm(), 1.0), buildExampleOp(1, *GetComm(), 2.0),
      buildExampleOp(1, *GetComm(), 3.0), buildExampleOp(2, *GetComm(), 4.0));
  Teko::BlockedLinearOp blkA = Teko::toBlockedLinearOp(A);

  Teko::LinearOp laplace = Teko::explicitMultiply(M2, M1);
  Teko::LinearOp Fp      = M1;
  Teko::LinearOp Qp      = buildExampleOp(3, *GetComm(), 1.0);

  std::string pcdStr      = Teko::NS::PCDStrategy::getPCDString();
  std::string presLapStr  = Teko::NS::PCDStrategy::getPressureLaplaceString();
  std::string presMassStr = Teko::NS::PCDStrategy::getPressureMassString();

  RCP<Teko::StaticRequestCallback<Teko::LinearOp> > pcdCb =
      rcp(new Teko::StaticRequestCallback<Teko::LinearOp>(pcdStr, Fp));
  RCP<Teko::StaticRequestCallback<Teko::LinearOp> > lapCb =
      rcp(new Teko::StaticRequestCallback<Teko::LinearOp>(presLapStr, laplace));
  RCP<Teko::StaticRequestCallback<Teko::LinearOp> > massCb =
      rcp(new Teko::StaticRequestCallback<Teko::LinearOp>(presMassStr, Qp));

  RCP<Teko::RequestHandler> rh = rcp(new Teko::RequestHandler);
  rh->addRequestCallback(pcdCb);
  rh->addRequestCallback(lapCb);
  rh->addRequestCallback(massCb);

  RCP<Teko::InverseLibrary> invLib      = Teko::InverseLibrary::buildFromStratimikos();
  RCP<Teko::InverseFactory> inverseFact = invLib->getInverseFactory("Amesos");

  // test the inverse of the 00 block
  {
    bool result;
    std::stringstream ss;
    Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss), "      |||");

    Teko::LinearOp invA00 = inverseFact->buildInverse(M2);
    Teko::BlockPreconditionerState state;
    Teko::NS::PCDStrategy strategy(inverseFact, inverseFact);
    strategy.setRequestHandler(rh);

    // build the state object
    // state.addLinearOp(presLapStr,  laplace);
    // state.addLinearOp(pcdStr,      Fp);
    // state.addLinearOp(presMassStr, Qp);

    // test the 0,0 block inverses: hat
    Teko::LinearOp hatInvA00 = strategy.getHatInvA00(blkA, state);
    ss.str("");
    result = tester.compare(*hatInvA00, *invA00, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tPCDStrategy::test_PCDStrategy " << toString(status)
                            << " : Comparing invHatA00 operators");
    if (not result || verbosity >= 10) os << ss.str();

    // test the 0,0 block inverses: tilde
    Teko::LinearOp tildeInvA00 = strategy.getTildeInvA00(blkA, state);
    ss.str("");
    result = tester.compare(*tildeInvA00, *invA00, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tPCDStrategy::test_PCDStrategy " << toString(status)
                            << " : Comparing invTildeA00 operators");
    if (not result || verbosity >= 10) os << ss.str();
  }

  // test the inverse of the Schur complement
  {
    bool result;
    std::stringstream ss;
    Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss), "      |||");

    Teko::LinearOp invLaplace = inverseFact->buildInverse(laplace);
    Teko::LinearOp iQp        = Teko::getInvDiagonalOp(Qp);
    Teko::LinearOp invSchur   = multiply(iQp, Fp, invLaplace);

    Teko::BlockPreconditionerState state;
    Teko::NS::PCDStrategy strategy(inverseFact, inverseFact);
    strategy.setRequestHandler(rh);

    // build the state object
    // state.addLinearOp(presLapStr,  laplace);
    // state.addLinearOp(pcdStr,      Fp);
    // state.addLinearOp(presMassStr, Qp);

    // test the 0,0 block inverses: hat
    Teko::LinearOp invS_strategy = strategy.getInvS(blkA, state);
    ss.str("");
    result = tester.compare(*invS_strategy, *invSchur, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tPCDStrategy::test_PCDStrategy " << toString(status)
                            << " : Comparing invS operators");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

}  // namespace Test
}  // end namespace Teko
