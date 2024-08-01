// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <MueLu_TestHelpers.hpp>
#include "MueLu_TestHelpersSmoothers.hpp"

#include <MueLu_BelosSmoother.hpp>
#include <MueLu_Utilities.hpp>

namespace MueLuTests {

// this namespace already has:  #include "MueLu_UseShortNames.hpp"
using namespace TestHelpers::Smoothers;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BelosSmoother, NotSetup, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    BelosSmoother smoother("Block CG", Teuchos::ParameterList());
    testApplyNoSetup(smoother, out, success);
  }
}

// Tests interface to Belos's Gauss-Seidel preconditioner.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BelosSmoother, HardCodedResult_BlockCG, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::ParameterList paramList;
    paramList.set("Maximum Iterations", 10);
    BelosSmoother smoother("Block CG", paramList);

    typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms = testApply_A125_X1_RHS0(smoother, out, success);

    RCP<const Teuchos::Comm<int> > comm                                  = TestHelpers::Parameters::getDefaultComm();
    const typename Teuchos::ScalarTraits<SC>::magnitudeType expectedNorm = 1.2856486930664495771e-01;
    TEST_FLOATING_EQUALITY(residualNorms, expectedNorm, 1e4 * Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<SC>::magnitudeType>::eps());
  }

}  // Block CG

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BelosSmoother, NotSetup, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BelosSmoother, HardCodedResult_BlockCG, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
