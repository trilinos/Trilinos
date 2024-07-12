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

#include <MueLu_StratimikosSmoother.hpp>
#include <MueLu_Utilities.hpp>

namespace MueLuTests {

// this namespace already has:  #include "MueLu_UseShortNames.hpp"
using namespace TestHelpers::Smoothers;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StratimikosSmoother, NotSetup, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  if (!TYPE_EQUAL(Scalar, double)) {
    out << "Skipping for SC != double" << std::endl;
    return;
  }

  StratimikosSmoother smoother("STRATIMIKOS", Teuchos::ParameterList());
  testApplyNoSetup(smoother, out, success);
}

// Tests interface to Stratimikos's Gauss-Seidel preconditioner.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StratimikosSmoother, HardCodedResult_BlockCG, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  if (!TYPE_EQUAL(Scalar, double)) {
    out << "Skipping for SC != double" << std::endl;
    return;
  }

  Teuchos::ParameterList paramList;
  paramList.set("Linear Solver Type", "Belos");
  paramList.sublist("Linear Solver Types").sublist("Belos").set("Solver Type", "Block CG");
  paramList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block CG").set("Convergence Tolerance", 1e-2);
  paramList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block CG").set("Maximum Iterations", 100);
  paramList.set("Preconditioner Type", "None");
  StratimikosSmoother smoother("STRATIMIKOS", paramList);

  typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms = testApply_A125_X0_RandomRHS(smoother, out, success);

  RCP<const Teuchos::Comm<int> > comm                                  = TestHelpers::Parameters::getDefaultComm();
  const typename Teuchos::ScalarTraits<SC>::magnitudeType expectedNorm = 0.02118661536508828144;
  TEST_FLOATING_EQUALITY(residualNorms, expectedNorm, 0.1);

}  // Block CG

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StratimikosSmoother, NotSetup, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StratimikosSmoother, HardCodedResult_BlockCG, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
