// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Amesos_config.h>
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include "MueLu_AmesosSmoother.hpp"

namespace MueLuTests {

typedef double Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;
typedef Xpetra::EpetraNode Node;
#include "MueLu_UseShortNames.hpp"

using namespace TestHelpers::Smoothers;

TEUCHOS_UNIT_TEST(AmesosSmoother, NotSetup) {
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {
    testApplyNoSetup(AmesosSmoother(), out, success);
  }
}

TEUCHOS_UNIT_TEST(AmesosSmoother, Apply_Correctness) {
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra) {
    Teuchos::RCP<AmesosSmoother> smoother;
#ifdef HAVE_AMESOS_KLU
    smoother = rcp(new AmesosSmoother("Klu"));
    testDirectSolver(*smoother, out, success);
#endif

#ifdef HAVE_AMESOS_SUPERLU
    smoother = rcp(new AmesosSmoother("Superlu"));
    testDirectSolver(*smoother, out, success);
#endif
  }
}

}  // namespace MueLuTests

// TODO
//  Test with invalid std::string
//  Test with invalid parameterList? (== a characterization test for Amesos)
//  Test if paramList takes into account
//  Check if all the defaults that are used by MueLu are tested
