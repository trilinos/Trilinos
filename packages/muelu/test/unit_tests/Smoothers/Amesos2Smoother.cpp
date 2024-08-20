// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Amesos2_config.h>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Amesos2Smoother.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

namespace MueLuTests {

// this namespace already has:  #include "MueLu_UseShortNames.hpp"
using namespace TestHelpers::Smoothers;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Amesos2Smoother, NotSetup, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
#if defined HAVE_AMESOS2_KLU2 || defined HAVE_AMESOS2_SUPERLU
    testApplyNoSetup(Amesos2Smoother(), out, success);
#endif
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Amesos2Smoother, Apply_Correctness, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra) {
    Teuchos::RCP<Amesos2Smoother> smoother;
#ifdef HAVE_AMESOS2_KLU2
    smoother = Teuchos::rcp(new Amesos2Smoother("Klu"));
    testDirectSolver(*smoother, out, success);
#endif

#ifdef HAVE_AMESOS2_SUPERLU
    smoother = Teuchos::rcp(new Amesos2Smoother("Superlu"));
    testDirectSolver(*smoother, out, success);
#endif
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Amesos2Smoother, NotSetup, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Amesos2Smoother, Apply_Correctness, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
