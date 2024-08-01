// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_CoupledRBMFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MueLu_CoupledRBMFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  auto A         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(100);
  const auto map = A->getRangeMap();

  MueLu::CoupledRBMFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> coupledRBMFactory("TestName");
  coupledRBMFactory.setNumPDEs(20);
  coupledRBMFactory.setLastAcousticDOF(10);

  {
    MueLu::Level level;
    level.SetLevelID(0);
    auto nullSpace = MultiVectorFactory::Build(map, 100);
    level.Set("TestName", nullSpace);
    level.Keep("Nullspace", &coupledRBMFactory);
    coupledRBMFactory.DeclareInput(level);

    TEST_EQUALITY(false, level.IsAvailable("Nullspace", &coupledRBMFactory));
    coupledRBMFactory.Build(level);
    TEST_EQUALITY(true, level.IsAvailable("Nullspace", &coupledRBMFactory));
  }

  {
    MueLu::Level level;
    level.SetLevelID(0);
    level.Set("A", A);
    level.Keep("Nullspace", &coupledRBMFactory);
    auto coords = MultiVectorFactory::Build(map, 100);
    level.Set("Coordinates", coords);

    TEST_EQUALITY(false, level.IsAvailable("Nullspace", &coupledRBMFactory));
    coupledRBMFactory.Build(level);
    TEST_EQUALITY(true, level.IsAvailable("Nullspace", &coupledRBMFactory));
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MueLu_CoupledRBMFactory, BuildRBM, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  MueLu::CoupledRBMFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> coupledRBMFactory(20);
  coupledRBMFactory.setLastAcousticDOF(10);

  auto A         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(100);
  const auto map = A->getRangeMap();

  auto coords = MultiVectorFactory::Build(map, 100);
  decltype(coords) nullSpace;

  coupledRBMFactory.BuildRBM(A, coords, nullSpace);

  TEST_INEQUALITY(nullSpace, Teuchos::null);
}

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node)                                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MueLu_CoupledRBMFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MueLu_CoupledRBMFactory, BuildRBM, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
