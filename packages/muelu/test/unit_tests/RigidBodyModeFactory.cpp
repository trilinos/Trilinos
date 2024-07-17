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

#include <MueLu_RigidBodyModeFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(RigidBodyModeFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  auto A         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(100);
  const auto map = A->getRangeMap();

  using RigidBodyModeFactory = MueLu::RigidBodyModeFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  {
    RigidBodyModeFactory rigidBodyModeFactory("TestName");
    rigidBodyModeFactory.setNumPDEs(2);

    MueLu::Level level;
    level.SetLevelID(0);
    auto nullSpace = MultiVectorFactory::Build(map, 100);
    level.Set("TestName", nullSpace);
    level.Keep("Nullspace", &rigidBodyModeFactory);
    rigidBodyModeFactory.DeclareInput(level);

    TEST_EQUALITY(false, level.IsAvailable("Nullspace", &rigidBodyModeFactory));
    rigidBodyModeFactory.Build(level);
    TEST_EQUALITY(true, level.IsAvailable("Nullspace", &rigidBodyModeFactory));
  }

  {
    RigidBodyModeFactory rigidBodyModeFactory(1);

    MueLu::Level level;
    level.SetLevelID(0);
    level.Set("A", A);
    level.Keep("Nullspace", &rigidBodyModeFactory);
    auto coords = MultiVectorFactory::Build(map, 100);
    level.Set("Coordinates", coords);

    TEST_EQUALITY(false, level.IsAvailable("Nullspace", &rigidBodyModeFactory));
    rigidBodyModeFactory.Build(level);
    TEST_EQUALITY(true, level.IsAvailable("Nullspace", &rigidBodyModeFactory));
  }
}

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(RigidBodyModeFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
