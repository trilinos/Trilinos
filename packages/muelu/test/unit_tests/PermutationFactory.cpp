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

#include <MueLu_PermutationFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(PermutationFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  auto A         = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(100);
  const auto map = A->getRangeMap();

  MueLu::PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> permFact;
  permFact.SetParameter("PermutationStrategy", Teuchos::ParameterEntry(std::string("Local")));
  permFact.SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(std::string("")));
  permFact.SetFactory("PermutationRowMapFactory", Teuchos::null);

  MueLu::Level level;
  level.SetLevelID(0);
  level.Set("A", A);

  level.Request("A", &permFact);
  permFact.Build(level);

  auto coords = MultiVectorFactory::Build(map, 100);
}

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(PermutationFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
