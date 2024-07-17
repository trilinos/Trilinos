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

#include <MueLu_AlgebraicPermutationStrategy.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MueLu_AlgebraicPermutationStrategy, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using TestFactory = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  auto A = TestFactory::Build1DPoisson(200);
  A->describe(out, Teuchos::VERB_EXTREME);
  auto map = A->getRowMap();

  Level currentLevel;
  MueLu::AlgebraicPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node> algebraicPermutationStrategy;
  algebraicPermutationStrategy.BuildPermutation(A, map, currentLevel, MueLu::NoFactory::get());

  TEST_EQUALITY(true, currentLevel.IsAvailable("A", MueLu::NoFactory::get()));
  TEST_EQUALITY(true, currentLevel.IsAvailable("permA", MueLu::NoFactory::get()));
  TEST_EQUALITY(true, currentLevel.IsAvailable("permP", MueLu::NoFactory::get()));
  TEST_EQUALITY(true, currentLevel.IsAvailable("permQT", MueLu::NoFactory::get()));
  TEST_EQUALITY(true, currentLevel.IsAvailable("permScaling", MueLu::NoFactory::get()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(MueLu_AlgebraicPermutationStrategy, Utils, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  std::vector<LocalOrdinal> vec;

  std::vector<Scalar> values{10, 12, 40, 50, 32, 78};
  MueLu::sortingPermutation(values, vec);

  TEST_EQUALITY(values.size(), vec.size());

  TEST_EQUALITY(5, vec[0]);  // 78
  TEST_EQUALITY(3, vec[1]);  // 50
  TEST_EQUALITY(2, vec[2]);  // 40
  TEST_EQUALITY(4, vec[3]);  // 32
  TEST_EQUALITY(1, vec[4]);  // 12
  TEST_EQUALITY(0, vec[5]);  // 10
}

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node)                                                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MueLu_AlgebraicPermutationStrategy, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(MueLu_AlgebraicPermutationStrategy, Utils, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
