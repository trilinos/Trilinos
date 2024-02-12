// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
