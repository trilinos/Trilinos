// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TEST_HELPERS_SMOOTHERS_H
#define MUELU_TEST_HELPERS_SMOOTHERS_H

#include <Teuchos_FancyOStream.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"

// Helper functions to test if derived classes conforms to the SmootherBase and SmootherPrototype interfaces

namespace MueLuTests {
namespace TestHelpers {
namespace Smoothers {

// SmootherPrototype test
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void testApplyNoSetup(const MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>& smoother, Teuchos::FancyOStream& out, bool& success) {
#include "MueLu_UseShortNames.hpp"

  GO numGlobalElements = 125;
  RCP<const Map> map   = MapFactory::Build(Parameters::getLib(), numGlobalElements, 0, Parameters::getDefaultComm());

  RCP<MultiVector> X   = MultiVectorFactory::Build(map, 1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map, 1);

  TEST_THROW(smoother.Apply(*X, *RHS), MueLu::Exceptions::RuntimeError);
}

// SmootherBase test
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
testApply(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
          const MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>& smoother,
          Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
          const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& RHS,
          Teuchos::FancyOStream& out, bool& success) {
#include "MueLu_UseShortNames.hpp"
  typedef Teuchos::ScalarTraits<SC> ST;

  Array<typename ST::magnitudeType> norms(1);

  RHS.norm2(norms);
  out << "||RHS|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  Teuchos::Array<typename ST::magnitudeType> initialNorms(1);
  X.norm2(initialNorms);
  out << "||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << initialNorms[0] << std::endl;

  smoother.Apply(X, RHS);  // TODO: bool const &InitialGuessIsZero=false

  Teuchos::Array<typename ST::magnitudeType> finalNorms(1);
  X.norm2(finalNorms);
  out << "||X_final|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(25) << norms[0] << std::endl;

  Teuchos::Array<typename ST::magnitudeType> residualNorms = Utilities::ResidualNorm(A, X, RHS);
  out << "||Residual|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << residualNorms[0] << std::endl;

  return residualNorms[0];
}

// SmootherBase test
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
testApply_X1_RHS0(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                  const MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>& smoother,
                  Teuchos::FancyOStream& out, bool& success) {
#include "MueLu_UseShortNames.hpp"

  RCP<MultiVector> X   = MultiVectorFactory::Build(A.getDomainMap(), 1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(A.getRangeMap(), 1);
  X->putScalar((SC)1.0);
  RHS->putScalar((SC)0.0);

  return testApply(A, smoother, *X, *RHS, out, success);
}

// SmootherBase test
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
testApply_X0_RandomRHS(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                       const MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>& smoother,
                       Teuchos::FancyOStream& out, bool& success) {
#include "MueLu_UseShortNames.hpp"

  RCP<MultiVector> X   = MultiVectorFactory::Build(A.getDomainMap(), 1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(A.getRangeMap(), 1);

  // Random X
  X->setSeed(846930886);
  X->randomize();

  // Normalize X
  typedef Teuchos::ScalarTraits<SC> ST;
  Array<typename ST::magnitudeType> norms(1);
  X->norm2(norms);
  X->scale(1 / norms[0]);

  // Compute RHS corresponding to X
  A.apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

  // Reset X to 0
  X->putScalar((SC)0.0);

  return testApply(A, smoother, *X, *RHS, out, success);
}

// SmootherPrototype helper function
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void setupSmoother(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                   MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>& smoother,
                   Teuchos::FancyOStream& out, bool& success) {
#include "MueLu_UseShortNames.hpp"

  Level level;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);
  smoother.Setup(level);
}

// SmootherPrototype test
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
testApply_A125_X1_RHS0(MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>& smoother,
                       Teuchos::FancyOStream& out, bool& success) {
#include "MueLu_UseShortNames.hpp"

  Teuchos::RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(125);

  setupSmoother(A, smoother, out, success);
  return testApply_X1_RHS0(*A, smoother, out, success);  // in MueLuTests::SmootherBase
}

// SmootherPrototype test
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
testApply_A125_X0_RandomRHS(MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>& smoother,
                            Teuchos::FancyOStream& out, bool& success) {
#include "MueLu_UseShortNames.hpp"

  Teuchos::RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(125);

  setupSmoother(A, smoother, out, success);
  return testApply_X0_RandomRHS(*A, smoother, out, success);  // in MueLuTests::SmootherBase
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void testDirectSolver(MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>& smoother,
                      Teuchos::FancyOStream& out, bool& success) {
#include "MueLu_UseShortNames.hpp"

  using ST             = Teuchos::ScalarTraits<SC>;
  using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  using MT             = Teuchos::ScalarTraits<magnitude_type>;

  magnitude_type residualNorms = testApply_A125_X0_RandomRHS(smoother, out, success);
  TEST_EQUALITY(residualNorms < 100 * MT::eps(), true);
}

}  // namespace Smoothers
}  // namespace TestHelpers
}  // namespace MueLuTests

#endif  // MUELU_TEST_HELPERS_SMOOTHERS_H
