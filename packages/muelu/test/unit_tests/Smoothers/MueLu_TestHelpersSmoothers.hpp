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
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_BlockedMultiVector.hpp>

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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void testDirectSolverBlocked(MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>& smoother,
                             Teuchos::FancyOStream& out, bool& success) {
#include "MueLu_UseShortNames.hpp"

  using ST             = Teuchos::ScalarTraits<SC>;
  using magnitude_type = typename ST::magnitudeType;
  using MT             = Teuchos::ScalarTraits<magnitude_type>;

  const GO n0 = 50;
  const GO n1 = 75;
  const GO n  = n0 + n1;

  auto comm = Parameters::getDefaultComm();

  // Full map
  RCP<const Map> fullMap = MapFactory::Build(Parameters::getLib(), n, 0, comm);

  // Build disjoint submaps from full GIDs
  Teuchos::Array<GO> gids0, gids1;
  const Teuchos::ArrayView<const GO> myGids = fullMap->getLocalElementList();

  for (size_t k = 0; k < (size_t)myGids.size(); ++k) {
    if (myGids[k] < n0)
      gids0.push_back(myGids[k]);
    else
      gids1.push_back(myGids[k]);
  }

  RCP<const Map> map0 = MapFactory::Build(Parameters::getLib(),
                                          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                          gids0(), 0, comm);

  RCP<const Map> map1 = MapFactory::Build(Parameters::getLib(),
                                          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                          gids1(), 0, comm);

  std::vector<RCP<const Map>> maps(2);
  maps[0] = map0;
  maps[1] = map1;

  RCP<const Xpetra::MapExtractor<SC, LO, GO, NO>> rangeMapExtractor =
      Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullMap, maps, false);

  RCP<const Xpetra::MapExtractor<SC, LO, GO, NO>> domainMapExtractor =
      Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullMap, maps, false);

  // Build simple diagonal matrices on the submaps
  RCP<CrsMatrixWrap> A00, A11;
  {
    RCP<CrsMatrix> crs0 = CrsMatrixFactory::Build(map0, 1);
    for (size_t l = 0; l < (size_t)map0->getLocalNumElements(); ++l) {
      GO gid = map0->getGlobalElement(l);
      crs0->insertGlobalValues(gid, Teuchos::tuple<GO>(gid), Teuchos::tuple<SC>(ST::one()));
    }
    crs0->fillComplete(map0, map0);
    A00 = rcp(new CrsMatrixWrap(crs0));
  }
  {
    RCP<CrsMatrix> crs1 = CrsMatrixFactory::Build(map1, 1);
    for (size_t l = 0; l < (size_t)map1->getLocalNumElements(); ++l) {
      GO gid = map1->getGlobalElement(l);
      crs1->insertGlobalValues(gid, Teuchos::tuple<GO>(gid), Teuchos::tuple<SC>(ST::one()));
    }
    crs1->fillComplete(map1, map1);
    A11 = rcp(new CrsMatrixWrap(crs1));
  }

  // Build 2x2 blocked matrix with only diagonal blocks
  RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> bA =
      rcp(new Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>(rangeMapExtractor, domainMapExtractor, 2));

  bA->setMatrix(0, 0, A00);
  bA->setMatrix(1, 1, A11);
  bA->fillComplete();

  RCP<Matrix> A = bA;
  setupSmoother(A, smoother, out, success);

  // Exact blocked solution
  RCP<MultiVector> X0 = MultiVectorFactory::Build(map0, 1);
  RCP<MultiVector> X1 = MultiVectorFactory::Build(map1, 1);
  X0->setSeed(846930886);
  X1->setSeed(846930887);
  X0->randomize();
  X1->randomize();

  {
    Teuchos::Array<magnitude_type> norms(1);
    X0->norm2(norms);
    if (norms[0] != MT::zero()) X0->scale(ST::one() / norms[0]);
    X1->norm2(norms);
    if (norms[0] != MT::zero()) X1->scale(ST::one() / norms[0]);
  }

  std::vector<RCP<MultiVector>> xBlocks(2);
  xBlocks[0] = X0;
  xBlocks[1] = X1;

  RCP<Xpetra::BlockedMultiVector<SC, LO, GO, NO>> Xexact =
      rcp(new Xpetra::BlockedMultiVector<SC, LO, GO, NO>(domainMapExtractor->getBlockedMap(), xBlocks));

  // For identity block matrix, RHS = Xexact
  RCP<MultiVector> RHSmerged = Xexact->Merge();

  RCP<Xpetra::BlockedMultiVector<SC, LO, GO, NO>> RHS =
      rcp(new Xpetra::BlockedMultiVector<SC, LO, GO, NO>(rangeMapExtractor, RHSmerged));

  // Initial guess
  RCP<Xpetra::BlockedMultiVector<SC, LO, GO, NO>> X =
      rcp(new Xpetra::BlockedMultiVector<SC, LO, GO, NO>(domainMapExtractor->getBlockedMap(), 1, true));

  smoother.Apply(*X, *RHS);

  Teuchos::Array<magnitude_type> residualNorms = Utilities::ResidualNorm(*bA, *X->Merge(), *RHS->Merge());
  out << "||Blocked Residual|| = "
      << std::setiosflags(std::ios::fixed) << std::setprecision(20)
      << residualNorms[0] << std::endl;

  TEST_EQUALITY(residualNorms[0] < 100 * MT::eps(), true);
}

}  // namespace Smoothers
}  // namespace TestHelpers
}  // namespace MueLuTests

#endif  // MUELU_TEST_HELPERS_SMOOTHERS_H
