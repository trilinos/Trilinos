// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA) && defined(HAVE_MUELU_THYRATPETRAADAPTERS) && defined(HAVE_MUELU_TEKO)

#include "MueLu_TestHelpers.hpp"
#include "MueLu_TekoMultiPhysBlockPreconditionerFactory.hpp"
#include "Stratimikos_MueLuHelpers.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_TpetraThyraConverter.hpp"

#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Tpetra_CrsMatrix.hpp"

#include "Xpetra_MultiVectorFactory.hpp"

namespace MueLuTests {

namespace {

using ST           = Teko::ST;
using LO           = Teko::LO;
using GO           = Teko::GO;
using NT           = Teko::NT;
using TpetraMap    = Tpetra::Map<LO, GO, NT>;
using TpetraMatrix = Tpetra::CrsMatrix<ST, LO, GO, NT>;
using TpetraMV     = Tpetra::MultiVector<ST, LO, GO, NT>;
using ThyraMV      = Thyra::MultiVectorBase<ST>;

// Block (i,j) of a 2x2 system: tridiagonal on the diagonal blocks, and a coupling block with entries
// at (k, k) and, if extraCoupling is set, at (k, k+1) on the off-diagonal blocks.
Teuchos::RCP<TpetraMatrix> buildBlock(const Teuchos::RCP<const TpetraMap>& rowMap,
                                      const Teuchos::RCP<const TpetraMap>& domainMap,
                                      bool diagonalBlock, ST scaling, bool extraCoupling) {
  const GO numGlobalCols = static_cast<GO>(domainMap->getGlobalNumElements());
  auto A                 = Teuchos::rcp(new TpetraMatrix(rowMap, 3));
  for (LO lid = 0; lid < static_cast<LO>(rowMap->getLocalNumElements()); ++lid) {
    const GO gid = rowMap->getGlobalElement(lid);
    auto insert  = [&](GO col, ST val) { A->insertGlobalValues(gid, Teuchos::tuple<GO>(col), Teuchos::tuple<ST>(val)); };
    if (diagonalBlock) {
      if (gid > 0) insert(gid - 1, -scaling);
      insert(gid, 2.5 * scaling);
      if (gid + 1 < numGlobalCols) insert(gid + 1, -scaling);
    } else {
      insert(gid, -0.1 * scaling);
      if (extraCoupling && gid + 1 < numGlobalCols) insert(gid + 1, -0.05 * scaling);
    }
  }
  A->fillComplete(domainMap, rowMap);
  return A;
}

Teko::BlockedLinearOp buildBlockedOperator(ST scaling = 1.0, bool extraCoupling = false) {
  auto comm              = MueLuTests::TestHelpers::Parameters::getDefaultComm();
  const GO numGlobalRows = 60;
  const auto map0        = Teuchos::rcp(new TpetraMap(numGlobalRows, 0, comm));
  const auto map1        = Teuchos::rcp(new TpetraMap(numGlobalRows, 0, comm));

  const std::vector<Teuchos::RCP<const TpetraMap> > maps = {map0, map1};

  auto blo = Teko::createBlockedOp();
  Teko::beginBlockFill(blo, 2, 2);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) {
      auto Aij = buildBlock(maps[i], maps[j], i == j, (i + 1) * scaling, extraCoupling);
      Teko::setBlock(i, j, blo, Thyra::createConstLinearOp<ST, LO, GO, NT>(Teuchos::rcp_implicit_cast<const Tpetra::Operator<ST, LO, GO, NT> >(Aij)));
    }
  Teko::endBlockFill(blo);
  return blo;
}

void initializeFactory(MueLu::TekoMultiPhysBlockPreconditionerFactory& factory) {
  Teuchos::ParameterList params;
  auto& amg = params.sublist("AMG Settings");
  amg.set("verbosity", "none");
  amg.set("smoother: type", "RELAXATION");
  amg.sublist("smoother: params").set("relaxation: type", "Jacobi");
  for (int block = 0; block < 2; ++block) {
    auto& sub = amg.sublist("subblockList" + std::to_string(block));
    sub.set("verbosity", "none");
    sub.set("coarse: max size", 20);
  }
  factory.initializeFromParameterList(params);
}

ST maxRelativeDifference(const ThyraMV& a, const ThyraMV& b) {
  const auto numVectors = a.domain()->dim();
  auto diff             = Thyra::createMembers(a.range(), numVectors);
  Thyra::assign(diff.ptr(), a);
  Thyra::update(-1.0, b, diff.ptr());
  Teuchos::Array<ST> diffNorms(numVectors), refNorms(numVectors);
  Thyra::norms_2(*diff, diffNorms());
  Thyra::norms_2(b, refNorms());
  ST maxDiff = 0.0;
  for (int k = 0; k < numVectors; ++k) maxDiff = std::max(maxDiff, diffNorms[k] / refNorms[k]);
  return maxDiff;
}

// Compare the monolithic matrix stored in the state against the blocked operator.
ST monolithicMismatch(const Teko::BlockedLinearOp& blo, const TpetraMatrix& monolithic) {
  auto x = Thyra::createMembers(blo->domain(), 2);
  Thyra::randomize(-1.0, 1.0, x.ptr());
  auto yRef = Thyra::createMembers(blo->range(), 2);
  Thyra::apply(*blo, Thyra::NOTRANS, *x, yRef.ptr());

  TpetraMV xFlat(monolithic.getDomainMap(), 2), yFlat(monolithic.getRangeMap(), 2);
  Teko::TpetraHelpers::blockThyraToTpetra(x, xFlat);
  monolithic.apply(xFlat, yFlat);
  auto y = Thyra::createMembers(blo->range(), 2);
  Teko::TpetraHelpers::blockTpetraToThyra(yFlat, y.ptr());
  return maxRelativeDifference(*y, *yRef);
}

}  // namespace

TEUCHOS_UNIT_TEST(TekoMultiPhysBlockPreconditionerFactory, Registration) {
  Stratimikos::registerMueLuTekoMultiPhysBlockAMG("Unit Test Block AMG");
  Teuchos::ParameterList settings;
  auto factory = Teko::PreconditionerFactory::buildPreconditionerFactory("Unit Test Block AMG", settings);
  TEST_ASSERT(!Teuchos::is_null(factory));
}

TEUCHOS_UNIT_TEST(TekoMultiPhysBlockPreconditionerFactory, RejectNonContiguousSubblocks) {
  MueLu::TekoMultiPhysBlockPreconditionerFactory factory;
  Teuchos::ParameterList params;
  auto& amg = params.sublist("AMG Settings");
  amg.sublist("subblockList0");
  amg.sublist("subblockList2");

  TEST_THROW(factory.initializeFromParameterList(params), MueLu::Exceptions::RuntimeError);
}

TEUCHOS_UNIT_TEST(TekoMultiPhysBlockPreconditionerFactory, AcceptTopLevelUserData) {
  using ST          = Teko::ST;
  using LO          = Teko::LO;
  using GO          = Teko::GO;
  using NT          = Teko::NT;
  using Map         = Xpetra::Map<LO, GO, NT>;
  using MultiVector = Xpetra::MultiVector<ST, LO, GO, NT>;

  MueLu::TekoMultiPhysBlockPreconditionerFactory factory;
  Teuchos::ParameterList params;
  auto& amg = params.sublist("AMG Settings");
  amg.sublist("subblockList0");
  amg.sublist("subblockList1");

  auto comm   = MueLuTests::TestHelpers::Parameters::getDefaultComm();
  auto map    = Xpetra::MapFactory<LO, GO, NT>::Build(Xpetra::UseTpetra, 4, 0, comm);
  auto coords = Xpetra::MultiVectorFactory<ST, LO, GO, NT>::Build(map, 1);
  params.sublist("user data").set<Teuchos::RCP<MultiVector> >("Coordinates0", coords);

  TEST_NOTHROW(factory.initializeFromParameterList(params));
}

TEUCHOS_UNIT_TEST(TekoMultiPhysBlockPreconditionerFactory, MonolithicConversionMatchesBlockedOperator) {
  using State = MueLu::TekoMultiPhysBlockPreconditionerState;

  MueLu::TekoMultiPhysBlockPreconditionerFactory factory;
  initializeFactory(factory);
  auto state     = factory.buildPreconditionerState();
  auto& blkState = Teuchos::dyn_cast<State>(*state);
  const ST tol   = 1e-12;

  // Initial conversion.
  auto blo = buildBlockedOperator();
  factory.buildPreconditionerOperator(blo, blkState);
  TEST_ASSERT(!Teuchos::is_null(blkState.monolithicMatrix));
  TEST_COMPARE(monolithicMismatch(blo, *blkState.monolithicMatrix), <, tol);

  // New values with the same sparsity pattern are updated in place.
  auto previous = blkState.monolithicMatrix;
  blo           = buildBlockedOperator(3.0);
  factory.buildPreconditionerOperator(blo, blkState);
  TEST_EQUALITY(blkState.monolithicMatrix.get(), previous.get());
  TEST_COMPARE(monolithicMismatch(blo, *blkState.monolithicMatrix), <, tol);

  // A changed sparsity pattern forces a rebuild.
  blo = buildBlockedOperator(1.0, true);
  factory.buildPreconditionerOperator(blo, blkState);
  TEST_INEQUALITY(blkState.monolithicMatrix.get(), previous.get());
  TEST_COMPARE(monolithicMismatch(blo, *blkState.monolithicMatrix), <, tol);
}

TEUCHOS_UNIT_TEST(TekoMultiPhysBlockPreconditionerFactory, ApplyHonorsScalingAndColumnCount) {
  MueLu::TekoMultiPhysBlockPreconditionerFactory factory;
  initializeFactory(factory);
  auto state   = factory.buildPreconditionerState();
  auto blo     = buildBlockedOperator();
  auto M       = factory.buildPreconditionerOperator(blo, Teuchos::dyn_cast<Teko::BlockPreconditionerState>(*state));
  const ST tol = 1e-10;

  // Single column first, then two columns whose first column matches, to exercise scratch reallocation.
  auto r1 = Thyra::createMembers(M->domain(), 1);
  Thyra::randomize(-1.0, 1.0, r1.ptr());
  auto z1 = Thyra::createMembers(M->range(), 1);
  Thyra::apply(*M, Thyra::NOTRANS, *r1, z1.ptr());

  auto r2 = Thyra::createMembers(M->domain(), 2);
  Thyra::randomize(-1.0, 1.0, r2.ptr());
  Thyra::V_V(r2->col(0).ptr(), *r1->col(0));
  auto z2 = Thyra::createMembers(M->range(), 2);
  Thyra::apply(*M, Thyra::NOTRANS, *r2, z2.ptr());
  TEST_COMPARE(maxRelativeDifference(*z2->col(0), *z1->col(0)), <, tol);

  // z = alpha * M r + beta * z
  const ST alpha = 2.0, beta = 0.5;
  auto zInit = Thyra::createMembers(M->range(), 2);
  Thyra::randomize(-1.0, 1.0, zInit.ptr());
  auto z = Thyra::createMembers(M->range(), 2);
  Thyra::assign(z.ptr(), *zInit);
  Thyra::apply(*M, Thyra::NOTRANS, *r2, z.ptr(), alpha, beta);

  auto expected = Thyra::createMembers(M->range(), 2);
  Thyra::assign(expected.ptr(), *zInit);
  Thyra::scale(beta, expected.ptr());
  Thyra::update(alpha, *z2, expected.ptr());
  TEST_COMPARE(maxRelativeDifference(*z, *expected), <, tol);

  // Pure scaling with beta == 0.
  Thyra::apply(*M, Thyra::NOTRANS, *r2, z.ptr(), alpha, 0.0);
  Thyra::assign(expected.ptr(), 0.0);
  Thyra::update(alpha, *z2, expected.ptr());
  TEST_COMPARE(maxRelativeDifference(*z, *expected), <, tol);
}

}  // namespace MueLuTests

#endif