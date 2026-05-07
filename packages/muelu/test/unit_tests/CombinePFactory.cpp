// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include "MueLu_CombinePFactory.hpp"
#include "MueLu_Exceptions.hpp"

#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include <Thyra_DefaultBlockedLinearOp.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CombinePFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);

  RCP<CombinePFactory> combinePFactory = rcp(new CombinePFactory());
  TEST_EQUALITY(combinePFactory != Teuchos::null, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CombinePFactory, CombineTwoSubblocks, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  const GO numElements0 = 200;
  const GO numElements1 = 200;
  const GO numElements  = numElements0 + numElements1;

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(1);

  // Build two subblock maps with contiguous global IDs:
  // block0: [0, ..., numElements0-1]
  // block1: [numElements0, ..., numElements-1]
  RCP<const Map> map0 = StridedMapFactory::Build(lib, numElements0, 0, stridingInfo, comm);
  RCP<const Map> map1 = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm);

  // Combined fine-level map
  RCP<const Map> bigMap = StridedMapFactory::Build(lib, numElements, 0, stridingInfo, comm);

  // Fine-level A only needs a compatible row map and communicator. A simple tridiagonal
  // on the combined map is sufficient.
  RCP<CrsMatrixWrap> A = Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(bigMap, bigMap->getGlobalNumElements(), 2.0, -1.0, -1.0);

  // Build two subblock prolongators. We use square matrices here for simplicity;
  // CombinePFactory only cares that they are valid matrices with compatible maps.
  RCP<CrsMatrixWrap> P0 = Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(map0, map0->getGlobalNumElements(), 2.0, -1.0, -1.0);

  RCP<CrsMatrixWrap> P1 = Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(map1, map1->getGlobalNumElements(), 3.0, -2.0, -1.0);

  TEST_EQUALITY(P0 != Teuchos::null, true);
  TEST_EQUALITY(P1 != Teuchos::null, true);

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);
  coarseLevel.SetFactoryManager(Teuchos::null);

  fineLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(A));
  coarseLevel.Set("Psubblock0", Teuchos::rcp_dynamic_cast<Matrix>(P0));
  coarseLevel.Set("Psubblock1", Teuchos::rcp_dynamic_cast<Matrix>(P1));

  RCP<CombinePFactory> PFact = rcp(new CombinePFactory());
  Teuchos::ParameterList params;
  params.set("combine: numBlks", 2);
  params.set("combine: useMaxLevels", false);
  PFact->SetParameterList(params);

  coarseLevel.Request("P", PFact.get());
  TEST_EQUALITY(coarseLevel.IsRequested("P", PFact.get()), true);

  RCP<Matrix> P = coarseLevel.Get<RCP<Matrix> >("P", PFact.get());
  TEST_EQUALITY(P != Teuchos::null, true);

  // Structural checks
  TEST_EQUALITY(P->getRowMap()->getGlobalNumElements(), numElements);
  TEST_EQUALITY(P->getDomainMap()->getGlobalNumElements(),
                P0->getDomainMap()->getGlobalNumElements() + P1->getDomainMap()->getGlobalNumElements());
  TEST_EQUALITY(P->getGlobalNumEntries(), P0->getGlobalNumEntries() + P1->getGlobalNumEntries());

  // Functional check: P * [1;1] should equal [P0*1; P1*1]
  RCP<Vector> x0 = VectorFactory::Build(P0->getDomainMap());
  RCP<Vector> x1 = VectorFactory::Build(P1->getDomainMap());
  x0->putScalar(1.0);
  x1->putScalar(1.0);

  RCP<Vector> y0 = VectorFactory::Build(P0->getRangeMap());
  RCP<Vector> y1 = VectorFactory::Build(P1->getRangeMap());
  P0->apply(*x0, *y0);
  P1->apply(*x1, *y1);

  RCP<Vector> x = VectorFactory::Build(P->getDomainMap());
  x->putScalar(1.0);

  RCP<Vector> y = VectorFactory::Build(P->getRangeMap());
  P->apply(*x, *y);

  // Since domain/range maps are concatenated in the same block order, the 1-norm
  // and infinity-norm of the combined result should match the blockwise assembled result.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType expectedNorm1 =
      y0->norm1() + y1->norm1();
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType actualNorm1 = y->norm1();

  TEST_EQUALITY(actualNorm1, expectedNorm1);
  TEST_COMPARE(actualNorm1, ==, expectedNorm1);

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType expectedInf =
      std::max(y0->normInf(), y1->normInf());
  TEST_EQUALITY(y->normInf(), expectedInf);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CombinePFactory, CombineSingleBlockDegenerate, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  const GO numElements = 200;

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(1);

  RCP<const Map> map0 = StridedMapFactory::Build(lib, numElements, 0, stridingInfo, comm);

  RCP<CrsMatrixWrap> A  = Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(map0, map0->getGlobalNumElements(), 2.0, -1.0, -1.0);
  RCP<CrsMatrixWrap> P0 = Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(map0, map0->getGlobalNumElements(), 3.0, -2.0, -1.0);

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);
  coarseLevel.SetFactoryManager(Teuchos::null);

  fineLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(A));
  coarseLevel.Set("Psubblock0", Teuchos::rcp_dynamic_cast<Matrix>(P0));

  RCP<CombinePFactory> PFact = rcp(new CombinePFactory());
  Teuchos::ParameterList params;
  params.set("combine: numBlks", 1);
  params.set("combine: useMaxLevels", false);
  PFact->SetParameterList(params);

  coarseLevel.Request("P", PFact.get());

  RCP<Matrix> P = coarseLevel.Get<RCP<Matrix> >("P", PFact.get());
  TEST_EQUALITY(P != Teuchos::null, true);

  TEST_EQUALITY(P->getRowMap()->getGlobalNumElements(), P0->getRowMap()->getGlobalNumElements());
  TEST_EQUALITY(P->getDomainMap()->getGlobalNumElements(), P0->getDomainMap()->getGlobalNumElements());
  TEST_EQUALITY(P->getGlobalNumEntries(), P0->getGlobalNumEntries());

  RCP<Vector> x = VectorFactory::Build(P0->getDomainMap());
  x->putScalar(1.0);

  RCP<Vector> yRef = VectorFactory::Build(P0->getRangeMap());
  RCP<Vector> y    = VectorFactory::Build(P->getRangeMap());

  P0->apply(*x, *yRef);
  P->apply(*x, *y);

  TEST_EQUALITY(y->norm1(), yRef->norm1());
  TEST_EQUALITY(y->normInf(), yRef->normInf());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CombinePFactory, CombineWithIdentityFallback, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  const GO numElements0 = 200;
  const GO numElements1 = 200;
  const GO numElements  = numElements0 + numElements1;

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(1);

  RCP<const Map> map0 = StridedMapFactory::Build(lib, numElements0, 0, stridingInfo, comm);
  RCP<const Map> map1 = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm);

  RCP<const Map> bigMap = StridedMapFactory::Build(lib, numElements, 0, stridingInfo, comm);

  RCP<CrsMatrixWrap> Acombo = Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(bigMap, bigMap->getGlobalNumElements(), 2.0, -1.0, -1.0);
  RCP<CrsMatrixWrap> P0     = Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(map0, map0->getGlobalNumElements(), 2.0, -1.0, -1.0);
  RCP<CrsMatrixWrap> A1     = Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(map1, map1->getGlobalNumElements(), 3.0, -2.0, -1.0);

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);
  coarseLevel.SetFactoryManager(Teuchos::null);

  fineLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(Acombo));
  coarseLevel.Set("Psubblock0", Teuchos::rcp_dynamic_cast<Matrix>(P0));

  // Provide Operatorsubblock1 but intentionally omit Psubblock1, so that
  // CombinePFactory constructs an identity prolongator for block 1.
  coarseLevel.Set("Operatorsubblock1",
                  Teuchos::rcp_dynamic_cast<Operator>(A1));

  RCP<CombinePFactory> PFact = rcp(new CombinePFactory());
  Teuchos::ParameterList params;
  params.set("combine: numBlks", 2);
  params.set("combine: useMaxLevels", true);
  PFact->SetParameterList(params);

  coarseLevel.Request("P", PFact.get());

  RCP<Matrix> P = coarseLevel.Get<RCP<Matrix> >("P", PFact.get());
  TEST_EQUALITY(P != Teuchos::null, true);

  TEST_EQUALITY(P->getRowMap()->getGlobalNumElements(), numElements);
  TEST_EQUALITY(P->getDomainMap()->getGlobalNumElements(),
                P0->getDomainMap()->getGlobalNumElements() + A1->getDomainMap()->getGlobalNumElements());

  // Action check with ones:
  // y = [P0*1; I*1] = [P0*1; 1]
  RCP<Vector> x = VectorFactory::Build(P->getDomainMap());
  x->putScalar(1.0);

  RCP<Vector> y = VectorFactory::Build(P->getRangeMap());
  P->apply(*x, *y);

  RCP<Vector> y0   = VectorFactory::Build(P0->getRangeMap());
  RCP<Vector> one1 = VectorFactory::Build(A1->getRangeMap());
  RCP<Vector> x0   = VectorFactory::Build(P0->getDomainMap());

  x0->putScalar(1.0);
  P0->apply(*x0, *y0);
  one1->putScalar(1.0);

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType expectedNorm1 =
      y0->norm1() + one1->norm1();
  TEST_EQUALITY(y->norm1(), expectedNorm1);
  TEST_COMPARE(y->normInf(), >=, one1->normInf());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CombinePFactory, CombineWithBlockedFineMatrix, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = MueLuTests::TestHelpers::Parameters::getLib();

  const GO numElements0 = 200;
  const GO numElements1 = 200;

  std::vector<size_t> stridingInfo(1, 1);

  // CombinePFactory expects each Psubblock to use index base 0
  RCP<const Map> subMap0 = StridedMapFactory::Build(lib, numElements0, 0, stridingInfo, comm);
  RCP<const Map> subMap1 = StridedMapFactory::Build(lib, numElements1, 0, stridingInfo, comm);

  // Fine blocked operator uses contiguous global IDs across its blocks
  RCP<const Map> rangeMap0 = StridedMapFactory::Build(lib, numElements0, 0, stridingInfo, comm);
  RCP<const Map> rangeMap1 = StridedMapFactory::Build(lib, numElements1, numElements0, stridingInfo, comm);

  // Build blocked fine operator subblocks
  RCP<CrsMatrixWrap> A00 =
      Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(rangeMap0, rangeMap0->getGlobalNumElements(), 2.0, -1.0, -1.0);
  RCP<CrsMatrixWrap> A01 =
      Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(rangeMap0, rangeMap1->getGlobalNumElements(), 1.0, 0.0, 0.0);
  RCP<CrsMatrixWrap> A10 =
      Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(rangeMap1, rangeMap0->getGlobalNumElements(), 1.0, 0.0, 0.0);
  RCP<CrsMatrixWrap> A11 =
      Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(rangeMap1, rangeMap1->getGlobalNumElements(), 3.0, -2.0, -1.0);

  TEST_EQUALITY(A00 != Teuchos::null, true);
  TEST_EQUALITY(A01 != Teuchos::null, true);
  TEST_EQUALITY(A10 != Teuchos::null, true);
  TEST_EQUALITY(A11 != Teuchos::null, true);

  RCP<BlockedCrsMatrix> bA;

#ifdef HAVE_XPETRA_THYRA
  using TpetraLO   = typename Tpetra::Map<>::local_ordinal_type;
  using TpetraGO   = typename Tpetra::Map<>::global_ordinal_type;
  using TpetraNode = typename Tpetra::Map<>::node_type;

  constexpr bool typesMatch =
      std::is_same_v<LocalOrdinal, TpetraLO> &&
      std::is_same_v<GlobalOrdinal, TpetraGO> &&
      std::is_same_v<Node, TpetraNode> &&
      std::is_same_v<Scalar, double>;

  if (typesMatch) {
    RCP<Thyra::LinearOpBase<Scalar> > thA00 =
        Xpetra::ThyraUtils<Scalar, LO, GO, Node>::toThyra(A00->getCrsMatrix());
    RCP<Thyra::LinearOpBase<Scalar> > thA01 =
        Xpetra::ThyraUtils<Scalar, LO, GO, Node>::toThyra(A01->getCrsMatrix());
    RCP<Thyra::LinearOpBase<Scalar> > thA10 =
        Xpetra::ThyraUtils<Scalar, LO, GO, Node>::toThyra(A10->getCrsMatrix());
    RCP<Thyra::LinearOpBase<Scalar> > thA11 =
        Xpetra::ThyraUtils<Scalar, LO, GO, Node>::toThyra(A11->getCrsMatrix());

    TEST_EQUALITY(thA00 != Teuchos::null, true);
    TEST_EQUALITY(thA01 != Teuchos::null, true);
    TEST_EQUALITY(thA10 != Teuchos::null, true);
    TEST_EQUALITY(thA11 != Teuchos::null, true);

    RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar> > thyraBlockedOp =
        Thyra::defaultBlockedLinearOp<Scalar>();

    thyraBlockedOp->beginBlockFill(2, 2);
    thyraBlockedOp->setNonconstBlock(0, 0, thA00);
    thyraBlockedOp->setNonconstBlock(0, 1, thA01);
    thyraBlockedOp->setNonconstBlock(1, 0, thA10);
    thyraBlockedOp->setNonconstBlock(1, 1, thA11);
    thyraBlockedOp->endBlockFill();

    RCP<const Thyra::BlockedLinearOpBase<Scalar> > thyraBlockedOpConst = thyraBlockedOp;

    bA = Teuchos::rcp(new BlockedCrsMatrix(thyraBlockedOpConst, comm));
    bA->fillComplete();
  }
#endif

#ifndef HAVE_XPETRA_THYRA
  // Without Thyra support, the fine matrix cannot be constructed through the Thyra-style constructor.
  bool threw = false;
  try {
    Level fineLevel, coarseLevel;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    fineLevel.SetFactoryManager(Teuchos::null);
    coarseLevel.SetFactoryManager(Teuchos::null);

    fineLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bA));
  } catch (const std::exception& e) {
    threw = true;
    out << "Caught expected exception without Xpetra/Thyra support: " << e.what() << std::endl;
  }
  TEST_EQUALITY(threw || bA == Teuchos::null, true);
  return;
#else
  if (!typesMatch) {
    TEST_EQUALITY(bA == Teuchos::null, true);
    return;
  }
#endif

  TEST_EQUALITY(bA != Teuchos::null, true);
  TEST_EQUALITY(bA->Rows(), 2);
  TEST_EQUALITY(bA->Cols(), 2);

  // Valid subblock prolongators
  RCP<CrsMatrixWrap> P0 =
      Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(subMap0, subMap0->getGlobalNumElements(), 2.0, -1.0, -1.0);
  RCP<CrsMatrixWrap> P1 =
      Galeri::Xpetra::TriDiag<Scalar, LO, GO, Map, CrsMatrixWrap>(subMap1, subMap1->getGlobalNumElements(), 3.0, -2.0, -1.0);

  TEST_EQUALITY(P0 != Teuchos::null, true);
  TEST_EQUALITY(P1 != Teuchos::null, true);

  Level fineLevel, coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);
  coarseLevel.SetFactoryManager(Teuchos::null);

  fineLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bA));
  coarseLevel.Set("Psubblock0", Teuchos::rcp_dynamic_cast<Matrix>(P0));
  coarseLevel.Set("Psubblock1", Teuchos::rcp_dynamic_cast<Matrix>(P1));

  RCP<CombinePFactory> PFact = rcp(new CombinePFactory());
  Teuchos::ParameterList params;
  params.set("combine: numBlks", 2);
  params.set("combine: useMaxLevels", false);
  PFact->SetParameterList(params);

  coarseLevel.Request("P", PFact.get());

  bool threw = false;
  RCP<Matrix> P;
  try {
    P = coarseLevel.Get<RCP<Matrix> >("P", PFact.get());
  } catch (const std::exception& e) {
    threw = true;
  }

  TEST_EQUALITY(threw, false);
  TEST_EQUALITY(P != Teuchos::null, true);

  RCP<BlockedCrsMatrix> bP = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(P);
  TEST_EQUALITY(bP != Teuchos::null, true);

  TEST_EQUALITY(bP->Rows(), 2);
  TEST_EQUALITY(bP->Cols(), 2);

  // Since BuildPBlocked sets only diagonal blocks, verify off-diagonals are null/empty
  TEST_EQUALITY(bP->getMatrix(0, 0) != Teuchos::null, true);
  TEST_EQUALITY(bP->getMatrix(1, 1) != Teuchos::null, true);

  // Off-diagonal operators are stored as 0 blocks
  RCP<Matrix> P01 = bP->getMatrix(0, 1);
  RCP<Matrix> P10 = bP->getMatrix(1, 0);

  TEST_EQUALITY(P01 != Teuchos::null, true);
  TEST_EQUALITY(P10 != Teuchos::null, true);

  TEST_EQUALITY(P01->getGlobalNumEntries(), 0);
  TEST_EQUALITY(P10->getGlobalNumEntries(), 0);
  TEST_EQUALITY(P01->getLocalNumEntries(), 0);
  TEST_EQUALITY(P10->getLocalNumEntries(), 0);

  // Structural checks on the two diagonal blocks
  TEST_EQUALITY(bP->getMatrix(0, 0)->getGlobalNumRows(), P0->getGlobalNumRows());
  TEST_EQUALITY(bP->getMatrix(1, 1)->getGlobalNumRows(), P1->getGlobalNumRows());
  TEST_EQUALITY(bP->getMatrix(0, 0)->getGlobalNumEntries(), P0->getGlobalNumEntries());
  TEST_EQUALITY(bP->getMatrix(1, 1)->getGlobalNumEntries(), P1->getGlobalNumEntries());
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CombinePFactory, Constructor, SC, LO, GO, Node)                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CombinePFactory, CombineTwoSubblocks, SC, LO, GO, Node)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CombinePFactory, CombineSingleBlockDegenerate, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CombinePFactory, CombineWithIdentityFallback, SC, LO, GO, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CombinePFactory, CombineWithBlockedFineMatrix, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
