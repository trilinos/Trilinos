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

#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>

#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_Aggregates.hpp>
#include <MueLu_BlockedCoarseMapFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedCoarseMapFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<BlockedCoarseMapFactory> mapFact = rcp(new BlockedCoarseMapFactory());
  TEST_EQUALITY(mapFact != Teuchos::null, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedCoarseMapFactory, BuildBlockedCoarseMapWithGIDOffset, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel;
  Level coarseLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);
  coarseLevel.SetFactoryManager(Teuchos::null);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(29);
  A->SetFixedBlockSize(1);
  fineLevel.Set("A", A);

  LO NSdim                   = 2;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->randomize();
  fineLevel.Set("Nullspace", nullSpace);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<UncoupledAggregationFactory> UncoupledAggFact = rcp(new UncoupledAggregationFactory());
  UncoupledAggFact->SetFactory("Graph", dropFact);
  UncoupledAggFact->SetFactory("DofsPerNode", dropFact);

  UncoupledAggFact->SetMinNodesPerAggregate(3);
  UncoupledAggFact->SetMaxNeighAlreadySelected(0);
  UncoupledAggFact->SetOrdering("natural");

  RCP<CoarseMapFactory> coarseMapFact = rcp(new CoarseMapFactory());
  coarseMapFact->SetFactory("Aggregates", UncoupledAggFact);

  RCP<BlockedCoarseMapFactory> blockedCoarseMapFact = rcp(new BlockedCoarseMapFactory());
  blockedCoarseMapFact->SetFactory("Aggregates", UncoupledAggFact);
  blockedCoarseMapFact->SetFactory("CoarseMap", coarseMapFact);

  // request input for BlockedCoarseMapFactory by hand
  fineLevel.Request("Aggregates", UncoupledAggFact.get());
  fineLevel.Request("CoarseMap", coarseMapFact.get());
  fineLevel.Request("CoarseMap", blockedCoarseMapFact.get());
  blockedCoarseMapFact->Build(fineLevel);
  RCP<const Map> map1 = fineLevel.Get<RCP<const Map>>("CoarseMap", coarseMapFact.get());
  RCP<const Map> map2 = fineLevel.Get<RCP<const Map>>("CoarseMap", blockedCoarseMapFact.get());

  // access aggregates
  RCP<Aggregates> aggregates         = fineLevel.Get<RCP<Aggregates>>("Aggregates", UncoupledAggFact.get());
  GO numAggs                         = aggregates->GetNumAggregates();
  GO numGlobalAggs                   = 0;
  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  MueLu_sumAll(comm, numAggs, numGlobalAggs);

  TEST_EQUALITY(map1->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(map1->getMaxAllGlobalIndex(), numGlobalAggs * Teuchos::as<GO>(NSdim) - 1);
  TEST_EQUALITY(map2->getMinAllGlobalIndex(), numGlobalAggs * Teuchos::as<GO>(NSdim));
  TEST_EQUALITY(map2->getMaxAllGlobalIndex(), 2 * numGlobalAggs * Teuchos::as<GO>(NSdim) - 1);
  TEST_EQUALITY(Teuchos::as<GO>(map1->getLocalNumElements()), numAggs * Teuchos::as<GO>(NSdim));
  TEST_EQUALITY(Teuchos::as<GO>(map2->getLocalNumElements()), numAggs * Teuchos::as<GO>(NSdim));
}

#define MUELU_ETI_GROUP(Scalar, LocalOrdinal, GlobalOrdinal, Node)                                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedCoarseMapFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedCoarseMapFactory, BuildBlockedCoarseMapWithGIDOffset, Scalar, LocalOrdinal, GlobalOrdinal, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
