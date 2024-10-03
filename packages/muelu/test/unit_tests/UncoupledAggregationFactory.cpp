// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_CoarseMapFactory.hpp>
#include <MueLu_config.hpp>
#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_UncoupledAggregationFactory.hpp>
#include <MueLu_BlockedCoarseMapFactory.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UncoupledAggregationFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>

  out << "version: " << MueLu::Version() << std::endl;
  RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
  TEST_EQUALITY(aggFact != Teuchos::null, true);
  TEST_THROW(aggFact->SetOrdering("unknown_ordering"), Teuchos::Exceptions::InvalidParameterValue);
  aggFact->SetOrdering("natural");
  TEST_EQUALITY(aggFact->GetOrdering() == "natural", true);
  aggFact->SetOrdering("graph");
  TEST_EQUALITY(aggFact->GetOrdering() == "graph", true);
  aggFact->SetOrdering("random");
  TEST_EQUALITY(aggFact->GetOrdering() == "random", true);

  aggFact->SetMaxNeighAlreadySelected(100);
  TEST_EQUALITY(aggFact->GetMaxNeighAlreadySelected() == 100, true);
  aggFact->SetMinNodesPerAggregate(10);
  TEST_EQUALITY(aggFact->GetMinNodesPerAggregate() == 10, true);

  auto validParamList = aggFact->GetValidParameterList();
  TEST_EQUALITY(validParamList != Teuchos::null, true);

}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UncoupledAggregationFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>

  out << "version: " << MueLu::Version() << std::endl;
  Teuchos::ParameterList galeriList;
  GlobalOrdinal nx = 20, ny = 20;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<const Map> map                  = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(TestHelpers::Parameters::getLib(), "Cartesian2D", comm, galeriList);
  map                                 = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, 2);  // expand map for 2 DOFs per node
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", map, galeriList);
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(2);

  Level aLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(aLevel);
  aLevel.Request("A");
  aLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());
  aggFact->SetFactory("Graph", dropFact);

  aLevel.Request(*aggFact);

  aggFact->Build(aLevel);

}  // Build
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UncoupledAggregationFactory, Build_ML, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>

  out << "version: " << MueLu::Version() << std::endl;
  Teuchos::ParameterList galeriList;
  GlobalOrdinal nx = 20, ny = 20;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<const Map> map                  = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(TestHelpers::Parameters::getLib(), "Cartesian2D", comm, galeriList);
  map                                 = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, 2);  // expand map for 2 DOFs per node
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", map, galeriList);
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(2);

  Level aLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(aLevel);
  aLevel.Request("A");
  aLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  RCP<UncoupledAggregationFactory> aggFact = rcp(new UncoupledAggregationFactory());

  // Test the ML style options
  aggFact->SetParameter("aggregation: match ML phase2a", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: match ML phase2b", Teuchos::ParameterEntry(true));
  aggFact->SetFactory("Graph", dropFact);

  aLevel.Request(*aggFact);

  aggFact->Build(aLevel);

}  // Build

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UncoupledAggregationFactory, Constructor, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UncoupledAggregationFactory, Build, Scalar, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UncoupledAggregationFactory, Build_ML, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
