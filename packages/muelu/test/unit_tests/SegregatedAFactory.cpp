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
#include <Teuchos_FancyOStream.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_SegregatedAFactory.hpp"
#include "MueLu_Utilities.hpp"

#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_IO.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SegregatedAFactory, Blockmap, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(SC, GO, NO);
  out << "version: " << MueLu::Version() << std::endl;

  using TST            = Teuchos::ScalarTraits<SC>;
  using magnitude_type = typename TST::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = MueLuTests::TestHelpers::Parameters::getLib();

  // test with simple example matrix and segregate first 4 rows from rest
  //     x x x x                x x x x
  //     x x x x                x x x x
  //     x x x x x x            x x x x
  //     x x x x x x      ==>   x x x x
  //         x x x x x                 x x x
  //         x x x x x                 x x x
  //             x x x                 x x x
  {
    RCP<Map> rowMap     = MapFactory::Build(lib, 7, 0, comm);
    RCP<CrsGraph> graph = CrsGraphFactory::Build(rowMap, 6);

    graph->insertGlobalIndices(0, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3));
    graph->insertGlobalIndices(1, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3));
    graph->insertGlobalIndices(2, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5));
    graph->insertGlobalIndices(3, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5));
    graph->insertGlobalIndices(4, Teuchos::tuple<GlobalOrdinal>(2, 3, 4, 5, 6));
    graph->insertGlobalIndices(5, Teuchos::tuple<GlobalOrdinal>(2, 3, 4, 5, 6));
    graph->insertGlobalIndices(6, Teuchos::tuple<GlobalOrdinal>(4, 5, 6));
    graph->fillComplete();

    RCP<Matrix> A = MatrixFactory::Build(graph.getConst());
    A->setAllToScalar(TST::one());
    A->fillComplete();

    TEST_EQUALITY(A.is_null(), false);
    TEST_EQUALITY(A->getGlobalNumEntries(), 33);
    TEST_FLOATING_EQUALITY(A->getFrobeniusNorm(), 5.744562646538029, 2e1 * TMT::eps());

    RCP<const Map> map = MapFactory::Build(lib, 4, 0, comm);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    std::string mapName = "dropMap1";
    level.Set(mapName, map);

    TEST_ASSERT(level.IsAvailable("A", MueLu::NoFactory::get()));
    TEST_ASSERT(level.IsAvailable(mapName, MueLu::NoFactory::get()));

    RCP<SegregatedAFactory> segregateFact = rcp(new SegregatedAFactory());
    segregateFact->SetFactory("A", MueLu::NoFactory::getRCP());
    segregateFact->SetParameter("droppingScheme", Teuchos::ParameterEntry(std::string("blockmap")));
    segregateFact->SetFactory("dropMap1", MueLu::NoFactory::getRCP());

    // request segregated operator
    level.Request("A", segregateFact.get());
    segregateFact->Build(level);

    RCP<Matrix> Aout = level.Get<RCP<Matrix>>("A", segregateFact.get());

    // Output (segregated Operator)
    // Node ID   Global Row  Num Entries(Index,Value)
    // 0         0           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         1           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         2           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         3           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         4           3 (4, 1)  (5, 1)  (6, 1)
    // 0         5           3 (4, 1)  (5, 1)  (6, 1)
    // 0         6           3 (4, 1)  (5, 1)  (6, 1)

    TEST_EQUALITY(Aout.is_null(), false);
    TEST_EQUALITY(Aout->getGlobalNumEntries(), 25);
    TEST_FLOATING_EQUALITY(Aout->getFrobeniusNorm(), 5.0, 2e1 * TMT::eps());
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SegregatedAFactory, MapPair, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(SC, GO, NO);
  out << "version: " << MueLu::Version() << std::endl;

  using TST            = Teuchos::ScalarTraits<SC>;
  using magnitude_type = typename TST::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = MueLuTests::TestHelpers::Parameters::getLib();

  // test with simple example matrix and segregate first 4 rows from rest
  //     x x x x                x x x x
  //     x x x x                x x x x
  //     x x x x x x            x x x x
  //     x x x x x x      ==>   x x x x
  //         x x x x x                 x x x
  //         x x x x x                 x x x
  //             x x x                 x x x

  {
    RCP<Map> rowMap     = MapFactory::Build(lib, 7, 0, comm);
    RCP<CrsGraph> graph = CrsGraphFactory::Build(rowMap, 6);

    graph->insertGlobalIndices(0, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3));
    graph->insertGlobalIndices(1, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3));
    graph->insertGlobalIndices(2, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5));
    graph->insertGlobalIndices(3, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5));
    graph->insertGlobalIndices(4, Teuchos::tuple<GlobalOrdinal>(2, 3, 4, 5, 6));
    graph->insertGlobalIndices(5, Teuchos::tuple<GlobalOrdinal>(2, 3, 4, 5, 6));
    graph->insertGlobalIndices(6, Teuchos::tuple<GlobalOrdinal>(4, 5, 6));
    graph->fillComplete();

    RCP<Matrix> A = MatrixFactory::Build(graph.getConst());
    A->setAllToScalar(TST::one());
    A->fillComplete();

    TEST_EQUALITY(A.is_null(), false);
    TEST_EQUALITY(A->getGlobalNumEntries(), 33);
    TEST_FLOATING_EQUALITY(A->getFrobeniusNorm(), 5.744562646538029, 2e1 * TMT::eps());

    std::vector<GO> dropMap1Entries{2, 3};
    RCP<const Map> dropMap1 = Xpetra::MapFactory<LO, GO, NO>::Build(lib, dropMap1Entries.size(), dropMap1Entries, 0, comm);
    std::vector<GO> dropMap2Entries{4, 5};
    RCP<const Map> dropMap2 = Xpetra::MapFactory<LO, GO, NO>::Build(lib, dropMap2Entries.size(), dropMap2Entries, 0, comm);

    Level level;
    TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    std::string map1Name = "dropMap1";
    std::string map2Name = "dropMap2";

    level.Set(map1Name, dropMap1);
    level.Set(map2Name, dropMap2);

    TEST_ASSERT(level.IsAvailable("A", MueLu::NoFactory::get()));
    TEST_ASSERT(level.IsAvailable(map1Name, MueLu::NoFactory::get()));
    TEST_ASSERT(level.IsAvailable(map2Name, MueLu::NoFactory::get()));

    RCP<SegregatedAFactory> segregateFact = rcp(new SegregatedAFactory());
    segregateFact->SetFactory("A", MueLu::NoFactory::getRCP());
    segregateFact->SetParameter("droppingScheme", Teuchos::ParameterEntry(std::string("map-pair")));
    segregateFact->SetFactory("dropMap1", MueLu::NoFactory::getRCP());
    segregateFact->SetFactory("dropMap2", MueLu::NoFactory::getRCP());

    // request segregated operator
    level.Request("A", segregateFact.get());
    segregateFact->Build(level);

    RCP<Matrix> Aout = level.Get<RCP<Matrix>>("A", segregateFact.get());

    // Output (segregated Operator)
    // Node ID   Global Row  Num Entries(Index,Value)
    // 0         0           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         1           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         2           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         3           4 (0, 1)  (1, 1)  (2, 1)  (3, 1)
    // 0         4           3 (4, 1)  (5, 1)  (6, 1)
    // 0         5           3 (4, 1)  (5, 1)  (6, 1)
    // 0         6           3 (4, 1)  (5, 1)  (6, 1)

    TEST_EQUALITY(Aout.is_null(), false);
    TEST_EQUALITY(Aout->getGlobalNumEntries(), 25);
    TEST_FLOATING_EQUALITY(Aout->getFrobeniusNorm(), 5.0, 2e1 * TMT::eps());
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SegregatedAFactory, importOffRankDroppingInfo, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  int numProcs                       = comm->getSize();
  int myRank                         = comm->getRank();
  Xpetra::UnderlyingLib lib          = MueLuTests::TestHelpers::Parameters::getLib();

  {
    // Create different map entries on each processor
    size_t globalDropMapSize       = 3;
    std::vector<GO> dropMapEntries = {};
    if (numProcs == 1) {
      dropMapEntries = {1, 7, 8};
    } else if (numProcs == 4) {
      if (myRank == 0) {
        dropMapEntries = {7, 8};
      } else if (myRank == 1) {
        dropMapEntries = {};
      } else if (myRank == 2) {
        dropMapEntries = {1};
      } else if (myRank == 3) {
        dropMapEntries = {};
      }
    } else {
      std::cout << "\nThis test was designed to run on exactly 1 or 4 ranks!\n"
                << std::endl;
      return;
    }
    RCP<const Map> dropMap = Xpetra::MapFactory<LO, GO, NO>::Build(lib, globalDropMapSize, dropMapEntries, 0, comm);

    int graphSize       = 13;
    RCP<Map> rowMap     = MapFactory::Build(lib, graphSize, 0, comm);
    RCP<CrsGraph> graph = CrsGraphFactory::Build(rowMap, graphSize);

    // Create graph for symmetric matrix with off-Rank entries
    //
    //     x x x x x x
    //     x x x x x x   x x
    //     x x x x x x       x x
    //     x x x x x x   x
    //     x x x x x x
    //     x x x x x x x   x     x
    //               x x x x x x x
    //       x   x     x x x x x x
    //       x       x x x x x x x
    //         x       x x x x x x
    //         x       x x x x x x
    //               x x x x x x x
    //                             x

    graph->insertGlobalIndices(0, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5));
    graph->insertGlobalIndices(1, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5, 7, 8));
    graph->insertGlobalIndices(2, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5, 9, 10));
    graph->insertGlobalIndices(3, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5, 7));
    graph->insertGlobalIndices(4, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5));
    graph->insertGlobalIndices(5, Teuchos::tuple<GlobalOrdinal>(0, 1, 2, 3, 4, 5, 6, 8, 11));
    graph->insertGlobalIndices(6, Teuchos::tuple<GlobalOrdinal>(5, 6, 7, 8, 9, 10, 11));
    graph->insertGlobalIndices(7, Teuchos::tuple<GlobalOrdinal>(1, 3, 6, 7, 8, 9, 10, 11));
    graph->insertGlobalIndices(8, Teuchos::tuple<GlobalOrdinal>(1, 5, 6, 7, 8, 9, 10, 11));
    graph->insertGlobalIndices(9, Teuchos::tuple<GlobalOrdinal>(2, 6, 7, 8, 9, 10, 11));
    graph->insertGlobalIndices(10, Teuchos::tuple<GlobalOrdinal>(2, 6, 7, 8, 9, 10, 11));
    graph->insertGlobalIndices(11, Teuchos::tuple<GlobalOrdinal>(5, 6, 7, 8, 9, 10, 11));
    graph->insertGlobalIndices(12, Teuchos::tuple<GlobalOrdinal>(12));
    graph->fillComplete();

    RCP<Matrix> A = MatrixFactory::Build(graph.getConst());
    A->setAllToScalar(Teuchos::ScalarTraits<SC>::one());
    A->fillComplete();

    Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> finalDropMap = MueLu::importOffRankDroppingInfo<SC, LO, GO, NO>(dropMap, A);

    size_t globalExpectedMapSize       = 0;
    std::vector<GO> expectedMapEntries = {};

    if (numProcs == 1) {
      globalExpectedMapSize = 3;
      expectedMapEntries    = {1, 7, 8};
    } else if (numProcs == 4) {
      globalExpectedMapSize = 11;
      if (myRank == 0) {
        expectedMapEntries = {1, 7, 8};
      } else if (myRank == 1) {
        expectedMapEntries = {1, 7, 8};
      } else if (myRank == 2) {
        expectedMapEntries = {7, 8, 1};
      } else if (myRank == 3) {
        expectedMapEntries = {7, 8};
      }
    }
    RCP<const Map> expectedMap = Xpetra::MapFactory<LO, GO, NO>::Build(lib, globalExpectedMapSize, expectedMapEntries, 0, comm);
    TEUCHOS_ASSERT(expectedMap->isSameAs(*finalDropMap));
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SegregatedAFactory, Blockmap, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SegregatedAFactory, MapPair, SC, LO, GO, NO)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SegregatedAFactory, importOffRankDroppingInfo, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>
}  // namespace MueLuTests
