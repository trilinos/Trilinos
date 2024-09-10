// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>

#include <MueLu_InitialBlockNumberFactory.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_PreDropFunctionConstVal.hpp>
#include <MueLu_LWGraph.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<CoalesceDropFactory> coalesceDropFact = rcp(new CoalesceDropFactory());
  coalesceDropFact->SetDefaultVerbLevel(MueLu::Extreme);
  TEST_EQUALITY(coalesceDropFact != Teuchos::null, true);

}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);
}  // Build

// TODO remove this
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, PreDrop, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(3);
  fineLevel.Set("A", A);
  A->describe(out, Teuchos::VERB_EXTREME);

  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetVerbLevel(MueLu::Extreme);
  dropFact.SetPreDropFunction(rcp(new PreDropFunctionConstVal(0.00001)));

  fineLevel.Request("Graph", &dropFact);

  dropFact.Build(fineLevel);

  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

  std::cout << graph->GetDomainMap()->getGlobalNumElements() << std::endl;
  graph->print(out, MueLu::Debug);

  //    TEST_EQUALITY(1 == 0, true);

}  // PreDrop

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationBasic, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 3.
  // lightweight wrap = false
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  int blockSize = 3;

  int nx        = blockSize * comm->getSize();
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx);
  A->SetFixedBlockSize(blockSize, 0);
  fineLevel.Set("A", A);
  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(false));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == blockSize, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == blockSize) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));  // u,v and p
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 2), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 1), true);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);
}  // AmalgamationBasic

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationStrided, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 3 using a strided map
  // lightweight wrap = false
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  int blockSize = 3;
  int nx        = blockSize * comm->getSize();
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx);

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(Teuchos::as<size_t>(blockSize));
  LocalOrdinal stridedBlockId = -1;

  RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      A->getRangeMap(),
      stridingInfo,
      stridedBlockId,
      0 /*offset*/
  );
  RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      A->getDomainMap(),
      stridingInfo,
      stridedBlockId,
      0 /*offset*/
  );

  if (A->IsView("stridedMaps") == true) A->RemoveView("stridedMaps");
  A->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

  fineLevel.Set("A", A);
  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(false));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);
  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == blockSize, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == blockSize) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 2), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 1), true);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);
}  // AmalgamationStrided

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationStrided2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 3 = (2,1). wrap block 0
  // lightweight wrap = false
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  // create strided map information
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(Teuchos::as<size_t>(2));
  stridingInfo.push_back(Teuchos::as<size_t>(1));
  LocalOrdinal stridedBlockId = 0;

  RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 3 * comm->getSize(), 0,
                                                                                                     stridingInfo, comm,
                                                                                                     stridedBlockId /*blockId*/, 0 /*offset*/);

  /////////////////////////////////////////////////////

  Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getRangeMap(),
      stridingInfo,
      stridedBlockId,
      0 /*offset*/
  );
  RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getDomainMap(),
      stridingInfo,
      stridedBlockId,
      0 /*offset*/
  );

  if (mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
  mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);  // should have holes in these maps

  fineLevel.Set("A", mtx);
  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(false));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

  LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 3, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == 3) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * 3 - 2), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * 3 - 1), true);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);
}  // AmalgamationStrided2

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationStridedOffset, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 9 = (2,3,4). wrap block 1.
  // lightweight wrap = false
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  // create strided map information
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(Teuchos::as<size_t>(2));
  stridingInfo.push_back(Teuchos::as<size_t>(3));
  stridingInfo.push_back(Teuchos::as<size_t>(4));
  LocalOrdinal stridedBlockId = 1;
  GlobalOrdinal offset        = 19;

  RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 9 * comm->getSize(), 0,
                                                                                                     stridingInfo, comm,
                                                                                                     stridedBlockId, offset);

  /////////////////////////////////////////////////////

  Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -3.0);

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<const Map> stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getRangeMap(),
      stridingInfo,
      stridedBlockId,
      offset);
  RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getDomainMap(),
      stridingInfo,
      stridedBlockId,
      offset);

  if (mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
  mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

  fineLevel.Set("A", mtx);
  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(false));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

  LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 9, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == 3) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * 3 - 2), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * 3 - 1), true);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);
}  // AmalgamationStridedOffset

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationLightweight, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 3
  // lightweight wrap = true
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  int blockSize = 3;

  RCP<const Map> dofMap    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, blockSize * comm->getSize(), 0, comm);
  Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);
  mtx->SetFixedBlockSize(blockSize, 0);
  fineLevel.Set("A", mtx);

  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == blockSize, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == blockSize) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 2), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 1), true);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);
}  // AmalgamationLightweight

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationLightweightDrop, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 1
  // lightweight wrap = true
  // drop small values
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<const Map> dofMap    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 3 * comm->getSize(), 0, comm);
  Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 1.0, -1.0, -0.0001);
  mtx->SetFixedBlockSize(1, 0);
  fineLevel.Set("A", mtx);

  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
  dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.5));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == 3 * comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0) {
      if (graph->getNeighborVertices(0).length == 1) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 3 * comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(3 * comm->getSize() + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 3 * comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(3 * comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 3), true);
}  // AmalgamationLightweightDrop

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationStridedLW, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 3 using a strided map
  // lightweight wrap = true
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  int blockSize = 3;

  int nx        = blockSize * comm->getSize();
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(nx);

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(Teuchos::as<size_t>(blockSize));
  LocalOrdinal stridedBlockId = -1;

  RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      A->getRangeMap(),
      stridingInfo,
      stridedBlockId,
      0 /*offset*/
  );
  RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      A->getDomainMap(),
      stridingInfo,
      stridedBlockId,
      0 /*offset*/
  );

  if (A->IsView("stridedMaps") == true) A->RemoveView("stridedMaps");
  A->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

  fineLevel.Set("A", A);
  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == blockSize, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == blockSize) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 2), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 1), true);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);
}  // AmalgamationStridedLW

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationStrided2LW, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 3 = (2,1). wrap block 0
  // lightweight wrap = true
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  // create strided map information
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(Teuchos::as<size_t>(2));
  stridingInfo.push_back(Teuchos::as<size_t>(1));
  LocalOrdinal stridedBlockId = 0;

  int blockSize = 3;

  RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, blockSize * comm->getSize(), 0,
                                                                                                     stridingInfo, comm,
                                                                                                     stridedBlockId /*blockId*/, 0 /*offset*/);

  /////////////////////////////////////////////////////

  Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getRangeMap(),
      stridingInfo,
      stridedBlockId,
      0 /*offset*/
  );
  RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getDomainMap(),
      stridingInfo,
      stridedBlockId,
      0 /*offset*/
  );
  if (mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
  mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

  fineLevel.Set("A", mtx);
  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

  LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == blockSize, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == blockSize) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 2), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(maxLocalIndex == numLocalRowMapElts * blockSize - 1), true);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);
}  // AmalgamationStrided2LW

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationStridedOffsetLW, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 9 = (2,3,4). wrap block 1.
  // lightweight wrap = true
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  // create strided map information
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(Teuchos::as<size_t>(2));
  stridingInfo.push_back(Teuchos::as<size_t>(3));
  stridingInfo.push_back(Teuchos::as<size_t>(4));
  LocalOrdinal stridedBlockId = 1;
  GlobalOrdinal offset        = 19;

  int blockSize = 9;

  RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, blockSize * comm->getSize(), 0,
                                                                                                     stridingInfo, comm,
                                                                                                     stridedBlockId, offset);

  /////////////////////////////////////////////////////

  Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -3.0);

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<const Map> stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getRangeMap(),
      stridingInfo,
      stridedBlockId,
      offset);
  RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getDomainMap(),
      stridingInfo,
      stridedBlockId,
      offset);

  if (mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
  mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

  fineLevel.Set("A", mtx);
  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

  LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == blockSize, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == 3) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);
}  // AmalgamationStridedOffsetLW

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationDroppingLW, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 3 = (3)
  // drop small entries
  // lightweight wrap = true
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<const Map> dofMap    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 3 * comm->getSize(), 0, comm);
  Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, 0.00001);
  mtx->SetFixedBlockSize(3, 0);
  fineLevel.Set("A", mtx);

  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
  dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0));
  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

  LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 3, true);
  TEST_EQUALITY(graph->getNeighborVertices(0).length, 1);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);

}  // AmalgamationDroppingLW

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AmalgamationStridedOffsetDropping2LW, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // unit test for block size 9 = (2,3,4). wrap block 1.
  // drop small entries
  // lightweight wrap = true
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  // create strided map information
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(Teuchos::as<size_t>(2));
  stridingInfo.push_back(Teuchos::as<size_t>(3));
  stridingInfo.push_back(Teuchos::as<size_t>(4));
  LocalOrdinal stridedBlockId = 1;
  GlobalOrdinal offset        = 19;

  RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 9 * comm->getSize(), 0,
                                                                                                     stridingInfo, comm,
                                                                                                     stridedBlockId, offset);

  /////////////////////////////////////////////////////

  Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, 1.0, 0.0001);

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<const Map> stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getRangeMap(),
      stridingInfo,
      stridedBlockId,
      offset);
  RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      mtx->getDomainMap(),
      stridingInfo,
      stridedBlockId,
      offset);

  if (mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
  mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

  fineLevel.Set("A", mtx);
  CoalesceDropFactory dropFact = CoalesceDropFactory();
  dropFact.SetDefaultVerbLevel(MueLu::Extreme);
  dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
  dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.3));

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  fineLevel.print(out);
  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

  LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
  TEST_EQUALITY(Teuchos::as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 9, true);
  bool bCorrectGraph = false;
  if (comm->getSize() == 1 && graph->getNeighborVertices(0).length == 1) {
    bCorrectGraph = true;
  } else {
    if (comm->getRank() == 0) {
      if (graph->getNeighborVertices(0).length == 1) bCorrectGraph = true;
    } else {
      if (graph->getNeighborVertices(0).length == 2) bCorrectGraph = true;
    }
  }
  TEST_EQUALITY(bCorrectGraph, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 1), true);
    } else {
      TEST_EQUALITY(Teuchos::as<bool>(numLocalImportElts == numLocalRowMapElts + 2), true);
    }
  }
  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), Teuchos::as<size_t>(comm->getSize()));
  TEST_EQUALITY(Teuchos::as<bool>(myDomainMap->getLocalNumElements() == 1), true);
}  // AmalgamationStridedOffsetDropping2LW

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, DistanceLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 40);

}  // DistanceLaplacian

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, DistanceLaplacianScaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  // Now we doctor the coordinates so that the off-diagonal pair row 0 will want to keep (0,1) and row 1 will want to drop (1,0)
  if (!comm->getRank()) {
    auto vals = coordinates->getDataNonConst(0);
    vals[0]   = vals[0] - 2000 * 36;
  }

  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.35));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(std::string("scaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // DistanceLaplacianCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, DistanceLaplacianUnscaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

  // Now we doctor the coordinates so that the off-diagonal pair row 0 will want to keep (0,1) and row 1 will want to drop (1,0)
  if (!comm->getRank()) {
    auto vals = coordinates->getDataNonConst(0);
    vals[0]   = vals[0] - 2000 * 36;
  }

  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.35));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(std::string("unscaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // DistanceLaplacianUnscaleCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, DistanceLaplacianCutSym, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

  // Now we doctor the coordinates so that the off-diagonal pair row 0 will want to keep (0,1) and row 1 will want to drop (1,0)
  if (!comm->getRank()) {
    auto vals = coordinates->getDataNonConst(0);
    vals[0]   = vals[0] - 2000 * 36;
  }

  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.35));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(std::string("scaled cut symmetric")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 106);

}  // DistanceLaplacianCutScaled

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, ClassicalScaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.1));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("scaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 72);

}  // SignaledClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, ClassicalUnScaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.1));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("unscaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 72);

}  // SignaledClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, SignaledClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // SignaledClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, SignaledScaledCutClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.1));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("scaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  TEST_THROW(coalesceDropFact.Build(fineLevel), MueLu::Exceptions::RuntimeError);

  //    RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  //    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  //    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  //    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  //    const RCP<const Map> myDomainMap = graph->GetDomainMap();

  //    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  //    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  //    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
  //    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(36 + (comm->getSize()-1)*2));

  //    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  //    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  //    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
  //    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),36);

  //    TEST_EQUALITY(graph->GetGlobalNumEdges(),36);

}  // SignaledScaledCutClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, SignaledUnscaledCutClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.1));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("unscaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  TEST_THROW(coalesceDropFact.Build(fineLevel), MueLu::Exceptions::RuntimeError);

  //    RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  //    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  //    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  //    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  //    const RCP<const Map> myDomainMap = graph->GetDomainMap();

  //    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  //    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  //    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
  //    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(36 + (comm->getSize()-1)*2));

  //    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  //    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  //    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
  //    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),36);

  //    TEST_EQUALITY(graph->GetGlobalNumEdges(),36);

}  // SignaledUnScaledCutClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, BlockDiagonalColoredSignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal colored signed classical")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // BlockDiagonalColoredSignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, BlockDiagonalNoColoredSignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  // this test is only compatible with rank higher than 1
  if (comm->getSize() == 1) {
    return;
  }

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561 = 3^8)
  // Nice size for 1D and perfect aggregation on small numbers of processors. (8748 = 4*3^7)
  Teuchos::CommandLineProcessor clp(false);
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);

  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  //    getCrsGraph()->getImporter()
  RCP<const Import> importer = ImportFactory::Build(A->getRowMap(), map);
  fineLevel.Set("Importer", importer);
  auto importerTest = A->getCrsGraph()->getImporter();  // NULL
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal colored signed classical")));
  coalesceDropFact.SetParameter("aggregation: coloring: localize color graph", Teuchos::ParameterEntry(false));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  // Need an importer
  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // BlockDiagonalNoColoredSignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, BlockDiagonalSignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal signed classical")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // BlockDiagonalSignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, BlockDiagonal, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(8.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, BlockDiagonalClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(8.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal classical")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, BlockDiagonalDistanceLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, BlockDiagonalDistanceDifferentCoordinatesLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  GO bnx = 15 * comm->getSize();
  Teuchos::ParameterList bMatrixList;
  matrixList.set("bnx", bnx);
  RCP<Matrix> B = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(bMatrixList, lib);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", B->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, BlockDiagonalDistanceLaplacianWeighted, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", (GO)10);
  matrixList.set("nz", (GO)10);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("3D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  std::vector<double> weights_v{100.0, 1.0, 1.0, 1.0, 100, 1.0, 1.0, 1.0, 100.0};
  Teuchos::Array<double> weights(weights_v);
  coalesceDropFact.SetParameter("aggregation: distance laplacian directional weights", Teuchos::ParameterEntry(weights));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, DistanceLaplacianWeighted, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", (GO)10);
  matrixList.set("nz", (GO)10);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixList, lib);

  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("3D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  std::vector<double> weights_v{100.0, 1.0, 1.0};
  Teuchos::Array<double> weights(weights_v);
  coalesceDropFact.SetParameter("aggregation: distance laplacian directional weights", Teuchos::ParameterEntry(weights));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, AggresiveDroppingIsMarkedAsBoundary, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // Test that when everything but the diagonal is dropped, the node is marked as boundary
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  RCP<const Map> dofMap    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 12 * comm->getSize(), 0, comm);
  Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);

  {
    Level fineLevel;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

    mtx->SetFixedBlockSize(1);
    fineLevel.Set("A", mtx);

    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetDefaultVerbLevel(MueLu::Extreme);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
    dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(4.1));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

    auto boundaryNodes         = graph->GetBoundaryNodeMap();
    bool allNodesAreOnBoundary = true;
    for (LO i = 0; i < Teuchos::as<LO>(boundaryNodes.size()); i++)
      allNodesAreOnBoundary &= boundaryNodes[i];
    TEST_EQUALITY(allNodesAreOnBoundary, true);
  }

  {
    Level fineLevel;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

    mtx->SetFixedBlockSize(2);
    fineLevel.Set("A", mtx);

    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetDefaultVerbLevel(MueLu::Extreme);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
    dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(4.1));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

    auto boundaryNodes         = graph->GetBoundaryNodeMap();
    bool allNodesAreOnBoundary = true;
    for (LO i = 0; i < Teuchos::as<LO>(boundaryNodes.size()); i++)
      allNodesAreOnBoundary &= boundaryNodes[i];
    TEST_EQUALITY(allNodesAreOnBoundary, true);
  }

  {
    Level fineLevel;
    TestHelpers::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

    mtx->SetFixedBlockSize(3);
    fineLevel.Set("A", mtx);

    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetDefaultVerbLevel(MueLu::Extreme);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("lightweight wrap", Teuchos::ParameterEntry(true));
    dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(4.1));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &dropFact);

    auto boundaryNodes         = graph->GetBoundaryNodeMap();
    bool allNodesAreOnBoundary = true;
    for (LO i = 0; i < Teuchos::as<LO>(boundaryNodes.size()); i++)
      allNodesAreOnBoundary &= boundaryNodes[i];
    TEST_EQUALITY(allNodesAreOnBoundary, true);
  }

}  // AggresiveDroppingIsMarkedAsBoundary

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory, SignedClassicalSA, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib           = TestHelpers::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", (GO)10);
  matrixList.set("nz", (GO)10);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixList, lib);

  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory coalesceDropFact;
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical sa")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph> graph = fineLevel.Get<RCP<LWGraph> >("Graph", &coalesceDropFact);
  LO myDofsPerNode   = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, Constructor, SC, LO, GO, Node)                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, Build, SC, LO, GO, Node)                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, PreDrop, SC, LO, GO, Node)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationBasic, SC, LO, GO, Node)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationStrided, SC, LO, GO, Node)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationStrided2, SC, LO, GO, Node)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationStridedOffset, SC, LO, GO, Node)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationLightweight, SC, LO, GO, Node)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationLightweightDrop, SC, LO, GO, Node)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationStridedLW, SC, LO, GO, Node)                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationStrided2LW, SC, LO, GO, Node)                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationStridedOffsetLW, SC, LO, GO, Node)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationDroppingLW, SC, LO, GO, Node)                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AmalgamationStridedOffsetDropping2LW, SC, LO, GO, Node)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, DistanceLaplacian, SC, LO, GO, Node)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, DistanceLaplacianScaledCut, SC, LO, GO, Node)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, DistanceLaplacianUnscaledCut, SC, LO, GO, Node)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, DistanceLaplacianCutSym, SC, LO, GO, Node)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, ClassicalScaledCut, SC, LO, GO, Node)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, ClassicalUnScaledCut, SC, LO, GO, Node)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, SignaledClassical, SC, LO, GO, Node)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, SignaledScaledCutClassical, SC, LO, GO, Node)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, SignaledUnscaledCutClassical, SC, LO, GO, Node)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, BlockDiagonalColoredSignedClassical, SC, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, BlockDiagonalNoColoredSignedClassical, SC, LO, GO, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, BlockDiagonalSignedClassical, SC, LO, GO, Node)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, BlockDiagonal, SC, LO, GO, Node)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, BlockDiagonalDistanceLaplacian, SC, LO, GO, Node)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, BlockDiagonalDistanceLaplacianWeighted, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, DistanceLaplacianWeighted, SC, LO, GO, Node)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, AggresiveDroppingIsMarkedAsBoundary, SC, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory, SignedClassicalSA, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
