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
#include <MueLu_Aggregates.hpp>
#include <MueLu_InterfaceAggregationFactory.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>

namespace MueLuTests {

// Helper function to convert from a node map to a dof map with striding information
template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
stridedDofMapHelper(
    const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& nodeMap,
    int dofsPerNode) {
  using GO = GlobalOrdinal;
  using LO = LocalOrdinal;
  using NO = Node;
  using Teuchos::ArrayView;
  using Teuchos::RCP;

  RCP<const Teuchos::Comm<int>> comm = nodeMap->getComm();
  const GO indexBase                 = 0;

  Teuchos::ArrayView<const GO> nodeGIDs = nodeMap->getLocalElementList();

  std::vector<GO> dofGIDs;
  dofGIDs.reserve(nodeGIDs.size() * dofsPerNode);
  for (GO nodeID : nodeGIDs)
    for (int d = 0; d < dofsPerNode; ++d)
      dofGIDs.push_back(nodeID * dofsPerNode + d);

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(dofsPerNode);

  size_t numGlobalDofs = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
  int stridedBlockId   = 0;
  GO offset            = 0;

  auto dofMap = Xpetra::StridedMapFactory<LO, GO, Node>::Build(
      nodeMap->lib(),
      numGlobalDofs,
      Teuchos::arrayViewFromVector(dofGIDs),
      indexBase,
      stridingInfo,
      nodeMap->getComm(),
      stridedBlockId,
      offset);

  return dofMap;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(InterfaceAggregationFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<InterfaceAggregationFactory> interfaceAggFact = rcp(new InterfaceAggregationFactory());
  TEST_EQUALITY(interfaceAggFact != Teuchos::null, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(InterfaceAggregationFactory, BuildBasedOnNodeMapping, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  using ST          = Teuchos::ScalarTraits<GO>;
  using device_type = typename Node::device_type;
  using LO_view     = Kokkos::View<LocalOrdinal*, device_type>;

  out << "version: " << MueLu::Version() << std::endl;

  Level level;
  // Set level ID to 0 so that we can provide the mapping to the InterfaceAggregationFactory
  level.SetLevelID(0);

  auto lib = TestHelpers::Parameters::getLib();

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  int rank                           = comm->getRank();
  int commSize                       = comm->getSize();

  ST::seedrandom(2654435761);

  // Problem definition for the construction of the A00 matrix later
  Teuchos::ParameterList matrixList;
  // We choose random problem dimensions of nx in [5,12] and ny in [14,32]
  GO nx, ny;
  if (rank == 0) {
    nx = 5 + ST::random() % 8;
    ny = 14 + ST::random() % 19;
  }
  Teuchos::broadcast(*comm, 0, &nx);
  Teuchos::broadcast(*comm, 0, &ny);
  matrixList.set("nx", nx);
  matrixList.set("ny", ny);
  // Distribution
  matrixList.set("mx", 1);
  matrixList.set("my", commSize);
  matrixList.set("matrixType", "Elasticity2D");

  // Number of dofs (degrees of freedom) per node for the primal field
  constexpr GO ndofnPrimal = 2;
  // Number of nodes for the primal field
  const GO nnodePrimal = nx * ny;
  LO mynnodePrimal;

  // We choose a fixed ndofn of 3 and a random nnode in [19,36] for the dual (interface) field
  constexpr GO ndofnDual = 3;
  GO nnodeDual;
  if (rank == 0) {
    nnodeDual = 19 + ST::random() % 18;
  }
  Teuchos::broadcast(*comm, 0, &nnodeDual);
  LO mynnodeDual;

  // Supply a random, injective global dual to primal mapping
  RCP<Array<GO>> dual2Primal = rcp(new Teuchos::Array<GO>());
  if (rank == 0) {
    RCP<Array<GO>> primalCandidates = rcp(new Teuchos::Array<GO>());
    for (GO i = 0; i < nnodePrimal; ++i)
      primalCandidates->push_back(i);
    for (GO dualNodeGID = 0; dualNodeGID < nnodeDual; ++dualNodeGID) {
      GO choose = ST::random() % (nnodePrimal - 1 - dualNodeGID);
      dual2Primal->push_back(primalCandidates->at(choose));
      primalCandidates->erase(primalCandidates->begin() + choose);
    }
  } else
    dual2Primal->resize(nnodeDual);
  Teuchos::broadcast(*comm, 0, nnodeDual, dual2Primal->getRawPtr());

  // Primal row map is uniform contiguous; nothing fancy
  RCP<const Map> nodeRowMapPrimal = Xpetra::MapFactory<LO, GO, NO>::createUniformContigMapWithNode(lib, nnodePrimal, comm);
  mynnodePrimal                   = nodeRowMapPrimal->getLocalNumElements();

  // We use the helper to make the corresponding strided row map
  RCP<const StridedMap> dofRowMapPrimal =
      Teuchos::rcp_dynamic_cast<const StridedMap>(
          stridedDofMapHelper<LO, GO, NO>(nodeRowMapPrimal, ndofnPrimal), true);

  // Distribute the dual2Primal mapping and dual row map
  Array<GO> myDualNodes;
  RCP<std::map<LO, LO>> myDual2Primal = rcp(new std::map<LO, LO>());
  mynnodeDual                         = 0;
  for (GO dualNodeGID = 0; dualNodeGID < nnodeDual; ++dualNodeGID) {
    GO primalGID = dual2Primal->at(dualNodeGID);
    if (nodeRowMapPrimal->isNodeGlobalElement(primalGID)) {
      (*myDual2Primal)[mynnodeDual] = nodeRowMapPrimal->getLocalElement(primalGID);
      myDualNodes.push_back(dualNodeGID);
      ++mynnodeDual;
    }
  }
  RCP<const Map> nodeRowMapDual = Xpetra::MapFactory<LO, GO, NO>::Build(
      lib, nnodeDual, myDualNodes, 0, comm);
  RCP<const Map> dofRowMapDual   = stridedDofMapHelper<LO, GO, NO>(nodeRowMapDual, ndofnDual);
  ArrayView<const GO> myDualDofs = dofRowMapDual->getLocalElementList();

  // Construct the upper left matrix A00
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
      Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixList.get("matrixType", ""), dofRowMapPrimal, matrixList);
  RCP<Matrix> A00 = Pr->BuildMatrix();
  A00->SetFixedBlockSize(ndofnPrimal);

  // Construct the lower right matrix A11 (empty, but with correct maps)
  RCP<Matrix> A11 = Teuchos::make_rcp<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(Xpetra::CrsMatrixFactory<SC, LO, GO, NO>::Build(dofRowMapDual, 0));
  A11->fillComplete(dofRowMapDual, dofRowMapDual);
  A11->SetFixedBlockSize(ndofnDual);

  // The true null space is not necessary for this test, so we construct a random one
  RCP<MultiVector> nullspace1 = MultiVectorFactory::Build(dofRowMapPrimal, ndofnPrimal);
  nullspace1->randomize();

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory> dropFact  = rcp(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Primal aggregation
  RCP<UncoupledAggregationFactory> uncoupledAggFact = rcp(new UncoupledAggregationFactory());
  uncoupledAggFact->SetFactory("Graph", dropFact);
  uncoupledAggFact->SetFactory("DofsPerNode", dropFact);
  level.Set("DofsPerNode", ndofnPrimal);
  level.Set("A", A00);
  level.Set("nullspace1", nullspace1);
  uncoupledAggFact->SetOrdering("natural");
  level.Request("Aggregates", uncoupledAggFact.get());
  uncoupledAggFact->Build(level);

  // Interface (dual) aggregation
  RCP<InterfaceAggregationFactory> interfaceAggFact = rcp(new InterfaceAggregationFactory());
  interfaceAggFact->SetFactory("Aggregates", uncoupledAggFact);
  interfaceAggFact->SetParameter("Dual/primal mapping strategy", Teuchos::ParameterEntry(std::string("node-based")));
  interfaceAggFact->SetParameter("number of DOFs per dual node", Teuchos::ParameterEntry(Teuchos::as<LO>(ndofnDual)));
  level.Set("DualNodeID2PrimalNodeID", myDual2Primal);
  level.Set("A", A11);
  interfaceAggFact->Build(level);

  // Access outputs and do some null checks
  level.Request("Aggregates", interfaceAggFact.get(), interfaceAggFact.get());
  level.Request("UnAmalgamationInfo", interfaceAggFact.get(), interfaceAggFact.get());
  RCP<Aggregates> primalAggs                   = level.Get<RCP<Aggregates>>("Aggregates", uncoupledAggFact.get());
  RCP<Aggregates> dualAggs                     = level.Get<RCP<Aggregates>>("Aggregates", interfaceAggFact.get());
  RCP<AmalgamationInfo> dualUnAmalgamationInfo = level.Get<RCP<AmalgamationInfo>>("UnAmalgamationInfo", interfaceAggFact.get());
  TEST_EQUALITY(dualAggs.is_null(), false);
  TEST_EQUALITY(dualUnAmalgamationInfo.is_null(), false);

  RCP<const Map> primalVertex2AggIdMap = primalAggs->GetVertex2AggId()->getMap();
  RCP<const Map> dualVertex2AggIdMap   = dualAggs->GetVertex2AggId()->getMap();

  const LO numPrimalAggs = primalAggs->GetNumAggregates();
  LO_view aggPtrPrimal;
  LO_view aggNodesPrimal;
  LO_view unaggregatedPrimal;
  primalAggs->ComputeNodesInAggregate(aggPtrPrimal, aggNodesPrimal, unaggregatedPrimal);

  const LO numDualAggs = dualAggs->GetNumAggregates();
  LO_view aggPtrDual;
  LO_view aggNodesDual;
  LO_view unaggregatedDual;
  dualAggs->ComputeNodesInAggregate(aggPtrDual, aggNodesDual, unaggregatedDual);

  // Write some outputs for a manual verification of the correct dual aggregates given the dual2Primal mapping
  std::stringstream myStream;
  myStream << std::endl;
  if (rank == 0) {
    myStream << "Global Dual to Primal mapping:" << std::endl;
    for (GO dualNodeId = 0; dualNodeId < nnodeDual; ++dualNodeId)
      myStream << "\t" << dualNodeId << " " << dual2Primal->at(dualNodeId) << std::endl;
    myStream << std::endl;
  }

  myStream << "++ pid " << rank << " ++ Aggregates" << std::endl;

  myStream << " - Primal Aggregates:" << std::endl;
  for (LO primalAggId = 0; primalAggId < numPrimalAggs; ++primalAggId) {
    myStream << "\tPrimal Aggregate with local AggID " << primalAggId << ":";
    for (LO ptrPos = aggPtrPrimal[primalAggId]; ptrPos < aggPtrPrimal[primalAggId + 1]; ++ptrPos) {
      // Global ID so we can compare with the global dual2Primal mapping
      const GO gNodeId = primalVertex2AggIdMap->getGlobalElement(aggNodesPrimal[ptrPos]);
      myStream << "  " << gNodeId;
    }
    myStream << std::endl;
  }

  myStream << " - Dual Aggregates:" << std::endl;
  for (LO dualAggId = 0; dualAggId < numDualAggs; ++dualAggId) {
    myStream << "\tDual Aggregate with local AggID " << dualAggId << ":";
    for (LO ptrPos = aggPtrDual[dualAggId]; ptrPos < aggPtrDual[dualAggId + 1]; ++ptrPos) {
      // Global ID so we can compare with the global dual2Primal mapping
      const GO gNodeId = dualVertex2AggIdMap->getGlobalElement(aggNodesDual[ptrPos]);
      myStream << "  " << gNodeId;
    }
    myStream << std::endl;
  }
  if (rank == 0)
    myStream << std::endl;

  for (int r = 0; r < commSize; ++r) {
    if (rank == r)
      out << myStream.str();
    comm->barrier();
  }

  // Note: even though GetVertex2AggId() includes ghost nodes, the locally owned node IDs always come first
  ArrayRCP<const LO> primalVertex2AggId = primalAggs->GetVertex2AggId()->getData(0);
  ArrayRCP<const LO> dualVertex2AggId   = dualAggs->GetVertex2AggId()->getData(0);

  // Check for unaggregated nodes (excluding ghost nodes)
  bool allAggregatedCheckPrimal = true;
  bool allAggregatedCheckDual   = true;
  for (LO i = 0; i < mynnodePrimal; ++i)
    if (primalVertex2AggId[i] == Teuchos::OrdinalTraits<LO>::invalid()) {
      allAggregatedCheckPrimal = false;
      break;
    }
  for (LO i = 0; i < mynnodeDual; ++i)
    if (dualVertex2AggId[i] == Teuchos::OrdinalTraits<LO>::invalid()) {
      allAggregatedCheckDual = false;
      break;
    }
  TEST_EQUALITY(allAggregatedCheckPrimal, true);
  TEST_EQUALITY(allAggregatedCheckDual, true);

  /* Formal check that the core idea is upheld
   * - dual aggregates should be the restriction of primal aggregates onto the interface
   * - hence, we check that the (partial) primal aggregates which are mapped back from the dual aggregates are a subset of the actual primal aggregates
   */
  bool aggStructureCheck = true;
  for (LO dualAggId = 0; dualAggId < numDualAggs && aggStructureCheck; ++dualAggId) {
    LO primalAggId_mapped_firstIter = -1;
    for (LO ptrPos = aggPtrDual[dualAggId]; ptrPos < aggPtrDual[dualAggId + 1]; ++ptrPos) {
      LO dualNodeId = aggNodesDual[ptrPos];

      // Sanity check for complete dual2Primal mapping
      if (myDual2Primal->find(dualNodeId) == myDual2Primal->end()) {
        aggStructureCheck = false;
        break;
      }
      LO primalNodeId_mapped = myDual2Primal->at(dualNodeId);

      // Sanity check for local primal node ID
      if (primalNodeId_mapped >= mynnodePrimal) {
        aggStructureCheck = false;
        break;
      }
      LO primalAggId_mapped = primalVertex2AggId[primalNodeId_mapped];

      // Check that all of the corresponding primal nodes are also in the same aggregate
      if (primalAggId_mapped_firstIter == -1)
        primalAggId_mapped_firstIter = primalAggId_mapped;
      else if (primalAggId_mapped != primalAggId_mapped_firstIter) {
        aggStructureCheck = false;
        break;
      }
    }
  }
  TEST_EQUALITY(aggStructureCheck, true);
}

#define MUELU_ETI_GROUP(SC, LO, GO, Node)                                                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(InterfaceAggregationFactory, Constructor, SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(InterfaceAggregationFactory, BuildBasedOnNodeMapping, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
