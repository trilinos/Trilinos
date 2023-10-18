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
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// Test
#include "contact_Helpers.hpp"

// MueLu
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_BlockedCoarseMapFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_CreateXpetraPreconditioner.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_InterfaceAggregationFactory.hpp"
#include "MueLu_MapTransferFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_SegregatedAFactory.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"

// Teuchos
#include <Teuchos_XMLParameterListHelpers.hpp>

// Xpetra
#include <Xpetra_BlockedMap.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using ST = Teuchos::ScalarTraits<Scalar>;
  using GST = Teuchos::ScalarTraits<GO>;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  std::string probName = ""; clp.setOption("probName", &probName, "Short name of the problem. Used to read-in the problem from files.");
  GO numGlobalDofPrimal = -GST::one(); clp.setOption("nPrimalDof", &numGlobalDofPrimal, "total number of primal Dof");
  GO numGlobalDofDual = -GST::one(); clp.setOption("nDualDof", &numGlobalDofDual, "total number of dual Dof");
  int numPrimalDofPerNode = -1; clp.setOption("numPrimalDofPerNode", &numPrimalDofPerNode, "number of primal Dof per mesh node");
  int numDualDofPerNode = -1; clp.setOption("numDualDofPerNode", &numDualDofPerNode, "number of dual Dof per interface node");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(probName=="", MueLu::Exceptions::InvalidArgument, "Please specify a valid problem name.");
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalDofPrimal==-GST::one(), MueLu::Exceptions::InvalidArgument, "Please specify the global number of primal Dof on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalDofDual==-GST::one(), MueLu::Exceptions::InvalidArgument, "Please specify the global number of dual Dof on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numPrimalDofPerNode==-1, MueLu::Exceptions::InvalidArgument, "Please specify the number of primal Dof per mesh node on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numDualDofPerNode==-1, MueLu::Exceptions::InvalidArgument, "Please specify the number of dual Dof per interface node on the command line.");

  const std::string matrixFileName = probName + "_block_matrix.mm";
  const std::string nullspace1FileName = probName + "_nullspace.mm";
  const std::string map1Name = "dropMap1";
  const std::string map2Name = "dropMap2";
  const std::string droppingScheme = "map-pair";

  //// Rowmaps
  // Primal dofmap
  std::vector<size_t> stridingInfoPrimal;
  stridingInfoPrimal.push_back(numPrimalDofPerNode);
  RCP<const StridedMap> StridedDofRowMapPrimal = StridedMapFactory::Build(lib, numGlobalDofPrimal, Teuchos::ScalarTraits<GO>::zero(), stridingInfoPrimal,comm,-1);

  // Dual dofmap
  std::vector<size_t> stridingInfoDual;
  stridingInfoDual.push_back(numDualDofPerNode);
  RCP<const StridedMap> StridedDofRowMapDual = StridedMapFactory::Build(lib, numGlobalDofDual, Teuchos::ScalarTraits<GO>::zero(), stridingInfoDual, comm, -1, numGlobalDofPrimal);

  // Construct the blocked map of the global system
  std::vector<RCP<const Map> > rowmaps;
  rowmaps.push_back(StridedDofRowMapPrimal);
  rowmaps.push_back(StridedDofRowMapDual);
  RCP<const Map> fullRowMap = MapUtils::concatenateMaps(rowmaps);
  RCP<const BlockedMap> blockedMap = rcp(new BlockedMap(fullRowMap, rowmaps));

  //// Matrix A
  // Read the matrix from file and transform it into a block matrix
  RCP<Matrix> mat = Xpetra::IO<SC,LO,GO,NO>::Read(matrixFileName, fullRowMap);
  RCP<MapExtractor> rangeMapExtractor = Xpetra::MapExtractorFactory<SC,LO,GO,NO>::Build(fullRowMap, rowmaps);
  RCP<BlockedCrsMatrix> blockedMatrix = Xpetra::MatrixUtils<SC,LO,GO,NO>::SplitMatrix(*mat,rangeMapExtractor, rangeMapExtractor);

  blockedMatrix->fillComplete();

  //// Nullspaces
  // Read the nullspace vector of the (0,0)-block from file
  RCP<MultiVector> nullspace1 = Xpetra::IO<SC,LO,GO,NO>::ReadMultiVector(nullspace1FileName, StridedDofRowMapPrimal);

  // Create the default nullspace vector of the (1,1)-block
  RCP<MultiVector> nullspace2 = MultiVectorFactory::Build(StridedDofRowMapDual, 2, true);
  const int dimNS = 2;
  for (int dim = 0; dim < dimNS; ++dim)
  {
    ArrayRCP<Scalar> nsValues = nullspace2->getDataNonConst(dim);
    const int numBlocks = nsValues.size() / dimNS;
    for (int j = 0; j < numBlocks; ++j)
      nsValues[j * dimNS + dim] = Teuchos::ScalarTraits<Scalar>::one();
  }

  // Create the interface slave and master dof row map (= map-pair for matrix filtering in segregatedAFactory
  std::vector<GO> slavemapEntries{42, 43, 44, 45, 52, 53};
  RCP<const Map> slavemap = Xpetra::MapFactory<LO, GO, NO>::Build(lib, slavemapEntries.size(), slavemapEntries, 0, comm);
  std::vector<GO> mastermapEntries{0, 1, 6, 7, 22, 23};
  RCP<const Map> mastermap = Xpetra::MapFactory<LO, GO, NO>::Build(lib, mastermapEntries.size(), mastermapEntries, 0, comm);

  // Distribute interface slave dof row map across all procs
  RCP<const Map> primalInterfaceDofMap = Teuchos::null;
  {
    std::vector<int> interfaceGlobalDofOnCurrentProc = std::vector<int>(slavemap->getGlobalNumElements());
    std::vector<int> interfaceGlobalDof = std::vector<int>(slavemap->getGlobalNumElements());
    if(comm->getRank() == 0)
    {
      for(size_t i = 0; i < slavemap->getGlobalNumElements(); ++i)
      {
        interfaceGlobalDofOnCurrentProc[i] = slavemap->getLocalElementList()[i];
      }
    }
    Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_MAX, interfaceGlobalDof.size(),
                            &interfaceGlobalDofOnCurrentProc[0], &interfaceGlobalDof[0]);

    Array<GlobalOrdinal> primalInterfaceDofMapDofOnCurProc;
    for(size_t i = 0; i < interfaceGlobalDof.size(); ++i)
    {
      if(StridedDofRowMapPrimal->isNodeGlobalElement((GlobalOrdinal) interfaceGlobalDof[i]))
        primalInterfaceDofMapDofOnCurProc.push_back((GlobalOrdinal) interfaceGlobalDof[i]);
    }
    primalInterfaceDofMap = MapFactory::Build(lib, slavemap->getGlobalNumElements(), primalInterfaceDofMapDofOnCurProc, Teuchos::ScalarTraits<GO>::zero(), comm);
  }

  ///////////// MueLu Factories Level0 /////////////
  //// Block (0,0)
  // define SubBlockAFactory for block (0,0) that provides Matrix A to SegregatedAFactory
  RCP<SubBlockAFactory> subBlockAFact00 = RCP(new SubBlockAFactory());
  subBlockAFact00->SetFactory("A", MueLu::NoFactory::getRCP());
  subBlockAFact00->SetParameter("block row", Teuchos::ParameterEntry(0));
  subBlockAFact00->SetParameter("block col", Teuchos::ParameterEntry(0));

  // define segregatedAFactory that provides Matrix A to CoalesceDropFactory
  RCP<SegregatedAFactory> segregatedAFact = RCP(new SegregatedAFactory());
  segregatedAFact->SetFactory("A", subBlockAFact00);
  segregatedAFact->SetParameter("droppingScheme", Teuchos::ParameterEntry(droppingScheme));
  segregatedAFact->SetParameter("Call ReduceAll on dropMap1", Teuchos::ParameterEntry(true));
  segregatedAFact->SetParameter("Call ReduceAll on dropMap2", Teuchos::ParameterEntry(true));

  // define coarse Map factory
  RCP<CoarseMapFactory> coarseMapFact00 = RCP(new CoarseMapFactory());

  // define amalgamation factory
  RCP<AmalgamationFactory> amalgFact00 = RCP(new AmalgamationFactory());
  amalgFact00->SetFactory("A", subBlockAFact00);

  // define CoalesceDropFactory that provides Graph
  RCP<CoalesceDropFactory> dropFact00 = RCP(new CoalesceDropFactory());
  dropFact00->SetFactory("UnAmalgamationInfo", amalgFact00);
  dropFact00->SetFactory("A", segregatedAFact);

  // Setup aggregation factory
  RCP<UncoupledAggregationFactory> uncoupledAggFact = rcp(new UncoupledAggregationFactory());
  uncoupledAggFact->SetFactory("Graph", dropFact00);
  uncoupledAggFact->SetFactory("DofsPerNode", dropFact00);
  uncoupledAggFact->SetOrdering("graph");
  uncoupledAggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(2));
  uncoupledAggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));

  // Setup TentativePFactory
  RCP<TentativePFactory> tentativePFact00 = rcp(new TentativePFactory());
  tentativePFact00->SetFactory("A", segregatedAFact);
  tentativePFact00->SetFactory("Aggregates", uncoupledAggFact);
  tentativePFact00->SetFactory("CoarseMap", coarseMapFact00);
  tentativePFact00->SetParameter("tentative: calculate qr", Teuchos::ParameterEntry(true));
  tentativePFact00->SetParameter("tentative: build coarse coordinates", Teuchos::ParameterEntry(false));

  // Setup nullspace Factory
  RCP<NullspaceFactory> nullspace1Fact = rcp(new NullspaceFactory());
  nullspace1Fact->SetParameter("Fine level nullspace", Teuchos::ParameterEntry(std::string("Nullspace1")));
  nullspace1Fact->SetFactory("Nullspace", tentativePFact00);

  // Setup Map Transfer Factory for map1
  RCP<MapTransferFactory> slaveMapTransferFact = Teuchos::rcp(new MapTransferFactory());
  slaveMapTransferFact->SetParameter("map: name", Teuchos::ParameterEntry(map1Name));
  slaveMapTransferFact->SetParameter("map: factory", Teuchos::ParameterEntry(map1Name));
  slaveMapTransferFact->SetFactory("P", tentativePFact00);

  RCP<MapTransferFactory> masterMapTransferFact = Teuchos::rcp(new MapTransferFactory());
  masterMapTransferFact->SetParameter("map: name", Teuchos::ParameterEntry(map2Name));
  masterMapTransferFact->SetParameter("map: factory", Teuchos::ParameterEntry(map2Name));
  masterMapTransferFact->SetFactory("P", tentativePFact00);

  segregatedAFact->SetFactory(map1Name, slaveMapTransferFact);
  segregatedAFact->SetFactory(map2Name, masterMapTransferFact);

  //// Block (1,1)
  // SubBlockAFact for block (1,1)
  RCP<SubBlockAFactory> subBlockAFact11 = rcp(new SubBlockAFactory());
  subBlockAFact11->SetFactory("A", MueLu::NoFactory::getRCP());
  subBlockAFact11->SetParameter("block row", Teuchos::ParameterEntry(1));
  subBlockAFact11->SetParameter("block col", Teuchos::ParameterEntry(1));

  //// Block (0,1) (off-diagonal) for InterfaceAggregationFactory
  // SubBlockAFact for block (0,1)
  RCP<SubBlockAFactory> subBlockAFact01 = rcp(new SubBlockAFactory);
  subBlockAFact01->SetFactory("A", MueLu::NoFactory::getRCP());
  subBlockAFact01->SetParameter("block row", Teuchos::ParameterEntry(0));
  subBlockAFact01->SetParameter("block col", Teuchos::ParameterEntry(1));

  // define BlockedCoarseMapFactory
  RCP<BlockedCoarseMapFactory> coarseMapFact01 = rcp(new BlockedCoarseMapFactory());
  coarseMapFact01->SetFactory("CoarseMap", coarseMapFact00);

  // define MapTransferFactory for "Primal interface DOF map"
  RCP<MapTransferFactory> primalInterfaceDofMapTransferFact = rcp (new MapTransferFactory());
  primalInterfaceDofMapTransferFact->SetParameter("map: factory", Teuchos::ParameterEntry(std::string("Primal interface DOF map")));
  primalInterfaceDofMapTransferFact->SetParameter("map: name", Teuchos::ParameterEntry(std::string("Primal interface DOF map")));
  primalInterfaceDofMapTransferFact->SetFactory("P", tentativePFact00);
  primalInterfaceDofMapTransferFact->SetParameter("nullspace vectors: limit to", Teuchos::ParameterEntry(std::string("translations")));

  // define factory for interface aggregates
  RCP<InterfaceAggregationFactory> interfaceAggFact = rcp(new InterfaceAggregationFactory());
  interfaceAggFact->SetFactory("A", subBlockAFact01);
  interfaceAggFact->SetFactory("Aggregates", uncoupledAggFact);
  interfaceAggFact->SetParameter("Dual/primal mapping strategy", Teuchos::ParameterEntry(std::string("dof-based")));
  interfaceAggFact->SetFactory("Primal interface DOF map", primalInterfaceDofMapTransferFact);

  // define tentativeP Factory
  RCP<TentativePFactory> tentativePFact01 = rcp(new TentativePFactory());
  tentativePFact01->SetFactory("A", subBlockAFact11);
  tentativePFact01->SetFactory("Aggregates", interfaceAggFact);
  tentativePFact01->SetFactory("CoarseMap", coarseMapFact01);
  tentativePFact01->SetFactory("UnAmalgamationInfo", interfaceAggFact);
  tentativePFact01->SetParameter("tentative: calculate qr", Teuchos::ParameterEntry(true));
  tentativePFact01->SetParameter("tentative: build coarse coordinates", Teuchos::ParameterEntry(false));

  // define Factory for nullspace2
  RCP<NullspaceFactory> nullspace2Fact = rcp(new NullspaceFactory());
  nullspace2Fact->SetParameter("Fine level nullspace", Teuchos::ParameterEntry(std::string("Nullspace2")));
  nullspace2Fact->SetFactory("Nullspace", tentativePFact01);

  // define CoarseMapFact01 dependencies
  coarseMapFact01->SetFactory("Aggregates", interfaceAggFact);
  coarseMapFact01->SetFactory("CoarseMap", coarseMapFact00);
  coarseMapFact01->SetFactory("Nullspace", nullspace2Fact);

  //// Factory Managers
  // First group
  RCP<FactoryManager> M1 = rcp(new FactoryManager());
  M1->SetFactory("A", subBlockAFact00);
  M1->SetFactory("P", tentativePFact00);
  M1->SetFactory("Aggregates", uncoupledAggFact);
  M1->SetFactory("Nullspace", nullspace1Fact);
  M1->SetFactory("CoarseMap", coarseMapFact00);
  M1->SetFactory("UnAmalgamationInfo", amalgFact00);
    M1->SetFactory(map1Name, slaveMapTransferFact);
    M1->SetFactory(map2Name, masterMapTransferFact);
  M1->SetKokkosRefactor(false);

  // Second group
  RCP<FactoryManager> M2 = rcp(new FactoryManager());
  M2->SetFactory("A", subBlockAFact11);
  M2->SetFactory("P", tentativePFact01);
  M2->SetFactory("Aggregates", interfaceAggFact);
  M2->SetFactory("Nullspace", nullspace2Fact);
  M2->SetFactory("CoarseMap", coarseMapFact01);
  M2->SetFactory("Primal interface DOF map", primalInterfaceDofMapTransferFact);
  M2->SetKokkosRefactor(false);

  //// Blocked Transfer Operators
  RCP<BlockedPFactory> blockedPFact = rcp(new BlockedPFactory);
  blockedPFact->AddFactoryManager(M1);
  blockedPFact->AddFactoryManager(M2);

  RCP<GenericRFactory> blockedRFact = rcp(new GenericRFactory());
  blockedRFact->SetFactory("P", blockedPFact);

  RCP<BlockedRAPFactory> blockedRAPFact = rcp(new BlockedRAPFactory());
  blockedRAPFact->SetFactory("P", blockedPFact);
  blockedRAPFact->SetFactory("R", blockedRFact);

  //// Main Factory Manager
  RCP<FactoryManager> manager = rcp(new FactoryManager());

  manager->SetFactory("A", blockedRAPFact);
  manager->SetFactory("P", blockedPFact);
  manager->SetFactory("R", blockedRFact);
  manager->SetFactory(map1Name, MueLu::NoFactory::getRCP());
  manager->SetFactory(map2Name, MueLu::NoFactory::getRCP());
  manager->SetKokkosRefactor(false);

  ///////////// Level Definition /////////////
  RCP<Level> levelZero = rcp(new Level());
  RCP<Level> levelOne = rcp(new Level());
  levelZero->SetLevelID(0);
  levelOne->SetLevelID(1);
  levelOne->SetPreviousLevel(levelZero);
  TEUCHOS_ASSERT_EQUALITY(levelZero->GetLevelID(), 0);
  TEUCHOS_ASSERT_EQUALITY(levelOne->GetLevelID(), 1);

  levelZero->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(blockedMatrix));
  levelZero->Set("Primal interface DOF map", primalInterfaceDofMap);
  levelZero->Set(map1Name, slavemap);
  levelZero->Set(map2Name, mastermap);
  levelZero->Set("Nullspace1", nullspace1);
  levelZero->Set("Nullspace2", nullspace2);

  TEUCHOS_ASSERT(levelZero->IsAvailable("A", MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable(map1Name, MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable(map2Name, MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable("Nullspace1", MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable("Nullspace2", MueLu::NoFactory::get()));

  levelZero->SetFactoryManager(manager);
  levelOne->SetFactoryManager(manager);

  ///////////// Get Level Information /////////////
  //// Level 0
  levelZero->Request("Aggregates", uncoupledAggFact.get());
  levelZero->Request(*uncoupledAggFact);
  TEUCHOS_ASSERT(levelZero->IsRequested("Aggregates", uncoupledAggFact.get()));
  uncoupledAggFact->Build(*levelZero);
  TEUCHOS_ASSERT(levelZero->IsAvailable("Aggregates", uncoupledAggFact.get()));
  RCP<Aggregates> primalAggsLevelZero = levelZero->Get<RCP<Aggregates>>("Aggregates", uncoupledAggFact.get());

  MueLuTests::checkAggregatesMapPair<SC, LO, GO, NO>(primalAggsLevelZero, stridingInfoPrimal,
                                                     levelZero->Get<RCP<const Map>>(map1Name, MueLu::NoFactory::get()),
                                                       levelZero->Get<RCP<const Map>>(map2Name, MueLu::NoFactory::get()));

  //// Level 1
  levelOne->Request("A", blockedRAPFact.get());
  levelOne->Request("P", tentativePFact00.get());
  TEUCHOS_ASSERT(levelOne->IsRequested("A", blockedRAPFact.get()));
  TEUCHOS_ASSERT(levelOne->IsRequested("P", tentativePFact00.get()));

  blockedRAPFact->Build(*levelZero, *levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable("A", blockedRAPFact.get()));
  TEUCHOS_ASSERT(levelOne->IsAvailable("P", tentativePFact00.get()));

  levelOne->Request(map1Name, MueLu::NoFactory::get());
  levelOne->Request(map2Name, MueLu::NoFactory::get());
  TEUCHOS_ASSERT(levelOne->IsRequested(map1Name, MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelOne->IsRequested(map2Name, MueLu::NoFactory::get()));

  levelOne->Request("P", tentativePFact00.get(), slaveMapTransferFact.get());
  levelOne->Request(*slaveMapTransferFact);

  TEUCHOS_ASSERT(levelOne->IsRequested("P", tentativePFact00.get()));
  RCP<Matrix> Ptent = levelOne->Get<RCP<Matrix>>("P", tentativePFact00.get());
  TEUCHOS_ASSERT(!Ptent.is_null());

  slaveMapTransferFact->Build(*levelZero, *levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable(map1Name, MueLu::NoFactory::get()));
  levelOne->Request(map2Name, MueLu::NoFactory::get());
  levelOne->Request("P", tentativePFact00.get(), masterMapTransferFact.get());
  levelOne->Request(*masterMapTransferFact);
  TEUCHOS_ASSERT(levelOne->IsRequested(map2Name, MueLu::NoFactory::get()));
  masterMapTransferFact->Build(*levelZero, *levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable(map2Name, MueLu::NoFactory::get()));

  // SubBLockAFactory -> A00
  RCP<SubBlockAFactory> subBlockAFact00LevelOne = rcp(new SubBlockAFactory());
  subBlockAFact00LevelOne->SetFactory("A", blockedRAPFact);
  subBlockAFact00LevelOne->SetParameter("block row", Teuchos::ParameterEntry(0));
  subBlockAFact00LevelOne->SetParameter("block col", Teuchos::ParameterEntry(0));

  RCP<SegregatedAFactory> segregatedAFactLevelOne = rcp(new SegregatedAFactory());
  segregatedAFactLevelOne->SetFactory("A", subBlockAFact00LevelOne);
  segregatedAFactLevelOne->SetParameter("droppingScheme", Teuchos::ParameterEntry(droppingScheme));
  segregatedAFact->SetParameter("Call ReduceAll on dropMap1", Teuchos::ParameterEntry(true));
  segregatedAFact->SetParameter("Call ReduceAll on dropMap2", Teuchos::ParameterEntry(true));

  // define amalgamation factory
  RCP<AmalgamationFactory> amalgFact00LevelOne = RCP(new AmalgamationFactory());
  amalgFact00LevelOne->SetFactory("A", subBlockAFact00LevelOne);

  // define CoalesceDropFactory that provides Graph
  RCP<CoalesceDropFactory> dropFact00LevelOne = RCP(new CoalesceDropFactory());
  dropFact00LevelOne->SetFactory("UnAmalgamationInfo", amalgFact00LevelOne);
  dropFact00LevelOne->SetFactory("A", segregatedAFactLevelOne);

  // Setup aggregation factory
  RCP<UncoupledAggregationFactory> uncoupledAggFactLevelOne = rcp(new UncoupledAggregationFactory());
  uncoupledAggFactLevelOne->SetFactory("Graph", dropFact00LevelOne);
  uncoupledAggFactLevelOne->SetFactory("DofsPerNode", dropFact00LevelOne);
  uncoupledAggFactLevelOne->SetOrdering("graph");
  uncoupledAggFactLevelOne->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(2));
  uncoupledAggFactLevelOne->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));

  // Request level info
  levelOne->Request("A", subBlockAFact00LevelOne.get());
  levelOne->Request(*subBlockAFact00LevelOne);
  subBlockAFact00LevelOne->Build(*levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable("A", subBlockAFact00LevelOne.get()));

  // SegregatedAFactory -> A00 filtered
  levelOne->Request("A", segregatedAFactLevelOne.get());
  levelOne->Request(*segregatedAFactLevelOne);
  segregatedAFactLevelOne->Build(*levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable("A", segregatedAFactLevelOne.get()));

  // AmalgamationFactory of block A00
  levelOne->Request("UnAmalgamationInfo", amalgFact00LevelOne.get());
  levelOne->Request(*amalgFact00LevelOne);
  amalgFact00LevelOne->Build(*levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable("UnAmalgamationInfo", amalgFact00LevelOne.get()));

  levelOne->Request("Aggregates", uncoupledAggFactLevelOne.get());
  TEUCHOS_ASSERT(levelOne->IsRequested("Aggregates", uncoupledAggFactLevelOne.get()));
  uncoupledAggFactLevelOne->Build(*levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable("Aggregates", uncoupledAggFactLevelOne.get()));

  RCP<Aggregates> primalAggsLevelOne = levelOne->Get<RCP<Aggregates>>("Aggregates", uncoupledAggFactLevelOne.get());

  MueLuTests::checkAggregatesMapPair<SC, LO, GO, NO>(primalAggsLevelOne, stridingInfoPrimal,
                                                     levelOne->Get<RCP<const Map>>(map1Name, MueLu::NoFactory::get()),
                                                     levelOne->Get<RCP<const Map>>(map2Name, MueLu::NoFactory::get()));
  return EXIT_SUCCESS;
}

//-----------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}