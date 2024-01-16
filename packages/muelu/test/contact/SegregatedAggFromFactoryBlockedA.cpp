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
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using ST  = Teuchos::ScalarTraits<Scalar>;
  using GST = Teuchos::ScalarTraits<GO>;

  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
  auto numRanks                      = comm->getSize();

  std::string probName = "";
  clp.setOption("probName", &probName, "Short name of the problem. Used to read-in the problem from files.");
  GO numGlobalDofPrimal = -GST::one();
  clp.setOption("nPrimalDof", &numGlobalDofPrimal, "total number of primal Dof");
  GO numGlobalDofDual = -GST::one();
  clp.setOption("nDualDof", &numGlobalDofDual, "total number of dual Dof");
  int numPrimalDofPerNode = -1;
  clp.setOption("numPrimalDofPerNode", &numPrimalDofPerNode, "number of primal Dof per mesh node");
  int numDualDofPerNode = -1;
  clp.setOption("numDualDofPerNode", &numDualDofPerNode, "number of dual Dof per interface node");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
      break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(probName == "", MueLu::Exceptions::InvalidArgument,
                             "Please specify a valid problem name.");
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalDofPrimal == -GST::one(), MueLu::Exceptions::InvalidArgument,
                             "Please specify the global number of primal Dof on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalDofDual == -GST::one(), MueLu::Exceptions::InvalidArgument,
                             "Please specify the global number of dual Dof on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numPrimalDofPerNode == -1, MueLu::Exceptions::InvalidArgument,
                             "Please specify the number of primal Dof per mesh node on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numDualDofPerNode == -1, MueLu::Exceptions::InvalidArgument,
                             "Please specify the number of dual Dof per interface node on the command line.");

  const std::string matrixFileName             = probName + "_block_matrix.mm";
  const std::string nullspace1FileName         = probName + "_nullspace.mm";
  const std::string primalInterfaceMapFileName = probName + "_interface_dof_map_MPI" + std::to_string(numRanks) + ".mm";
  const std::string dropMap1Name               = "dropMap1";
  const std::string dropMap2Name               = "dropMap2";
  const std::string droppingScheme             = "map-pair";

  /// Rowmaps
  // Primal dofmap
  std::vector<size_t> stridingInfoPrimal;
  stridingInfoPrimal.push_back(numPrimalDofPerNode);
  RCP<const StridedMap> StridedDofRowMapPrimal = StridedMapFactory::Build(lib, numGlobalDofPrimal,
                                                                          Teuchos::ScalarTraits<GO>::zero(),
                                                                          stridingInfoPrimal, comm, -1);

  // Dual dofmap
  std::vector<size_t> stridingInfoDual;
  stridingInfoDual.push_back(numDualDofPerNode);
  RCP<const StridedMap> StridedDofRowMapDual = StridedMapFactory::Build(lib, numGlobalDofDual,
                                                                        Teuchos::ScalarTraits<GO>::zero(),
                                                                        stridingInfoDual, comm, -1, numGlobalDofPrimal);

  // Construct the blocked map of the global system
  std::vector<RCP<const Map>> rowmaps;
  rowmaps.push_back(StridedDofRowMapPrimal);
  rowmaps.push_back(StridedDofRowMapDual);
  RCP<const Map> fullRowMap        = MapUtils::concatenateMaps(rowmaps);
  RCP<const BlockedMap> blockedMap = rcp(new BlockedMap(fullRowMap, rowmaps));

  /// Matrix A
  // Read the matrix from file and transform it into a block matrix
  RCP<Matrix> mat                     = Xpetra::IO<SC, LO, GO, NO>::Read(matrixFileName, fullRowMap);
  RCP<MapExtractor> rangeMapExtractor = Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullRowMap, rowmaps);
  RCP<BlockedCrsMatrix> blockedMatrix = Xpetra::MatrixUtils<SC, LO, GO, NO>::SplitMatrix(*mat, rangeMapExtractor,
                                                                                         rangeMapExtractor);

  blockedMatrix->fillComplete();

  /// Nullspaces
  // Read the nullspace vector of the (0,0)-block from file
  RCP<MultiVector> nullspace1 = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(nullspace1FileName, StridedDofRowMapPrimal);

  // Create the default nullspace vector of the (1,1)-block
  RCP<MultiVector> nullspace2 = MultiVectorFactory::Build(StridedDofRowMapDual, 2, true);
  const int dimNS             = 2;
  for (int dim = 0; dim < dimNS; ++dim) {
    ArrayRCP<Scalar> nsValues = nullspace2->getDataNonConst(dim);
    const int numBlocks       = nsValues.size() / dimNS;
    for (int j = 0; j < numBlocks; ++j)
      nsValues[j * dimNS + dim] = Teuchos::ScalarTraits<Scalar>::one();
  }

  // Read the interface slave dof row map from file
  TEUCHOS_ASSERT(!primalInterfaceMapFileName.empty());
  RCP<const Map> primalInterfaceDofMap = Xpetra::IO<SC, LO, GO, NO>::ReadMap(primalInterfaceMapFileName, lib, comm);
  // Create the interface master dof row map (for map-pair for matrix filtering in segregatedAFactory)
  std::vector<GO> mastermapEntries{0, 1, 6, 7, 22, 23};
  RCP<const Map> mastermap = Xpetra::MapFactory<LO, GO, NO>::Build(lib, mastermapEntries.size(), mastermapEntries, 0,
                                                                   comm);

  ///////////// Factories Definitions /////////////
  /// Factories: level 0
  /// Block (0,0)
  // define SubBlockAFactory for block (0,0) that provides Matrix A to SegregatedAFactory
  RCP<SubBlockAFactory> subBlockAFact00LevelZero = rcp(new SubBlockAFactory());
  subBlockAFact00LevelZero->SetFactory("A", MueLu::NoFactory::getRCP());
  subBlockAFact00LevelZero->SetParameter("block row", Teuchos::ParameterEntry(0));
  subBlockAFact00LevelZero->SetParameter("block col", Teuchos::ParameterEntry(0));

  // define segregatedAFactory that provides Matrix A to CoalesceDropFactory
  RCP<SegregatedAFactory> segregatedAFactLevelZero = rcp(new SegregatedAFactory());
  segregatedAFactLevelZero->SetFactory("A", subBlockAFact00LevelZero);
  segregatedAFactLevelZero->SetParameter("droppingScheme", Teuchos::ParameterEntry(droppingScheme));
  segregatedAFactLevelZero->SetParameter("Call ReduceAll on dropMap1", Teuchos::ParameterEntry(true));
  segregatedAFactLevelZero->SetParameter("Call ReduceAll on dropMap2", Teuchos::ParameterEntry(true));

  // define coarse Map factory
  RCP<CoarseMapFactory> coarseMapFact00LevelZero = rcp(new CoarseMapFactory());

  // define amalgamation factory
  RCP<AmalgamationFactory> amalgFact00LevelZero = rcp(new AmalgamationFactory());
  amalgFact00LevelZero->SetFactory("A", subBlockAFact00LevelZero);

  // define CoalesceDropFactory that provides Graph
  RCP<CoalesceDropFactory> dropFact00LevelZero = rcp(new CoalesceDropFactory());
  dropFact00LevelZero->SetFactory("UnAmalgamationInfo", amalgFact00LevelZero);
  dropFact00LevelZero->SetFactory("A", segregatedAFactLevelZero);

  // Setup aggregation factory
  RCP<UncoupledAggregationFactory> uncoupledAggFactLevelZero = rcp(new UncoupledAggregationFactory());
  uncoupledAggFactLevelZero->SetFactory("Graph", dropFact00LevelZero);
  uncoupledAggFactLevelZero->SetFactory("DofsPerNode", dropFact00LevelZero);
  uncoupledAggFactLevelZero->SetOrdering("graph");
  uncoupledAggFactLevelZero->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(2));
  uncoupledAggFactLevelZero->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));

  // Setup TentativePFactory
  RCP<TentativePFactory> tentativePFact00LevelZero = rcp(new TentativePFactory());
  tentativePFact00LevelZero->SetFactory("A", segregatedAFactLevelZero);
  tentativePFact00LevelZero->SetFactory("Aggregates", uncoupledAggFactLevelZero);
  tentativePFact00LevelZero->SetFactory("CoarseMap", coarseMapFact00LevelZero);
  tentativePFact00LevelZero->SetParameter("tentative: calculate qr", Teuchos::ParameterEntry(true));
  tentativePFact00LevelZero->SetParameter("tentative: build coarse coordinates", Teuchos::ParameterEntry(false));

  // Setup nullspace Factory
  RCP<NullspaceFactory> nullspace1FactLevelZero = rcp(new NullspaceFactory());
  nullspace1FactLevelZero->SetParameter("Fine level nullspace", Teuchos::ParameterEntry(std::string("Nullspace1")));
  nullspace1FactLevelZero->SetFactory("Nullspace1", tentativePFact00LevelZero);

  // Setup Map Transfer Factory for map1
  RCP<MapTransferFactory> slaveMapTransferFactLevelZero = Teuchos::rcp(new MapTransferFactory());
  slaveMapTransferFactLevelZero->SetParameter("map: name", Teuchos::ParameterEntry(dropMap1Name));
  slaveMapTransferFactLevelZero->SetParameter("map: factory", Teuchos::ParameterEntry(dropMap1Name));
  slaveMapTransferFactLevelZero->SetFactory("P", tentativePFact00LevelZero);

  RCP<MapTransferFactory> masterMapTransferFactLevelZero = Teuchos::rcp(new MapTransferFactory());
  masterMapTransferFactLevelZero->SetParameter("map: name", Teuchos::ParameterEntry(dropMap2Name));
  masterMapTransferFactLevelZero->SetParameter("map: factory", Teuchos::ParameterEntry(dropMap2Name));
  masterMapTransferFactLevelZero->SetFactory("P", tentativePFact00LevelZero);

  segregatedAFactLevelZero->SetFactory(dropMap1Name, slaveMapTransferFactLevelZero);
  segregatedAFactLevelZero->SetFactory(dropMap2Name, masterMapTransferFactLevelZero);

  /// Block (1,1)
  // SubBlockAFact for block (1,1)
  RCP<SubBlockAFactory> subBlockAFact11LevelZero = rcp(new SubBlockAFactory());
  subBlockAFact11LevelZero->SetFactory("A", MueLu::NoFactory::getRCP());
  subBlockAFact11LevelZero->SetParameter("block row", Teuchos::ParameterEntry(1));
  subBlockAFact11LevelZero->SetParameter("block col", Teuchos::ParameterEntry(1));

  /// Block (0,1) (off-diagonal) for InterfaceAggregationFactory
  // SubBlockAFact for block (0,1)
  RCP<SubBlockAFactory> subBlockAFact01LevelZero = rcp(new SubBlockAFactory);
  subBlockAFact01LevelZero->SetFactory("A", MueLu::NoFactory::getRCP());
  subBlockAFact01LevelZero->SetParameter("block row", Teuchos::ParameterEntry(0));
  subBlockAFact01LevelZero->SetParameter("block col", Teuchos::ParameterEntry(1));

  // define BlockedCoarseMapFactory
  RCP<BlockedCoarseMapFactory> coarseMapFact01LevelZero = rcp(new BlockedCoarseMapFactory());
  coarseMapFact01LevelZero->SetFactory("CoarseMap", coarseMapFact00LevelZero);

  // define MapTransferFactory for "Primal interface DOF map"
  RCP<MapTransferFactory> primalInterfaceDofMapTransferFactLevelZero = rcp(new MapTransferFactory());
  primalInterfaceDofMapTransferFactLevelZero->SetParameter("map: factory", Teuchos::ParameterEntry(
                                                                               std::string("Primal interface DOF map")));
  primalInterfaceDofMapTransferFactLevelZero->SetParameter("map: name", Teuchos::ParameterEntry(
                                                                            std::string("Primal interface DOF map")));
  primalInterfaceDofMapTransferFactLevelZero->SetFactory("P", tentativePFact00LevelZero);
  primalInterfaceDofMapTransferFactLevelZero->SetParameter("nullspace vectors: limit to",
                                                           Teuchos::ParameterEntry(std::string("translations")));

  // define factory for interface aggregates
  RCP<InterfaceAggregationFactory> interfaceAggFactLevelZero = rcp(new InterfaceAggregationFactory());
  interfaceAggFactLevelZero->SetFactory("A", subBlockAFact01LevelZero);
  interfaceAggFactLevelZero->SetFactory("Aggregates", uncoupledAggFactLevelZero);
  interfaceAggFactLevelZero->SetParameter("Dual/primal mapping strategy",
                                          Teuchos::ParameterEntry(std::string("dof-based")));
  interfaceAggFactLevelZero->SetFactory("Primal interface DOF map", primalInterfaceDofMapTransferFactLevelZero);

  // define tentativeP Factory
  RCP<TentativePFactory> tentativePFact01LevelZero = rcp(new TentativePFactory());
  tentativePFact01LevelZero->SetFactory("A", subBlockAFact11LevelZero);
  tentativePFact01LevelZero->SetFactory("Aggregates", interfaceAggFactLevelZero);
  tentativePFact01LevelZero->SetFactory("CoarseMap", coarseMapFact01LevelZero);
  tentativePFact01LevelZero->SetFactory("UnAmalgamationInfo", interfaceAggFactLevelZero);
  tentativePFact01LevelZero->SetParameter("tentative: calculate qr", Teuchos::ParameterEntry(true));
  tentativePFact01LevelZero->SetParameter("tentative: build coarse coordinates", Teuchos::ParameterEntry(false));

  // define Factory for nullspace2
  RCP<NullspaceFactory> nullspace2FactLevelZero = rcp(new NullspaceFactory());
  nullspace2FactLevelZero->SetParameter("Fine level nullspace", Teuchos::ParameterEntry(std::string("Nullspace2")));
  nullspace2FactLevelZero->SetFactory("Nullspace2", tentativePFact01LevelZero);

  // define CoarseMapFact01 dependencies
  coarseMapFact01LevelZero->SetFactory("Aggregates", interfaceAggFactLevelZero);
  coarseMapFact01LevelZero->SetFactory("CoarseMap", coarseMapFact00LevelZero);
  coarseMapFact01LevelZero->SetFactory("Nullspace", nullspace2FactLevelZero);

  /// Factory Managers
  // First group
  RCP<FactoryManager> M1LevelZero = rcp(new FactoryManager());
  M1LevelZero->SetFactory("A", subBlockAFact00LevelZero);
  M1LevelZero->SetFactory("P", tentativePFact00LevelZero);
  M1LevelZero->SetFactory("Aggregates", uncoupledAggFactLevelZero);
  M1LevelZero->SetFactory("Nullspace", nullspace1FactLevelZero);
  M1LevelZero->SetFactory("CoarseMap", coarseMapFact00LevelZero);
  M1LevelZero->SetFactory("UnAmalgamationInfo", amalgFact00LevelZero);
  M1LevelZero->SetFactory(dropMap1Name, slaveMapTransferFactLevelZero);
  M1LevelZero->SetFactory(dropMap2Name, masterMapTransferFactLevelZero);
  M1LevelZero->SetKokkosRefactor(false);

  // Second group
  RCP<FactoryManager> M2LevelZero = rcp(new FactoryManager());
  M2LevelZero->SetFactory("A", subBlockAFact11LevelZero);
  M2LevelZero->SetFactory("P", tentativePFact01LevelZero);
  M2LevelZero->SetFactory("Aggregates", interfaceAggFactLevelZero);
  M2LevelZero->SetFactory("Nullspace", nullspace2FactLevelZero);
  M2LevelZero->SetFactory("CoarseMap", coarseMapFact01LevelZero);
  M2LevelZero->SetFactory("Primal interface DOF map", primalInterfaceDofMapTransferFactLevelZero);
  M2LevelZero->SetKokkosRefactor(false);

  /// Blocked Transfer Operators
  RCP<BlockedPFactory> blockedPFactLevelZero = rcp(new BlockedPFactory);
  blockedPFactLevelZero->AddFactoryManager(M1LevelZero);
  blockedPFactLevelZero->AddFactoryManager(M2LevelZero);

  RCP<GenericRFactory> blockedRFactLevelZero = rcp(new GenericRFactory());
  blockedRFactLevelZero->SetFactory("P", blockedPFactLevelZero);

  RCP<BlockedRAPFactory> blockedRAPFactLevelZero = rcp(new BlockedRAPFactory());
  blockedRAPFactLevelZero->SetFactory("P", blockedPFactLevelZero);
  blockedRAPFactLevelZero->SetFactory("R", blockedRFactLevelZero);

  /// Main Factory Manager
  RCP<FactoryManager> managerLevelZero = rcp(new FactoryManager());
  managerLevelZero->SetFactory("P", blockedPFactLevelZero);
  managerLevelZero->SetFactory("R", blockedRFactLevelZero);
  managerLevelZero->SetFactory("A", blockedRAPFactLevelZero);
  managerLevelZero->SetFactory(dropMap1Name, MueLu::NoFactory::getRCP());
  managerLevelZero->SetFactory(dropMap2Name, MueLu::NoFactory::getRCP());
  managerLevelZero->SetKokkosRefactor(false);

  ///////////// Level Definition /////////////
  RCP<Level> levelZero = rcp(new Level());
  RCP<Level> levelOne  = rcp(new Level());

  levelZero->SetLevelID(0);
  levelOne->SetLevelID(1);
  levelOne->SetPreviousLevel(levelZero);

  TEUCHOS_ASSERT_EQUALITY(levelZero->GetLevelID(), 0);
  TEUCHOS_ASSERT_EQUALITY(levelOne->GetLevelID(), 1);

  levelZero->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(blockedMatrix));
  levelZero->Set("Nullspace1", nullspace1);
  levelZero->Set("Nullspace2", nullspace2);
  levelZero->Set("Primal interface DOF map", primalInterfaceDofMap);
  levelZero->Set(dropMap1Name, primalInterfaceDofMap);
  levelZero->Set(dropMap2Name, mastermap);

  TEUCHOS_ASSERT(levelZero->IsAvailable("A", MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable("Nullspace1", MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable("Nullspace2", MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable("Primal interface DOF map", MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable(dropMap1Name, MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable(dropMap2Name, MueLu::NoFactory::get()));

  levelZero->SetFactoryManager(managerLevelZero);
  //  levelOne->SetFactoryManager(managerLevelOne);

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  ///////////// Request and build level info /////////////
  // Check level 0
  levelZero->Request("Aggregates", uncoupledAggFactLevelZero.get());
  levelOne->Request("A", blockedRAPFactLevelZero.get());
  levelOne->Request("P", blockedPFactLevelZero.get());
  TEUCHOS_ASSERT(levelZero->IsRequested("Aggregates", uncoupledAggFactLevelZero.get()));
  TEUCHOS_ASSERT(levelOne->IsRequested("A", blockedRAPFactLevelZero.get()));
  TEUCHOS_ASSERT(levelOne->IsRequested("P", blockedPFactLevelZero.get()));
  blockedRAPFactLevelZero->Build(*levelZero, *levelOne);

  TEUCHOS_ASSERT(levelZero->IsAvailable("Aggregates", uncoupledAggFactLevelZero.get()));
  RCP<Aggregates> primalAggsLevelZero = levelZero->Get<RCP<Aggregates>>("Aggregates", uncoupledAggFactLevelZero.get());
  MueLuTests::checkAggregatesMapPair<SC, LO, GO, NO>(primalAggsLevelZero, stridingInfoPrimal,
                                                     levelZero->Get<RCP<const Map>>(dropMap1Name, MueLu::NoFactory::get()),
                                                     levelZero->Get<RCP<const Map>>(dropMap2Name, MueLu::NoFactory::get()));

  return EXIT_SUCCESS;
}

//-----------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}