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
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CreateXpetraPreconditioner.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_MapTransferFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_SegregatedAFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"

// Teuchos
#include <Teuchos_XMLParameterListHelpers.hpp>

// Xpetra
#include <Xpetra_IO.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  using GST = Teuchos::ScalarTraits<GO>;

  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  std::string probName = ""; clp.setOption("probName", &probName, "Short name of the problem. Used to read-in the problem from files.");
  GO numTotalDofs = -GST::one(); clp.setOption("nDof", &numTotalDofs, "Total number of DOFs");
  int nDofPerNode = -1; clp.setOption("nDofPerNode", &nDofPerNode, "Number of DOFS per Node");

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

  TEUCHOS_TEST_FOR_EXCEPTION(probName=="", MueLu::Exceptions::InvalidArgument, "Please specify the problem name on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numTotalDofs==-GST::one(), MueLu::Exceptions::InvalidArgument, "Please specify the global number of DOF on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(nDofPerNode==-1, MueLu::Exceptions::InvalidArgument, "Please specify the number of DOF per mesh node on the command line.");

  const std::string matrixFileName = probName + "_matrix.mm";
  const std::string nullspaceFileName = probName + "_nullspace.mm";
  const std::string dropMap1Name = "dropMap1";
  const std::string dropMap2Name = "dropMap2";
  const std::string droppingScheme = "map-pair";

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(nDofPerNode);

  TEUCHOS_ASSERT(!matrixFileName.empty());
  RCP<Matrix> mat = Xpetra::IO<SC, LO, GO, NO>::Read(matrixFileName, lib, comm);
  Xpetra::MatrixUtils<SC, LO, GO, NO>::convertMatrixToStridedMaps(mat, stridingInfo, stridingInfo);
  RCP<const Map> rowmap = mat->getRowMap();
  TEUCHOS_ASSERT(!nullspaceFileName.empty());
  RCP<MultiVector> nullspace = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(nullspaceFileName, rowmap);

  // Create the interface dof row maps
  std::vector<GO> slavemapEntries{42, 43, 44, 45, 52, 53};
  RCP<const Map> dropMap1 = Xpetra::MapFactory<LO, GO, NO>::Build(lib, slavemapEntries.size(), slavemapEntries, 0,
                                                                            comm);
  std::vector<GO> mastermapEntries{0, 1, 6, 7, 22, 23};
  RCP<const Map> dropMap2 = Xpetra::MapFactory<LO, GO, NO>::Build(lib, mastermapEntries.size(), mastermapEntries, 0,
                                                                             comm);

  // define segregatedAFactory that provides Matrix A to CoalesceDropFactory
  RCP<SegregatedAFactory> segregatedAFact = RCP(new SegregatedAFactory());
  segregatedAFact->SetFactory("A", MueLu::NoFactory::getRCP());
  segregatedAFact->SetParameter("droppingScheme", Teuchos::ParameterEntry(droppingScheme));
  segregatedAFact->SetParameter("Call ReduceAll on dropMap1", Teuchos::ParameterEntry(true));
  segregatedAFact->SetParameter("Call ReduceAll on dropMap2", Teuchos::ParameterEntry(true));

  // define amalgamation factory
  RCP<AmalgamationFactory> amalgFact = RCP(new AmalgamationFactory());

  // define CoalesceDropFactory that provides Graph
  RCP<CoalesceDropFactory> dropFact = RCP(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  dropFact->SetFactory("A", segregatedAFact);

  // Setup aggregation factory
  RCP<UncoupledAggregationFactory> uncoupledAggFact = rcp(new UncoupledAggregationFactory());
  uncoupledAggFact->SetFactory("Graph", dropFact);
  uncoupledAggFact->SetFactory("DofsPerNode", dropFact);
  uncoupledAggFact->SetOrdering("graph");
  uncoupledAggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(2));
  uncoupledAggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));

  // Setup TentativePFactory
  RCP<TentativePFactory> tentativePFact = rcp(new TentativePFactory());
  tentativePFact->SetFactory("A", segregatedAFact);
  tentativePFact->SetFactory("Aggregates", uncoupledAggFact);

  // Setup transPFact for generating R
  RCP<TransPFactory> transPFact = rcp(new TransPFactory());
  transPFact->SetFactory("P", tentativePFact);

  // Setup nullspace Factory
  RCP<NullspaceFactory> nullspaceFact = rcp(new NullspaceFactory());
  nullspaceFact->SetParameter("Fine level nullspace", Teuchos::ParameterEntry(std::string("Nullspace")));
  nullspaceFact->SetFactory("Nullspace", tentativePFact);

  // Setup RAP Factory
  RCP<RAPFactory> rapFact = rcp(new RAPFactory());
  rapFact->SetFactory("A", MueLu::NoFactory::getRCP());
  rapFact->SetFactory("P", tentativePFact);
  rapFact->SetFactory("R", transPFact);

  // Setup Map Transfer Factory for dropMap1
  RCP<MapTransferFactory> slaveMapTransferFact = Teuchos::rcp(new MapTransferFactory());
  slaveMapTransferFact->SetParameter("map: factory", Teuchos::ParameterEntry(dropMap1Name));
  slaveMapTransferFact->SetParameter("map: name", Teuchos::ParameterEntry(dropMap1Name));
  slaveMapTransferFact->SetFactory("P", tentativePFact);
  rapFact->AddTransferFactory(slaveMapTransferFact);
  segregatedAFact->SetFactory(dropMap1Name, slaveMapTransferFact);

  // Setup Map Transfer Factory for dropMap2
  RCP<MapTransferFactory> masterMapTransferFact = Teuchos::rcp(new MapTransferFactory());
  masterMapTransferFact->SetParameter("map: factory", Teuchos::ParameterEntry(dropMap2Name));
  masterMapTransferFact->SetParameter("map: name", Teuchos::ParameterEntry(dropMap2Name));
  masterMapTransferFact->SetFactory("P", tentativePFact);
  rapFact->AddTransferFactory(masterMapTransferFact);
  segregatedAFact->SetFactory(dropMap2Name, masterMapTransferFact);

  // Main Factory manager
  RCP<FactoryManager> M = rcp(new FactoryManager());
  M->SetFactory("A", rapFact);
  M->SetFactory("Graph", dropFact);
  M->SetFactory("UnAmalgamationInfo", amalgFact);
  M->SetFactory("Nullspace", nullspaceFact);
  M->SetFactory("P", tentativePFact);
  M->SetFactory("R", transPFact);
  M->SetFactory("Aggregates", uncoupledAggFact);
  M->SetFactory(dropMap1Name, MueLu::NoFactory::getRCP());
  M->SetFactory(dropMap2Name, MueLu::NoFactory::getRCP());

  ///////////// Level Definition /////////////
  RCP<Level> levelZero = rcp(new Level());
  RCP<Level> levelOne = rcp(new Level());
  RCP<Level> levelTwo = rcp(new Level());

  levelZero->SetLevelID(0);
  levelOne->SetLevelID(1);
  levelTwo->SetLevelID(2);

  levelOne->SetPreviousLevel(levelZero);
  levelTwo->SetPreviousLevel(levelOne);

  TEUCHOS_ASSERT_EQUALITY(levelZero->GetLevelID(), 0);
  TEUCHOS_ASSERT_EQUALITY(levelOne->GetLevelID(), 1);
  TEUCHOS_ASSERT_EQUALITY(levelTwo->GetLevelID(), 2);

  levelZero->Set("A", mat);
  levelZero->Set(dropMap1Name, dropMap1);
  levelZero->Set(dropMap2Name, dropMap2);
  levelZero->Set("Nullspace", nullspace);

  TEUCHOS_ASSERT(levelZero->IsAvailable("A", MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable("Nullspace", MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable(dropMap1Name, MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelZero->IsAvailable(dropMap2Name, MueLu::NoFactory::get()));

  levelZero->SetFactoryManager(M);
  levelOne->SetFactoryManager(M);

  //// level 0
  levelZero->Request("Aggregates", uncoupledAggFact.get());
  levelZero->Request(*uncoupledAggFact);
  TEUCHOS_ASSERT(levelZero->IsRequested("Aggregates", uncoupledAggFact.get()));
  uncoupledAggFact->Build(*levelZero);
  TEUCHOS_ASSERT(levelZero->IsAvailable("Aggregates", uncoupledAggFact.get()));
  RCP<Aggregates> aggsLevelZero = levelZero->Get<RCP<Aggregates>>("Aggregates", uncoupledAggFact.get());

  MueLuTests::checkAggregatesMapPair<SC, LO, GO, NO>(aggsLevelZero, stridingInfo,
                                                     levelZero->Get<RCP<const Map>>(dropMap1Name, MueLu::NoFactory::get()),
                                                     levelZero->Get<RCP<const Map>>(dropMap2Name, MueLu::NoFactory::get()));

  //// level 1
  // Define factories for level one
  // define Segregated A factory
  RCP<SegregatedAFactory> segregatedAFactLevelOne = rcp(new SegregatedAFactory());
  segregatedAFactLevelOne->SetFactory("A", rapFact);
  segregatedAFactLevelOne->SetParameter("droppingScheme", Teuchos::ParameterEntry(droppingScheme));
  segregatedAFactLevelOne->SetParameter("Call ReduceAll on dropMap1", Teuchos::ParameterEntry(true));
  segregatedAFactLevelOne->SetParameter("Call ReduceAll on dropMap2", Teuchos::ParameterEntry(true));

  // define amalgamation factory
  RCP<AmalgamationFactory> amalgFactLevelOne = RCP(new AmalgamationFactory());

  // define CoalesceDropFactory that provides Graph
  RCP<CoalesceDropFactory> dropFactLevelOne = RCP(new CoalesceDropFactory());
  dropFactLevelOne->SetFactory("UnAmalgamationInfo", amalgFactLevelOne);
  dropFactLevelOne->SetFactory("A", segregatedAFactLevelOne);

  // Setup aggregation factory
  RCP<UncoupledAggregationFactory> uncoupledAggFactLevelOne = rcp(new UncoupledAggregationFactory());
  uncoupledAggFactLevelOne->SetFactory("Graph", dropFactLevelOne);
  uncoupledAggFactLevelOne->SetFactory("DofsPerNode", dropFactLevelOne);
  uncoupledAggFactLevelOne->SetOrdering("graph");
  uncoupledAggFactLevelOne->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(2));
  uncoupledAggFactLevelOne->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));

  // Request level info
  levelOne->Request("A", rapFact.get());
  TEUCHOS_ASSERT(levelOne->IsRequested("A", rapFact.get()));
  rapFact->Build(*levelZero, *levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable("A", rapFact.get()));

  levelOne->Request(dropMap1Name, MueLu::NoFactory::get());
  levelOne->Request(dropMap2Name, MueLu::NoFactory::get());
  TEUCHOS_ASSERT(levelOne->IsRequested(dropMap1Name, MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelOne->IsRequested(dropMap2Name, MueLu::NoFactory::get()));
  slaveMapTransferFact->Build(*levelZero, *levelOne);
  masterMapTransferFact->Build(*levelZero, *levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable(dropMap1Name, MueLu::NoFactory::get()));
  TEUCHOS_ASSERT(levelOne->IsAvailable(dropMap2Name, MueLu::NoFactory::get()));

  levelOne->Request("A", segregatedAFactLevelOne.get());
  levelOne->Request(*segregatedAFactLevelOne);
  segregatedAFactLevelOne->Build(*levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable("A", segregatedAFactLevelOne.get()));

  // AmalgamationFactory of block A00
  levelOne->Request("UnAmalgamationInfo", amalgFactLevelOne.get());
  levelOne->Request(*amalgFactLevelOne);
  amalgFactLevelOne->Build(*levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable("UnAmalgamationInfo", amalgFactLevelOne.get()));

  levelOne->Request("Aggregates", uncoupledAggFactLevelOne.get());
  TEUCHOS_ASSERT(levelOne->IsRequested("Aggregates", uncoupledAggFactLevelOne.get()));
  uncoupledAggFactLevelOne->Build(*levelOne);
  TEUCHOS_ASSERT(levelOne->IsAvailable("Aggregates", uncoupledAggFactLevelOne.get()));

  RCP<Aggregates> primalAggsLevelOne = levelOne->Get<RCP<Aggregates>>("Aggregates", uncoupledAggFactLevelOne.get());

  MueLuTests::checkAggregatesMapPair<SC, LO, GO, NO>(aggsLevelZero, stridingInfo,
                                                     levelZero->Get<RCP<const Map>>(dropMap1Name, MueLu::NoFactory::get()),
                                                     levelZero->Get<RCP<const Map>>(dropMap2Name, MueLu::NoFactory::get()));
  return EXIT_SUCCESS;
}

//-----------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}