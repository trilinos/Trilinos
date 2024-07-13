// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test
#include "contact_Helpers.hpp"

// MueLu
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CreateXpetraPreconditioner.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"

// Teuchos
#include <Teuchos_XMLParameterListHelpers.hpp>

// Xpetra
#include <Xpetra_IO.hpp>
#include <Xpetra_MatrixUtils.hpp>

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include "MueLu_UseShortNames.hpp"

  //  The SegregatedAFactory tests only work with real Scalar types,
  if (Teuchos::ScalarTraits<Scalar>::isComplex) return EXIT_SUCCESS;

  using GST = Teuchos::ScalarTraits<GO>;

  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  std::string probName = "";
  clp.setOption("probName", &probName, "Short name of the problem. Used to read-in the problem from files.");
  GO numTotalDofs = -GST::one();
  clp.setOption("nDof", &numTotalDofs, "Total number of DOFs");
  int nDofPerNode = -1;
  clp.setOption("nDofPerNode", &nDofPerNode, "Number of DOFS per Node");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(probName == "", MueLu::Exceptions::InvalidArgument, "Please specify a valid problem name.");
  TEUCHOS_TEST_FOR_EXCEPTION(numTotalDofs == -GST::one(), MueLu::Exceptions::InvalidArgument,
                             "Please specify the global number of DOFs on the command line.");

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(nDofPerNode);

  // Create the solid block dof row map
  std::string map1Name = "blockmap";
  std::vector<GO> blockmapEntries{30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
  RCP<const Map> blockmap = Xpetra::MapFactory<LO, GO, NO>::Build(lib, blockmapEntries.size(), blockmapEntries, 30, comm);

  // Create the interface dof row maps
  std::vector<GO> integration_side_map_entries{42, 43, 44, 45, 52, 53};
  RCP<const Map> integration_side_map = Xpetra::MapFactory<LO, GO, NO>::Build(lib, integration_side_map_entries.size(), integration_side_map_entries, 0,
                                                                              comm);
  std::vector<GO> projection_side_map_entries{0, 1, 6, 7, 22, 23};
  RCP<const Map> projection_side_map = Xpetra::MapFactory<LO, GO, NO>::Build(lib, projection_side_map_entries.size(), projection_side_map_entries, 0,
                                                                             comm);

  // Read the matrix from file
  const std::string matrixFileName = probName + "_matrix.mm";
  TEUCHOS_ASSERT(!matrixFileName.empty());
  RCP<Matrix> mat = Xpetra::IO<SC, LO, GO, NO>::Read(matrixFileName, lib, comm);
  Xpetra::MatrixUtils<SC, LO, GO, NO>::convertMatrixToStridedMaps(mat, stridingInfo, stridingInfo);

  // define amalgamation factory
  RCP<AmalgamationFactory> amalgFact = RCP(new AmalgamationFactory());

  // define CoalesceDropFactory that provides Graph
  RCP<CoalesceDropFactory> dropFact = RCP(new CoalesceDropFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory
  RCP<UncoupledAggregationFactory> uncoupledAggFact = rcp(new UncoupledAggregationFactory());
  uncoupledAggFact->SetFactory("Graph", dropFact);
  uncoupledAggFact->SetOrdering("graph");

  RCP<FactoryManager> M   = rcp(new FactoryManager());
  RCP<MueLu::Level> level = rcp(new MueLu::Level());
  level->SetLevelID(0);
  TEUCHOS_ASSERT_EQUALITY(level->GetLevelID(), 0);
  level->Set("A", mat);
  TEUCHOS_ASSERT(level->IsAvailable("A", MueLu::NoFactory::get()));
  level->SetFactoryManager(M);

  level->Request("Aggregates", uncoupledAggFact.get());
  level->Request(*uncoupledAggFact);
  uncoupledAggFact->Build(*level);
  RCP<Aggregates> aggregates = level->Get<RCP<Aggregates>>("Aggregates", uncoupledAggFact.get());

  typename Aggregates::LO_view aggPtr;
  typename Aggregates::LO_view aggNodes;
  typename Aggregates::LO_view unaggregated;

  aggregates->ComputeNodesInAggregate(aggPtr, aggNodes, unaggregated);

  MueLuTests::checkAggregatesBlockmap<SC, LO, GO, NO>(aggregates, stridingInfo, blockmap);
  MueLuTests::checkAggregatesMapPair<SC, LO, GO, NO>(aggregates, stridingInfo, integration_side_map, projection_side_map);
  return EXIT_SUCCESS;
}

//-----------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
