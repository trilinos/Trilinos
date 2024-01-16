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
  std::vector<GO> slavemapEntries{42, 43, 44, 45, 52, 53};
  RCP<const Map> slavemap = Xpetra::MapFactory<LO, GO, NO>::Build(lib, slavemapEntries.size(), slavemapEntries, 0,
                                                                  comm);
  std::vector<GO> mastermapEntries{0, 1, 6, 7, 22, 23};
  RCP<const Map> mastermap = Xpetra::MapFactory<LO, GO, NO>::Build(lib, mastermapEntries.size(), mastermapEntries, 0,
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
  MueLuTests::checkAggregatesMapPair<SC, LO, GO, NO>(aggregates, stridingInfo, slavemap, mastermap);
  return EXIT_SUCCESS;
}

//-----------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
