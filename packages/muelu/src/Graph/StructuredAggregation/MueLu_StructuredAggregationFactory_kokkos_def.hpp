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
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_STRUCTUREDAGGREGATIONFACTORY_KOKKOS_DEF_HPP
#define MUELU_STRUCTUREDAGGREGATIONFACTORY_KOKKOS_DEF_HPP

// Xpetra includes
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsGraph.hpp>

// MueLu generic includes
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

// MueLu specific includes (kokkos version)
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_IndexManager_kokkos.hpp"
#include "MueLu_AggregationStructuredAlgorithm_kokkos.hpp"

#include "MueLu_StructuredAggregationFactory_kokkos_decl.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredAggregationFactory_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    StructuredAggregationFactory_kokkos()
  : bDefinitionPhase_(true) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> StructuredAggregationFactory_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: preserve Dirichlet points");
  SET_VALID_ENTRY("aggregation: allow user-specified singletons");
  SET_VALID_ENTRY("aggregation: error on nodes with no on-rank neighbors");
  SET_VALID_ENTRY("aggregation: phase3 avoid singletons");
#undef SET_VALID_ENTRY

  // general variables needed in StructuredAggregationFactory
  validParamList->set<std::string>("aggregation: output type", "Aggregates",
                                   "Type of object holding the aggregation data: Aggregtes or CrsGraph");
  validParamList->set<std::string>("aggregation: coarsening rate", "{3}",
                                   "Coarsening rate per spatial dimensions");
  validParamList->set<int>("aggregation: coarsening order", 0,
                           "The interpolation order used to construct grid transfer operators based off these aggregates.");
  validParamList->set<RCP<const FactoryBase> >("Graph", Teuchos::null,
                                               "Graph of the matrix after amalgamation but without dropping.");
  validParamList->set<RCP<const FactoryBase> >("DofsPerNode", Teuchos::null,
                                               "Number of degrees of freedom per mesh node, provided by the coalsce drop factory.");
  validParamList->set<RCP<const FactoryBase> >("numDimensions", Teuchos::null,
                                               "Number of spatial dimension provided by CoordinatesTransferFactory.");
  validParamList->set<RCP<const FactoryBase> >("lNodesPerDim", Teuchos::null,
                                               "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");

  return validParamList;
}  // GetValidParameterList()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredAggregationFactory_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "Graph");
  Input(currentLevel, "DofsPerNode");

  // Request the local number of nodes per dimensions
  if (currentLevel.GetLevelID() == 0) {
    if (currentLevel.IsAvailable("numDimensions", NoFactory::get())) {
      currentLevel.DeclareInput("numDimensions", NoFactory::get(), this);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("numDimensions", NoFactory::get()),
                                 Exceptions::RuntimeError,
                                 "numDimensions was not provided by the user on level0!");
    }
    if (currentLevel.IsAvailable("lNodesPerDim", NoFactory::get())) {
      currentLevel.DeclareInput("lNodesPerDim", NoFactory::get(), this);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("lNodesPerDim", NoFactory::get()),
                                 Exceptions::RuntimeError,
                                 "lNodesPerDim was not provided by the user on level0!");
    }
  } else {
    Input(currentLevel, "lNodesPerDim");
    Input(currentLevel, "numDimensions");
  }
}  // DeclareInput()

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredAggregationFactory_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
    Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_STRUCTUREDAGGREGATION_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  using device_type     = typename LWGraph_kokkos::local_graph_type::device_type;
  using execution_space = typename LWGraph_kokkos::local_graph_type::device_type::execution_space;
  using memory_space    = typename LWGraph_kokkos::local_graph_type::device_type::memory_space;

  *out << "Entering structured aggregation" << std::endl;

  ParameterList pL  = GetParameterList();
  bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

  // General problem informations are gathered from data stored in the problem matix.
  RCP<const LWGraph_kokkos> graph = Get<RCP<LWGraph_kokkos> >(currentLevel, "Graph");
  RCP<const Map> fineMap          = graph->GetDomainMap();
  const int myRank                = fineMap->getComm()->getRank();
  const LO dofsPerNode            = Get<LO>(currentLevel, "DofsPerNode");

  // Since we want to operate on nodes and not dof, we need to modify the rowMap in order to
  // obtain a nodeMap.
  const int interpolationOrder = pL.get<int>("aggregation: coarsening order");
  std::string outputType       = pL.get<std::string>("aggregation: output type");
  const bool outputAggregates  = (outputType == "Aggregates" ? true : false);
  Array<LO> lFineNodesPerDir(3);
  int numDimensions;
  if (currentLevel.GetLevelID() == 0) {
    // On level 0, data is provided by applications and has no associated factory.
    lFineNodesPerDir = currentLevel.Get<Array<LO> >("lNodesPerDim", NoFactory::get());
    numDimensions    = currentLevel.Get<int>("numDimensions", NoFactory::get());
  } else {
    // On level > 0, data is provided directly by generating factories.
    lFineNodesPerDir = Get<Array<LO> >(currentLevel, "lNodesPerDim");
    numDimensions    = Get<int>(currentLevel, "numDimensions");
  }

  // First make sure that input parameters are set logically based on dimension
  for (int dim = 0; dim < 3; ++dim) {
    if (dim >= numDimensions) {
      lFineNodesPerDir[dim] = 1;
    }
  }

  // Get the coarsening rate
  std::string coarseningRate = pL.get<std::string>("aggregation: coarsening rate");
  Teuchos::Array<LO> coarseRate;
  try {
    coarseRate = Teuchos::fromStringToArray<LO>(coarseningRate);
  } catch (const Teuchos::InvalidArrayStringRepresentation& e) {
    GetOStream(Errors, -1) << " *** \"aggregation: coarsening rate\" must be a string convertible into an array! *** "
                           << std::endl;
    throw e;
  }
  TEUCHOS_TEST_FOR_EXCEPTION((coarseRate.size() > 1) && (coarseRate.size() < numDimensions),
                             Exceptions::RuntimeError,
                             "\"aggregation: coarsening rate\" must have at least as many"
                             " components as the number of spatial dimensions in the problem.");

  // Now that we have extracted info from the level, create the IndexManager
  RCP<IndexManager_kokkos> geoData = rcp(new IndexManager_kokkos(numDimensions,
                                                                 interpolationOrder, myRank,
                                                                 lFineNodesPerDir,
                                                                 coarseRate));

  *out << "The index manager has now been built" << std::endl;
  TEUCHOS_TEST_FOR_EXCEPTION(fineMap->getLocalNumElements() != static_cast<size_t>(geoData->getNumLocalFineNodes()),
                             Exceptions::RuntimeError,
                             "The local number of elements in the graph's map is not equal to "
                             "the number of nodes given by: lNodesPerDim!");

  // Now we are ready for the big loop over the fine node that will assign each
  // node on the fine grid to an aggregate and a processor.
  RCP<AggregationStructuredAlgorithm_kokkos> myStructuredAlgorithm = rcp(new AggregationStructuredAlgorithm_kokkos());

  if (interpolationOrder == 0 && outputAggregates) {
    RCP<Aggregates> aggregates = rcp(new Aggregates(graph->GetDomainMap()));
    aggregates->setObjectLabel("ST");
    aggregates->SetIndexManagerKokkos(geoData);
    aggregates->AggregatesCrossProcessors(false);
    aggregates->SetNumAggregates(geoData->getNumCoarseNodes());

    LO numNonAggregatedNodes = geoData->getNumLocalFineNodes();
    Kokkos::View<unsigned*, device_type> aggStat("aggStat", numNonAggregatedNodes);
    Kokkos::parallel_for(
        "StructuredAggregation: initialize aggStat",
        Kokkos::RangePolicy<execution_space>(0, numNonAggregatedNodes),
        KOKKOS_LAMBDA(const LO nodeIdx) { aggStat(nodeIdx) = READY; });

    myStructuredAlgorithm->BuildAggregates(pL, *graph, *aggregates, aggStat,
                                           numNonAggregatedNodes);

    *out << "numNonAggregatedNodes: " << numNonAggregatedNodes << std::endl;

    TEUCHOS_TEST_FOR_EXCEPTION(numNonAggregatedNodes, Exceptions::RuntimeError,
                               "MueLu::StructuredAggregationFactory::Build: Leftover nodes found! Error!");
    aggregates->ComputeAggregateSizes(true /*forceRecompute*/);
    GetOStream(Statistics1) << aggregates->description() << std::endl;
    Set(currentLevel, "Aggregates", aggregates);

  } else {
    // Create Coarse Data
    RCP<CrsGraph> myGraph;
    myStructuredAlgorithm->BuildGraph(*graph, geoData, dofsPerNode, myGraph);
    Set(currentLevel, "prolongatorGraph", myGraph);
  }

  Set(currentLevel, "lCoarseNodesPerDim", geoData->getCoarseNodesPerDirArray());
  Set(currentLevel, "indexManager", geoData);
  Set(currentLevel, "structuredInterpolationOrder", interpolationOrder);
  Set(currentLevel, "numDimensions", numDimensions);

}  // Build()

}  // namespace MueLu

#endif /* MUELU_STRUCTUREDAGGREGATIONFACTORY_KOKKOS_DEF_HPP */
