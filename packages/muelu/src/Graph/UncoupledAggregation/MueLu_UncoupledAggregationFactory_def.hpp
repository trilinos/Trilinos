// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_
#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_

#include <climits>

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_UncoupledAggregationFactory_decl.hpp"

#include "MueLu_InterfaceAggregationAlgorithm.hpp"
#include "MueLu_OnePtAggregationAlgorithm.hpp"
#include "MueLu_PreserveDirichletAggregationAlgorithm.hpp"

#include "MueLu_AggregationPhase1Algorithm.hpp"
#include "MueLu_AggregationPhase2aAlgorithm.hpp"
#include "MueLu_AggregationPhase2bAlgorithm.hpp"
#include "MueLu_AggregationPhase3Algorithm.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

#include "KokkosGraph_Distance2ColorHandle.hpp"
#include "KokkosGraph_Distance2Color.hpp"
#include "KokkosGraph_MIS2.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::UncoupledAggregationFactory()
  : bDefinitionPhase_(true) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::~UncoupledAggregationFactory() = default;

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  // Aggregation parameters (used in aggregation algorithms)
  // TODO introduce local member function for each aggregation algorithm such that each aggregation algorithm can define its own parameters

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: max agg size");
  SET_VALID_ENTRY("aggregation: min agg size");
  SET_VALID_ENTRY("aggregation: max selected neighbors");
  SET_VALID_ENTRY("aggregation: ordering");
  validParamList->getEntry("aggregation: ordering").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("natural", "graph", "random"))));
  SET_VALID_ENTRY("aggregation: deterministic");
  SET_VALID_ENTRY("aggregation: coloring algorithm");
  SET_VALID_ENTRY("aggregation: enable phase 1");
  SET_VALID_ENTRY("aggregation: enable phase 2a");
  SET_VALID_ENTRY("aggregation: enable phase 2b");
  SET_VALID_ENTRY("aggregation: enable phase 3");
  SET_VALID_ENTRY("aggregation: match ML phase1");
  SET_VALID_ENTRY("aggregation: match ML phase2a");
  SET_VALID_ENTRY("aggregation: match ML phase2b");
  SET_VALID_ENTRY("aggregation: phase2a agg factor");
  SET_VALID_ENTRY("aggregation: preserve Dirichlet points");
  SET_VALID_ENTRY("aggregation: allow user-specified singletons");
  SET_VALID_ENTRY("aggregation: use interface aggregation");
  SET_VALID_ENTRY("aggregation: error on nodes with no on-rank neighbors");
  SET_VALID_ENTRY("aggregation: phase3 avoid singletons");
  SET_VALID_ENTRY("aggregation: compute aggregate qualities");
  SET_VALID_ENTRY("aggregation: phase 1 algorithm");
#undef SET_VALID_ENTRY

  // general variables needed in AggregationFactory
  validParamList->set<RCP<const FactoryBase>>("Graph", null, "Generating factory of the graph");
  validParamList->set<RCP<const FactoryBase>>("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");
  validParamList->set<RCP<const FactoryBase>>("AggregateQualities", null, "Generating factory for variable \'AggregateQualities\'");

  // special variables necessary for OnePtAggregationAlgorithm
  validParamList->set<std::string>("OnePt aggregate map name", "", "Name of input map for single node aggregates. (default='')");
  validParamList->set<std::string>("OnePt aggregate map factory", "", "Generating factory of (DOF) map for single node aggregates.");
  // validParamList->set< RCP<const FactoryBase> >("OnePt aggregate map factory",    NoFactory::getRCP(), "Generating factory of (DOF) map for single node aggregates.");

  // InterfaceAggregation parameters
  // validParamList->set< bool >                  ("aggregation: use interface aggregation", "false", "Flag to trigger aggregation along an interface using specified aggregate seeds.");
  validParamList->set<std::string>("Interface aggregate map name", "", "Name of input map for interface aggregates. (default='')");
  validParamList->set<std::string>("Interface aggregate map factory", "", "Generating factory of (DOF) map for interface aggregates.");
  validParamList->set<RCP<const FactoryBase>>("nodeOnInterface", Teuchos::null, "Array specifying whether or not a node is on the interface (1 or 0).");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "Graph");
  Input(currentLevel, "DofsPerNode");

  const ParameterList& pL = GetParameterList();

  // request special data necessary for OnePtAggregationAlgorithm
  std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name");
  if (mapOnePtName.length() > 0) {
    std::string mapOnePtFactName = pL.get<std::string>("OnePt aggregate map factory");
    if (mapOnePtFactName == "" || mapOnePtFactName == "NoFactory") {
      currentLevel.DeclareInput(mapOnePtName, NoFactory::get());
    } else {
      RCP<const FactoryBase> mapOnePtFact = GetFactory(mapOnePtFactName);
      currentLevel.DeclareInput(mapOnePtName, mapOnePtFact.get());
    }
  }

  // request special data necessary for InterfaceAggregation
  if (pL.get<bool>("aggregation: use interface aggregation") == true) {
    if (currentLevel.GetLevelID() == 0) {
      if (currentLevel.IsAvailable("nodeOnInterface", NoFactory::get())) {
        currentLevel.DeclareInput("nodeOnInterface", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("nodeOnInterface", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "nodeOnInterface was not provided by the user on level0!");
      }
    } else {
      Input(currentLevel, "nodeOnInterface");
    }
  }

  if (pL.get<bool>("aggregation: compute aggregate qualities")) {
    Input(currentLevel, "AggregateQualities");
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  ParameterList pL  = GetParameterList();
  bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

  if (pL.get<int>("aggregation: max agg size") == -1)
    pL.set("aggregation: max agg size", INT_MAX);

  // define aggregation algorithms
  RCP<const FactoryBase> graphFact = GetFactory("Graph");

  // TODO Can we keep different aggregation algorithms over more Build calls?
  algos_.clear();
  algos_.push_back(rcp(new PreserveDirichletAggregationAlgorithm(graphFact)));
  if (pL.get<bool>("aggregation: use interface aggregation") == true) algos_.push_back(rcp(new InterfaceAggregationAlgorithm(graphFact)));
  if (pL.get<bool>("aggregation: allow user-specified singletons") == true) algos_.push_back(rcp(new OnePtAggregationAlgorithm(graphFact)));
  if (pL.get<bool>("aggregation: enable phase 1") == true) algos_.push_back(rcp(new AggregationPhase1Algorithm(graphFact)));
  if (pL.get<bool>("aggregation: enable phase 2a") == true) algos_.push_back(rcp(new AggregationPhase2aAlgorithm(graphFact)));
  if (pL.get<bool>("aggregation: enable phase 2b") == true) algos_.push_back(rcp(new AggregationPhase2bAlgorithm(graphFact)));
  if (pL.get<bool>("aggregation: enable phase 3") == true) algos_.push_back(rcp(new AggregationPhase3Algorithm(graphFact)));

  std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name");
  RCP<Map> OnePtMap        = Teuchos::null;
  if (mapOnePtName.length()) {
    std::string mapOnePtFactName = pL.get<std::string>("OnePt aggregate map factory");
    if (mapOnePtFactName == "" || mapOnePtFactName == "NoFactory") {
      OnePtMap = currentLevel.Get<RCP<Map>>(mapOnePtName, NoFactory::get());
    } else {
      RCP<const FactoryBase> mapOnePtFact = GetFactory(mapOnePtFactName);
      OnePtMap                            = currentLevel.Get<RCP<Map>>(mapOnePtName, mapOnePtFact.get());
    }
  }

  // Set map for interface aggregates
  std::string mapInterfaceName = pL.get<std::string>("Interface aggregate map name");
  RCP<Map> InterfaceMap        = Teuchos::null;

  RCP<const LWGraph> graph;
  RCP<const LWGraph_kokkos> graph_kokkos;
  RCP<Aggregates> aggregates;
  RCP<const Teuchos::Comm<int>> comm;
  LO numRows;
  bool runOnHost;
  if (IsType<RCP<LWGraph>>(currentLevel, "Graph")) {
    graph      = Get<RCP<LWGraph>>(currentLevel, "Graph");
    aggregates = rcp(new Aggregates(*graph));
    comm       = graph->GetComm();
    numRows    = graph->GetNodeNumVertices();
    runOnHost  = true;
  } else {
    graph_kokkos = Get<RCP<LWGraph_kokkos>>(currentLevel, "Graph");
    aggregates   = rcp(new Aggregates(*graph_kokkos));
    comm         = graph_kokkos->GetComm();
    numRows      = graph_kokkos->GetNodeNumVertices();
    runOnHost    = false;

    TEUCHOS_TEST_FOR_EXCEPTION(pL.get<bool>("aggregation: use interface aggregation"), std::invalid_argument, "Option: 'aggregation: use interface aggregation' is not supported in the Kokkos version of uncoupled aggregation");
    // Sanity Checking: match ML behavior is not supported in UncoupledAggregation_Kokkos in Phase 1 or Phase 2b, but is in 2a
    TEUCHOS_TEST_FOR_EXCEPTION(pL.get<bool>("aggregation: match ML phase1"), std::invalid_argument, "Option: 'aggregation: match ML phase1' is not supported in the Kokkos version of uncoupled aggregation");
    TEUCHOS_TEST_FOR_EXCEPTION(pL.get<bool>("aggregation: match ML phase2b"), std::invalid_argument, "Option: 'aggregation: match ML phase2b' is not supported in the Kokkos version of uncoupled aggregation");
  }

  // Build
  aggregates->setObjectLabel("UC");

  // construct aggStat information
  using AggStatHostType = typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType;
  using AggStatType     = typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType;
  AggStatHostType aggStatHost;
  AggStatType aggStat;

  if (runOnHost) {
    aggStatHost = AggStatHostType(Kokkos::ViewAllocateWithoutInitializing("aggregation status"), numRows);
    Kokkos::deep_copy(aggStatHost, READY);
  } else {
    aggStat = AggStatType(Kokkos::ViewAllocateWithoutInitializing("aggregation status"), numRows);
    Kokkos::deep_copy(aggStat, READY);
  }

  // interface
  if (pL.get<bool>("aggregation: use interface aggregation") == true) {
    Teuchos::Array<LO> nodeOnInterface = Get<Array<LO>>(currentLevel, "nodeOnInterface");
    for (LO i = 0; i < numRows; i++) {
      if (nodeOnInterface[i])
        aggStatHost[i] = INTERFACE;
    }
  }

  // Dirichlet nodes
  {
    if (runOnHost) {
      auto dirichletBoundaryMap = graph->GetBoundaryNodeMap();
      Kokkos::parallel_for(
          "MueLu - UncoupledAggregation: tagging boundary nodes in aggStat",
          Kokkos::RangePolicy<LocalOrdinal, typename LWGraph::execution_space>(0, numRows),
          KOKKOS_LAMBDA(const LocalOrdinal nodeIdx) {
            if (dirichletBoundaryMap(nodeIdx) == true) {
              aggStatHost(nodeIdx) = BOUNDARY;
            }
          });
    } else {
      auto dirichletBoundaryMap = graph_kokkos->GetBoundaryNodeMap();
      Kokkos::parallel_for(
          "MueLu - UncoupledAggregation: tagging boundary nodes in aggStat",
          Kokkos::RangePolicy<LocalOrdinal, typename LWGraph_kokkos::execution_space>(0, numRows),
          KOKKOS_LAMBDA(const LocalOrdinal nodeIdx) {
            if (dirichletBoundaryMap(nodeIdx) == true) {
              aggStat(nodeIdx) = BOUNDARY;
            }
          });
    }
  }

  if (OnePtMap != Teuchos::null) {
    LO nDofsPerNode = Get<LO>(currentLevel, "DofsPerNode");

    if (runOnHost) {
      GO indexBase = graph->GetDomainMap()->getIndexBase();
      for (LO i = 0; i < numRows; i++) {
        // reconstruct global row id (FIXME only works for contiguous maps)
        GO grid = (graph->GetDomainMap()->getGlobalElement(i) - indexBase) * nDofsPerNode + indexBase;

        for (LO kr = 0; kr < nDofsPerNode; kr++)
          if (OnePtMap->isNodeGlobalElement(grid + kr))
            aggStatHost(i) = ONEPT;
      }
    } else {
      GO indexBase               = graph_kokkos->GetDomainMap()->getIndexBase();
      auto lclDomainMap          = graph_kokkos->GetDomainMap()->getLocalMap();
      auto lclOnePtMap           = OnePtMap->getLocalMap();
      const LocalOrdinal INVALID = Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid();
      Kokkos::parallel_for(
          "MueLu - UncoupledAggregation: tagging OnePt map",
          Kokkos::RangePolicy<LocalOrdinal, typename LWGraph_kokkos::execution_space>(0, numRows),
          KOKKOS_LAMBDA(const LocalOrdinal i) {
            // reconstruct global row id (FIXME only works for contiguous maps)
            GO grid = (lclDomainMap.getGlobalElement(i) - indexBase) * nDofsPerNode + indexBase;

            for (LO kr = 0; kr < nDofsPerNode; kr++)
              if (lclOnePtMap.getLocalElement(grid + kr) != INVALID)
                aggStat(i) = ONEPT;
          });
    }
  }

  LO numNonAggregatedNodes = numRows;
  std::string aggAlgo      = pL.get<std::string>("aggregation: coloring algorithm");
  if (aggAlgo == "mis2 coarsening" || aggAlgo == "mis2 aggregation") {
    TEUCHOS_ASSERT(!runOnHost);

    SubFactoryMonitor sfm(*this, "Algo \"MIS2\"", currentLevel);
    using graph_t      = typename LWGraph_kokkos::local_graph_type;
    using device_t     = typename graph_t::device_type;
    using exec_space   = typename device_t::execution_space;
    using rowmap_t     = typename graph_t::row_map_type;
    using colinds_t    = typename graph_t::entries_type;
    using lno_t        = typename colinds_t::non_const_value_type;
    rowmap_t aRowptrs  = graph_kokkos->getRowPtrs();
    colinds_t aColinds = graph_kokkos->getEntries();
    lno_t numAggs      = 0;
    typename colinds_t::non_const_type labels;

    if (aggAlgo == "mis2 coarsening") {
      if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: MIS-2 coarsening" << std::endl;
      labels = KokkosGraph::graph_mis2_coarsen<device_t, rowmap_t, colinds_t>(aRowptrs, aColinds, numAggs);
    } else if (aggAlgo == "mis2 aggregation") {
      if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: MIS-2 aggregation" << std::endl;
      labels = KokkosGraph::graph_mis2_aggregate<device_t, rowmap_t, colinds_t>(aRowptrs, aColinds, numAggs);
    }
    auto vertex2AggId = aggregates->GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
    auto procWinner   = aggregates->GetProcWinner()->getDeviceLocalView(Xpetra::Access::OverwriteAll);
    int rank          = comm->getRank();
    Kokkos::parallel_for(
        Kokkos::RangePolicy<exec_space>(0, numRows),
        KOKKOS_LAMBDA(lno_t i) {
          procWinner(i, 0) = rank;
          if (aggStat(i) == READY) {
            aggStat(i)         = AGGREGATED;
            vertex2AggId(i, 0) = labels(i);
          }
        });
    numNonAggregatedNodes = 0;
    aggregates->SetNumAggregates(numAggs);
  } else {
    if (!runOnHost) {
      DoGraphColoring(currentLevel, aggAlgo, pL.get<bool>("aggregation: deterministic"), graph_kokkos, aggregates);
      if (IsPrint(Statistics1)) {
        GetOStream(Statistics1) << "  num colors: " << aggregates->GetGraphNumColors() << std::endl;
      }
    }

    GO numGlobalRows           = 0;
    GO numGlobalAggregatedPrev = 0, numGlobalAggsPrev = 0;
    if (IsPrint(Statistics1))
      MueLu_sumAll(comm, as<GO>(numRows), numGlobalRows);
    for (size_t a = 0; a < algos_.size(); a++) {
      std::string phase = algos_[a]->description();
      SubFactoryMonitor sfm2(*this, "Algo \"" + phase + "\"", currentLevel);

      int oldRank = algos_[a]->SetProcRankVerbose(this->GetProcRankVerbose());
      if (runOnHost)
        algos_[a]->BuildAggregatesNonKokkos(pL, *graph, *aggregates, aggStatHost, numNonAggregatedNodes);
      else
        algos_[a]->BuildAggregates(pL, *graph_kokkos, *aggregates, aggStat, numNonAggregatedNodes);
      algos_[a]->SetProcRankVerbose(oldRank);

      if (IsPrint(Statistics1)) {
        GO numLocalAggregated = numRows - numNonAggregatedNodes, numGlobalAggregated = 0;
        GO numLocalAggs = aggregates->GetNumAggregates(), numGlobalAggs = 0;
        MueLu_sumAll(comm, numLocalAggregated, numGlobalAggregated);
        MueLu_sumAll(comm, numLocalAggs, numGlobalAggs);

        double aggPercent = 100 * as<double>(numGlobalAggregated) / as<double>(numGlobalRows);
        if (aggPercent > 99.99 && aggPercent < 100.00) {
          // Due to round off (for instance, for 140465733/140466897), we could
          // get 100.00% display even if there are some remaining nodes. This
          // is bad from the users point of view. It is much better to change
          // it to display 99.99%.
          aggPercent = 99.99;
        }
        GetOStream(Statistics1) << "  aggregated : " << (numGlobalAggregated - numGlobalAggregatedPrev) << " (phase), " << std::fixed
                                << std::setprecision(2) << numGlobalAggregated << "/" << numGlobalRows << " [" << aggPercent << "%] (total)\n"
                                << "  remaining  : " << numGlobalRows - numGlobalAggregated << "\n"
                                << "  aggregates : " << numGlobalAggs - numGlobalAggsPrev << " (phase), " << numGlobalAggs << " (total)" << std::endl;
        numGlobalAggregatedPrev = numGlobalAggregated;
        numGlobalAggsPrev       = numGlobalAggs;
      }
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(numNonAggregatedNodes, Exceptions::RuntimeError, "MueLu::UncoupledAggregationFactory::Build: Leftover nodes found! Error!");

  aggregates->AggregatesCrossProcessors(false);
  aggregates->ComputeAggregateSizes(true /*forceRecompute*/);

  Set(currentLevel, "Aggregates", aggregates);

  if (pL.get<bool>("aggregation: compute aggregate qualities")) {
    RCP<Xpetra::MultiVector<DefaultScalar, LO, GO, Node>> aggQualities = Get<RCP<Xpetra::MultiVector<DefaultScalar, LO, GO, Node>>>(currentLevel, "AggregateQualities");
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void UncoupledAggregationFactory<LocalOrdinal, GlobalOrdinal, Node>::
    DoGraphColoring(Level& currentLevel,
                    const std::string& aggAlgo,
                    const bool deterministic,
                    const RCP<const LWGraph_kokkos> graph,
                    RCP<Aggregates> aggregates) const {
  SubFactoryMonitor sfm(*this, "Algo \"Graph Coloring\"", currentLevel);

  // LBV on Sept 06 2019: the note below is a little worrisome,
  // can we guarantee that MueLu is never used on a non-symmetric
  // graph?
  // note: just using colinds_view in place of scalar_view_t type
  // (it won't be used at all by symbolic SPGEMM)
  using graph_t      = typename LWGraph_kokkos::local_graph_type;
  using KernelHandle = KokkosKernels::Experimental::
      KokkosKernelsHandle<typename graph_t::row_map_type::value_type,
                          typename graph_t::entries_type::value_type,
                          typename graph_t::entries_type::value_type,
                          typename graph_t::device_type::execution_space,
                          typename graph_t::device_type::memory_space,
                          typename graph_t::device_type::memory_space>;
  KernelHandle kh;
  // leave gc algorithm choice as the default
  kh.create_distance2_graph_coloring_handle();

  // get the distance-2 graph coloring handle
  auto coloringHandle = kh.get_distance2_graph_coloring_handle();

  const LO numRows = graph->GetNodeNumVertices();

  // Set the distance-2 graph coloring algorithm to use.
  // Options:
  //     COLORING_D2_DEFAULT        - Let the kernel handle pick the variation
  //     COLORING_D2_SERIAL         - Use the legacy serial-only implementation
  //     COLORING_D2_VB             - Use the parallel vertex based direct method
  //     COLORING_D2_VB_BIT         - Same as VB but using the bitvector forbidden array
  //     COLORING_D2_VB_BIT_EF      - Add experimental edge-filtering to VB_BIT
  //     COLORING_D2_NB_BIT         - Net-based coloring (generally the fastest)
  if (deterministic) {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_D2_SERIAL);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: serial" << std::endl;
  } else if (aggAlgo == "serial") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_D2_SERIAL);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: serial" << std::endl;
  } else if (aggAlgo == "default") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_D2_DEFAULT);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: default" << std::endl;
  } else if (aggAlgo == "vertex based") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_D2_VB);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based" << std::endl;
  } else if (aggAlgo == "vertex based bit set") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_D2_VB_BIT);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based bit set" << std::endl;
  } else if (aggAlgo == "edge filtering") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_D2_VB_BIT_EF);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: edge filtering" << std::endl;
  } else if (aggAlgo == "net based bit set") {
    coloringHandle->set_algorithm(KokkosGraph::COLORING_D2_NB_BIT);
    if (IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: net based bit set" << std::endl;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unrecognized distance 2 coloring algorithm, valid options are: serial, default, matrix squared, vertex based, vertex based bit set, edge filtering")
  }

  // Create device views for graph rowptrs/colinds
  typename graph_t::row_map_type aRowptrs = graph->getRowPtrs();
  typename graph_t::entries_type aColinds = graph->getEntries();

  // run d2 graph coloring
  // graph is symmetric so row map/entries and col map/entries are the same
  {
    SubFactoryMonitor sfm2(*this, "Algo \"Graph Coloring\": KokkosGraph Call", currentLevel);  // CMS HACK
    KokkosGraph::Experimental::graph_color_distance2(&kh, numRows, aRowptrs, aColinds);
  }

  // extract the colors and store them in the aggregates
  aggregates->SetGraphColors(coloringHandle->get_vertex_colors());
  aggregates->SetGraphNumColors(static_cast<LO>(coloringHandle->get_num_colors()));

  // clean up coloring handle
  kh.destroy_distance2_graph_coloring_handle();
}

}  // namespace MueLu

#endif /* MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_ */
