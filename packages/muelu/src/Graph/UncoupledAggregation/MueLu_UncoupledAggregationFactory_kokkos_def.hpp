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
#ifndef MUELU_UNCOUPLEDAGGREGATIONFACTORY_KOKKOS_DEF_HPP_
#define MUELU_UNCOUPLEDAGGREGATIONFACTORY_KOKKOS_DEF_HPP_

#ifdef HAVE_MUELU_KOKKOS_REFACTOR

#include <climits>

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_UncoupledAggregationFactory_kokkos_decl.hpp"

#include "MueLu_OnePtAggregationAlgorithm_kokkos.hpp"
#include "MueLu_PreserveDirichletAggregationAlgorithm_kokkos.hpp"
#include "MueLu_IsolatedNodeAggregationAlgorithm_kokkos.hpp"

#include "MueLu_AggregationPhase1Algorithm_kokkos.hpp"
#include "MueLu_AggregationPhase2aAlgorithm_kokkos.hpp"
#include "MueLu_AggregationPhase2bAlgorithm_kokkos.hpp"
#include "MueLu_AggregationPhase3Algorithm_kokkos.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Utilities.hpp" // for sum_all and similar stuff...

#include "KokkosGraph_Distance2ColorHandle.hpp"
#include "KokkosGraph_Distance2Color.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  UncoupledAggregationFactory_kokkos<LocalOrdinal, GlobalOrdinal, Node>::UncoupledAggregationFactory_kokkos()
  : bDefinitionPhase_(true)
  { }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> UncoupledAggregationFactory_kokkos<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    // Aggregation parameters (used in aggregation algorithms)
    // TODO introduce local member function for each aggregation algorithm such that each aggregation algorithm can define its own parameters

    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: max agg size");
    SET_VALID_ENTRY("aggregation: min agg size");
    SET_VALID_ENTRY("aggregation: max selected neighbors");
    SET_VALID_ENTRY("aggregation: ordering");
    validParamList->getEntry("aggregation: ordering").setValidator(
      rcp(new validatorType(Teuchos::tuple<std::string>("natural", "graph", "random"), "aggregation: ordering")));
    SET_VALID_ENTRY("aggregation: enable phase 1");
    SET_VALID_ENTRY("aggregation: phase 1 algorithm");
    SET_VALID_ENTRY("aggregation: deterministic");
    SET_VALID_ENTRY("aggregation: enable phase 2a");
    SET_VALID_ENTRY("aggregation: enable phase 2b");
    SET_VALID_ENTRY("aggregation: enable phase 3");
    SET_VALID_ENTRY("aggregation: preserve Dirichlet points");
    SET_VALID_ENTRY("aggregation: allow user-specified singletons");
#undef  SET_VALID_ENTRY

    // general variables needed in AggregationFactory
    validParamList->set< RCP<const FactoryBase> >("Graph",       null, "Generating factory of the graph");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");

    // special variables necessary for OnePtAggregationAlgorithm
    validParamList->set< std::string >           ("OnePt aggregate map name",         "", "Name of input map for single node aggregates. (default='')");
    validParamList->set< std::string >           ("OnePt aggregate map factory",      "", "Generating factory of (DOF) map for single node aggregates.");
    //validParamList->set< RCP<const FactoryBase> >("OnePt aggregate map factory",    NoFactory::getRCP(), "Generating factory of (DOF) map for single node aggregates.");

    return validParamList;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void UncoupledAggregationFactory_kokkos<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
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
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void UncoupledAggregationFactory_kokkos<LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    ParameterList pL = GetParameterList();
    bDefinitionPhase_ = false;  // definition phase is finished, now all aggregation algorithm information is fixed

    if (pL.get<int>("aggregation: max agg size") == -1)
      pL.set("aggregation: max agg size", INT_MAX);

    // define aggregation algorithms
    RCP<const FactoryBase> graphFact = GetFactory("Graph");

    // TODO Can we keep different aggregation algorithms over more Build calls?
    algos_.clear();
    algos_.push_back(rcp(new PreserveDirichletAggregationAlgorithm_kokkos(graphFact)));
    if (pL.get<bool>("aggregation: allow user-specified singletons") == true)   algos_.push_back(rcp(new OnePtAggregationAlgorithm_kokkos             (graphFact)));
    if (pL.get<bool>("aggregation: enable phase 1" )                 == true)   algos_.push_back(rcp(new AggregationPhase1Algorithm_kokkos            (graphFact)));
    if (pL.get<bool>("aggregation: enable phase 2a")                 == true)   algos_.push_back(rcp(new AggregationPhase2aAlgorithm_kokkos           (graphFact)));
    if (pL.get<bool>("aggregation: enable phase 2b")                 == true)   algos_.push_back(rcp(new AggregationPhase2bAlgorithm_kokkos           (graphFact)));
    if (pL.get<bool>("aggregation: enable phase 3" )                 == true)   algos_.push_back(rcp(new AggregationPhase3Algorithm_kokkos            (graphFact)));

    std::string mapOnePtName = pL.get<std::string>("OnePt aggregate map name");
    RCP<Map> OnePtMap = Teuchos::null;
    if (mapOnePtName.length()) {
      std::string mapOnePtFactName = pL.get<std::string>("OnePt aggregate map factory");
      if (mapOnePtFactName == "" || mapOnePtFactName == "NoFactory") {
        OnePtMap = currentLevel.Get<RCP<Map> >(mapOnePtName, NoFactory::get());
      } else {
        RCP<const FactoryBase> mapOnePtFact = GetFactory(mapOnePtFactName);
        OnePtMap = currentLevel.Get<RCP<Map> >(mapOnePtName, mapOnePtFact.get());
      }
    }

    RCP<const LWGraph_kokkos> graph = Get< RCP<LWGraph_kokkos> >(currentLevel, "Graph");

    // Build
    RCP<Aggregates_kokkos> aggregates = rcp(new Aggregates_kokkos(*graph));
    aggregates->setObjectLabel("UC");

    const LO numRows = graph->GetNodeNumVertices();

    // construct aggStat information
    Kokkos::View<unsigned*, typename LWGraph_kokkos::memory_space> aggStat("aggregation status",
                                                                           numRows);
    Kokkos::deep_copy(aggStat, READY);

    // LBV on Sept 06 2019: re-commenting out the dirichlet boundary map
    // even if the map is correctly extracted from the graph, aggStat is
    // now a Kokkos::View<unsinged*, memory_space> and filling it will
    // require a parallel_for or to copy it to the Host which is not really
    // good from a performance point of view.
    // If dirichletBoundaryMap was an actual Xpetra::Map, one could call
    // getLocalMap to have a Kokkos::View on the appropriate memory_space
    // instead of an ArrayRCP.
    // // TODO
    // //ArrayRCP<const bool> dirichletBoundaryMap = graph->GetBoundaryNodeMap();
    // ArrayRCP<const bool> dirichletBoundaryMap;

    // if (dirichletBoundaryMap != Teuchos::null)
    //   for (LO i = 0; i < numRows; i++)
    //     if (dirichletBoundaryMap[i] == true)
    //       aggStat[i] = BOUNDARY;

    LO nDofsPerNode = Get<LO>(currentLevel, "DofsPerNode");
    GO indexBase = graph->GetDomainMap()->getIndexBase();
    if (OnePtMap != Teuchos::null) {
      typename Kokkos::View<unsigned*,typename LWGraph_kokkos::memory_space>::HostMirror aggStatHost
        = Kokkos::create_mirror_view(aggStat);
      Kokkos::deep_copy(aggStatHost, aggStat);

      for (LO i = 0; i < numRows; i++) {
        // reconstruct global row id (FIXME only works for contiguous maps)
        GO grid = (graph->GetDomainMap()->getGlobalElement(i)-indexBase) * nDofsPerNode + indexBase;

        for (LO kr = 0; kr < nDofsPerNode; kr++)
          if (OnePtMap->isNodeGlobalElement(grid + kr))
            aggStatHost(i) = ONEPT;
      }

      Kokkos::deep_copy(aggStat, aggStatHost);
    }


    const RCP<const Teuchos::Comm<int> > comm = graph->GetComm();
    GO numGlobalRows = 0;
    if (IsPrint(Statistics1))
      MueLu_sumAll(comm, as<GO>(numRows), numGlobalRows);

    std::string algoParamName = "aggregation: phase 1 algorithm";
    if(pL.isParameter(algoParamName) && pL.get<std::string>(algoParamName) == "Distance2") {
      SubFactoryMonitor sfm(*this, "Algo \"Graph Coloring\"", currentLevel);

      // LBV on Sept 06 2019: the note below is a little worrisome,
      // can we guarantee that MueLu is never used on a non-symmetric
      // graph?
      // note: just using colinds_view in place of scalar_view_t type
      // (it won't be used at all by symbolic SPGEMM)
      using graph_t = typename LWGraph_kokkos::local_graph_type;
      using KernelHandle = KokkosKernels::Experimental::
        KokkosKernelsHandle<typename graph_t::row_map_type::value_type,
                            typename graph_t::entries_type::value_type,
                            typename graph_t::entries_type::value_type,
                            typename graph_t::device_type::execution_space,
                            typename graph_t::device_type::memory_space,
                            typename graph_t::device_type::memory_space>;
      KernelHandle kh;
      //leave gc algorithm choice as the default
      kh.create_distance2_graph_coloring_handle();

      // get the distance-2 graph coloring handle
      auto coloringHandle = kh.get_distance2_graph_coloring_handle();

      // Set the distance-2 graph coloring algorithm to use.
      // Options:
      //     COLORING_D2_DEFAULT        - Let the kernel handle pick the variation
      //     COLORING_D2_SERIAL         - Use the legacy serial-only implementation
      //     COLORING_D2_MATRIX_SQUARED - Use the SPGEMM + D1GC method
      //     COLORING_D2_SPGEMM         - Same as MATRIX_SQUARED
      //     COLORING_D2_VB             - Use the parallel vertex based direct method
      //     COLORING_D2_VB_BIT         - Same as VB but using the bitvector forbidden array
      //     COLORING_D2_VB_BIT_EF      - Add experimental edge-filtering to VB_BIT
      coloringHandle->set_algorithm( KokkosGraph::COLORING_D2_SERIAL );

      //Create device views for graph rowptrs/colinds
      typename graph_t::row_map_type aRowptrs = graph->getRowPtrs();
      typename graph_t::entries_type aColinds = graph->getEntries();

      //run d2 graph coloring
      //graph is symmetric so row map/entries and col map/entries are the same
      KokkosGraph::Experimental::graph_compute_distance2_color(&kh, numRows, numRows,
                                                               aRowptrs, aColinds,
                                                               aRowptrs, aColinds);

      // extract the colors and store them in the aggregates
      aggregates->SetGraphColors(coloringHandle->get_vertex_colors());
      aggregates->SetGraphNumColors(static_cast<LO>(coloringHandle->get_num_colors()));

      //clean up coloring handle
      kh.destroy_distance2_graph_coloring_handle();
    }

    LO numNonAggregatedNodes = numRows;
    GO numGlobalAggregatedPrev = 0, numGlobalAggsPrev = 0;
    for (size_t a = 0; a < algos_.size(); a++) {
      std::string phase = algos_[a]->description();
      SubFactoryMonitor sfm(*this, "Algo \"" + phase + "\"", currentLevel);

      int oldRank = algos_[a]->SetProcRankVerbose(this->GetProcRankVerbose());
      algos_[a]->BuildAggregates(pL, *graph, *aggregates, aggStat, numNonAggregatedNodes);
      algos_[a]->SetProcRankVerbose(oldRank);

      if (IsPrint(Statistics1)) {
        GO numLocalAggregated = numRows - numNonAggregatedNodes, numGlobalAggregated = 0;
        GO numLocalAggs       = aggregates->GetNumAggregates(),  numGlobalAggs = 0;
        MueLu_sumAll(comm, numLocalAggregated, numGlobalAggregated);
        MueLu_sumAll(comm, numLocalAggs,       numGlobalAggs);

        double aggPercent = 100*as<double>(numGlobalAggregated)/as<double>(numGlobalRows);
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
                                << "  aggregates : " << numGlobalAggs-numGlobalAggsPrev << " (phase), " << numGlobalAggs << " (total)" << std::endl;
        numGlobalAggregatedPrev = numGlobalAggregated;
        numGlobalAggsPrev       = numGlobalAggs;
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(numNonAggregatedNodes, Exceptions::RuntimeError, "MueLu::UncoupledAggregationFactory::Build: Leftover nodes found! Error!");

    aggregates->AggregatesCrossProcessors(false);
    aggregates->ComputeAggregateSizes(true/*forceRecompute*/);

    Set(currentLevel, "Aggregates", aggregates);

    GetOStream(Statistics1) << aggregates->description() << std::endl;
  }

} //namespace MueLu

#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif /* MUELU_UNCOUPLEDAGGREGATIONFACTORY_DEF_HPP_ */
