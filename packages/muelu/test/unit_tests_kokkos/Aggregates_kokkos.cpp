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
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_VerboseObject.hpp>

#include "Kokkos_StaticCrsGraph.hpp"
#include "KokkosGraph_Distance2ColorHandle.hpp"
#include "KokkosGraph_Distance2Color.hpp"

#include <Xpetra_Matrix.hpp>
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AggregationPhase1Algorithm_kokkos.hpp"
#include "MueLu_AggregationPhase2aAlgorithm_kokkos.hpp"
#include "MueLu_AggregationPhase2bAlgorithm_kokkos.hpp"
#include "MueLu_AggregationPhase3Algorithm_kokkos.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos.hpp"

//#include "MueLu_UseDefaultTypes.hpp"

namespace MueLuTests {

// Little utility to generate uncoupled aggregates.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void gimmeUncoupledAggregates(const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                              RCP<MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>>& graph,
                              Teuchos::RCP<MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node>>& aggregates,
                              bool bPhase1 = true, bool bPhase2a = true, bool bPhase2b = true, bool bPhase3 = true) {
#include "MueLu_UseShortNames.hpp"
  Level level;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  level.Request("Graph", dropFact.get());
  level.Request(*dropFact);
  dropFact->Build(level);
  graph             = level.Get<RCP<LWGraph_kokkos>>("Graph", dropFact.get());
  const LO numNodes = graph->GetNodeNumVertices();
  aggregates        = rcp(new Aggregates(*graph));
  aggregates->setObjectLabel("UC");

  using graph_t      = typename LWGraph_kokkos::local_graph_type;
  using device_type  = typename graph_t::device_type;
  using KernelHandle = KokkosKernels::Experimental::
      KokkosKernelsHandle<typename graph_t::row_map_type::value_type,
                          typename graph_t::entries_type::value_type,
                          typename graph_t::entries_type::value_type,
                          typename graph_t::device_type::execution_space,
                          typename graph_t::device_type::memory_space,
                          typename graph_t::device_type::memory_space>;
  KernelHandle kh;
  // Leave gc algorithm choice as the default:
  // COLORING_D2_SERIAL for Serial execspace, and COLORING_D2_NB_BIT otherwise.
  kh.create_distance2_graph_coloring_handle();

  // get the distance-2 graph coloring handle
  auto coloringHandle = kh.get_distance2_graph_coloring_handle();

  // Create device views for graph rowptrs/colinds
  typename graph_t::row_map_type aRowptrs = graph->getLocalLWGraph().getRowPtrs();
  typename graph_t::entries_type aColinds = graph->getLocalLWGraph().getEntries();

  // run d2 graph coloring
  // graph is symmetric so row map/entries and col map/entries are the same
  KokkosGraph::Experimental::graph_color_distance2(&kh, numNodes, aRowptrs, aColinds);

  // extract the colors and store them in the aggregates
  aggregates->SetGraphColors(coloringHandle->get_vertex_colors());
  aggregates->SetGraphNumColors(static_cast<LO>(coloringHandle->get_num_colors()));

  LO numNonAggregatedNodes = 0;
  Kokkos::View<unsigned*, device_type> aggStat("aggStat", numNodes);
  Kokkos::deep_copy(aggStat, MueLu::READY);
  Teuchos::ParameterList params;
  params.set<int>("aggregation: min agg size", 1);
  params.set<int>("aggregation: max agg size", 3);
  params.set<bool>("aggregation: deterministic", false);

  params.set<bool>("aggregation: match ML phase2a", true);
  params.set<bool>("aggregation: error on nodes with no on-rank neighbors", false);
  params.set<bool>("aggregation: phase3 avoid singletons", false);

  if (bPhase1) {
    RCP<MueLu::AggregationAlgorithmBase_kokkos<LO, GO, NO>> phase1 = rcp(new AggregationPhase1Algorithm_kokkos(dropFact));
    phase1->BuildAggregates(params, *graph, *aggregates, aggStat, numNonAggregatedNodes);
  }
  if (bPhase2a) {
    RCP<MueLu::AggregationAlgorithmBase_kokkos<LO, GO, NO>> phase2a = rcp(new AggregationPhase2aAlgorithm_kokkos(dropFact));
    phase2a->BuildAggregates(params, *graph, *aggregates, aggStat, numNonAggregatedNodes);
  }
  if (bPhase2b) {
    RCP<MueLu::AggregationAlgorithmBase_kokkos<LO, GO, NO>> phase2b = rcp(new AggregationPhase2bAlgorithm_kokkos(dropFact));
    phase2b->BuildAggregates(params, *graph, *aggregates, aggStat, numNonAggregatedNodes);
  }
  if (bPhase3) {
    RCP<MueLu::AggregationAlgorithmBase_kokkos<LO, GO, NO>> phase3 = rcp(new AggregationPhase3Algorithm_kokkos(dropFact));
    phase3->BuildAggregates(params, *graph, *aggregates, aggStat, numNonAggregatedNodes);
  }
  aggregates->AggregatesCrossProcessors(false);
  aggregates->ComputeAggregateSizes(true /*forceRecompute*/);
  level.Release("Graph", dropFact.get());
}

// Little utility to generate uncoupled aggregates.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node>>
gimmeUncoupledAggregates(const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                         bool bPhase1 = true, bool bPhase2a = true, bool bPhase2b = true, bool bPhase3 = true) {
#include "MueLu_UseShortNames.hpp"
  RCP<LWGraph_kokkos> graph;
  RCP<Aggregates> aggregates;

  gimmeUncoupledAggregates(A, graph, aggregates, bPhase1, bPhase2a, bPhase2b, bPhase3);

  return aggregates;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal checkAggregatesContiguous(MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node> graph,
                                       MueLu::Aggregates<LocalOrdinal, GlobalOrdinal, Node> aggregates) {
  using LO              = LocalOrdinal;
  using GO              = GlobalOrdinal;
  using Aggregates      = MueLu::Aggregates<LO, GO, Node>;
  using LWGraph_kokkos  = MueLu::LWGraph_kokkos<LO, GO, Node>;
  using execution_space = typename LWGraph_kokkos::execution_space;
  using memory_space    = typename LWGraph_kokkos::memory_space;
  using device_type     = typename LWGraph_kokkos::device_type;

  const LO numNodes = graph.GetNodeNumVertices();
  auto vertex2AggId = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadOnly);
  auto aggSizes     = aggregates.ComputeAggregateSizes(true);

  auto lclLWGraph = graph.getLocalLWGraph();

  Kokkos::View<LO*, device_type> discontiguousAggs("discontiguous aggregates",
                                                   aggregates.GetNumAggregates());
  Kokkos::parallel_for(
      "Mark discontiguous aggregates",
      Kokkos::RangePolicy<LO, execution_space>(0, numNodes),
      KOKKOS_LAMBDA(const LO nodeIdx) {
        const LO myAggId = vertex2AggId(nodeIdx, 0);
        // Check that the node is actually aggregated
        if (myAggId == -1) {
          return;
        }
        const LO myAggSize = aggSizes(myAggId);

        if (myAggSize == 1) {
          // Can't have a discontiguous singleton
          return;
        } else {
          auto neighbors = lclLWGraph.getNeighborVertices(nodeIdx);
          for (LO neigh = 0; neigh < neighbors.length; ++neigh) {
            const LO neighIdx   = neighbors(neigh);
            const LO neighAggId = vertex2AggId(neighIdx, 0);
            if ((nodeIdx != neighIdx) && (neighAggId == myAggId)) {
              // This aggregate might be discontiguous
              // but at least not because of this node
              return;
            }
          }
          discontiguousAggs(myAggId) = 1;
        }
      });

  LO numDiscontiguousAggregates = 0;
  Kokkos::parallel_reduce(
      "Count discontiguous aggregates",
      Kokkos::RangePolicy<LO, execution_space>(0, aggregates.GetNumAggregates()),
      KOKKOS_LAMBDA(const LO aggIdx, LO& numDiscontiguous) {
        if (discontiguousAggs(aggIdx) == 1) {
          ++numDiscontiguous;
        }
      },
      numDiscontiguousAggregates);

  return numDiscontiguousAggregates;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, JustUncoupledAggregationFactory, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  // Setup aggregation factory (use default factory for graph)
  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetOrdering("graph");
  TEST_EQUALITY(aggFact->GetOrdering() == "graph", true);
  aggFact->SetOrdering("natural");
  TEST_EQUALITY(aggFact->GetOrdering() == "natural", true);
  aggFact->SetOrdering("random");
  TEST_EQUALITY(aggFact->GetOrdering() == "random", true);

  aggFact->SetMaxNeighAlreadySelected(12);
  TEST_EQUALITY(aggFact->GetMaxNeighAlreadySelected() == 12, true);

  aggFact->SetMaxNeighAlreadySelected(0);
  TEST_EQUALITY(aggFact->GetMaxNeighAlreadySelected() == 0, true);

  aggFact->SetMinNodesPerAggregate(0);
  TEST_EQUALITY(aggFact->GetMinNodesPerAggregate() == 0, true);

  aggFact->SetMinNodesPerAggregate(3);
  TEST_EQUALITY(aggFact->GetMinNodesPerAggregate() == 3, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, JustUncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);

  RCP<Aggregates> aggregates = gimmeUncoupledAggregates(A);

  TEST_EQUALITY(aggregates != Teuchos::null, true);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, JustDist2UncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // TODO bmk: A lot of test code duplicated here from gimmeUncoupledAggregates
  // because it can't take a custom parameter list, add that as parameter?
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);
  aggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));
  aggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(3));
  aggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(0));
  aggFact->SetParameter("aggregation: preserve Dirichlet points", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 1", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 3", Teuchos::ParameterEntry(true));

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);
  aggFact->Build(level);
  RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());  // fix me
  TEST_INEQUALITY(aggregates, Teuchos::null);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
  amalgInfo = level.Get<RCP<AmalgamationInfo>>("UnAmalgamationInfo", amalgFact.get());  // fix me
  level.Release("UnAmalgamationInfo", amalgFact.get());
  level.Release("Aggregates", aggFact.get());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, JustDist2PreserveUncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // TODO bmk: A lot of test code duplicated here from gimmeUncoupledAggregates
  // because it can't take a custom parameter list, add that as parameter?
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);
  aggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));
  aggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(3));
  aggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(0));
  aggFact->SetParameter("aggregation: enable phase 1", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 3", Teuchos::ParameterEntry(true));

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);
  aggFact->Build(level);
  RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());  // fix me
  TEST_INEQUALITY(aggregates, Teuchos::null);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
  amalgInfo = level.Get<RCP<AmalgamationInfo>>("UnAmalgamationInfo", amalgFact.get());  // fix me
  level.Release("UnAmalgamationInfo", amalgFact.get());
  level.Release("Aggregates", aggFact.get());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, JustOnePt2UncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // TODO bmk: A lot of test code duplicated here from gimmeUncoupledAggregates
  // because it can't take a custom parameter list, add that as parameter?
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);
  aggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));
  aggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(3));
  aggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(0));
  aggFact->SetParameter("aggregation: allow user-specified singletons", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 1", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 3", Teuchos::ParameterEntry(true));

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);
  aggFact->Build(level);
  RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());  // fix me
  TEST_INEQUALITY(aggregates, Teuchos::null);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
  amalgInfo = level.Get<RCP<AmalgamationInfo>>("UnAmalgamationInfo", amalgFact.get());  // fix me
  level.Release("UnAmalgamationInfo", amalgFact.get());
  level.Release("Aggregates", aggFact.get());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, JustDist2DeterUncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // TODO bmk: A lot of test code duplicated here from gimmeUncoupledAggregates
  // because it can't take a custom parameter list, add that as parameter?
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);
  aggFact->SetParameter("aggregation: max agg size", Teuchos::ParameterEntry(3));
  aggFact->SetParameter("aggregation: min agg size", Teuchos::ParameterEntry(3));
  aggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(0));
  aggFact->SetParameter("aggregation: deterministic", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 1", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 2a", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 2b", Teuchos::ParameterEntry(true));
  aggFact->SetParameter("aggregation: enable phase 3", Teuchos::ParameterEntry(true));

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);

  // Set some node to ONEPT state

  aggFact->Build(level);
  RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());  // fix me
  TEST_INEQUALITY(aggregates, Teuchos::null);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
  amalgInfo = level.Get<RCP<AmalgamationInfo>>("UnAmalgamationInfo", amalgFact.get());  // fix me
  level.Release("UnAmalgamationInfo", amalgFact.get());
  level.Release("Aggregates", aggFact.get());
}

// A test that creates discontiguous aggregates to make sure the detection algorithm works well
/// pretty much testing the test...
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, DiscontiguousAggregates, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  using device_type  = typename Aggregates::device_type;
  using memory_space = typename device_type::memory_space;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers_kokkos::Parameters::getDefaultComm();
  RCP<Matrix> A                      = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(3 * comm->getSize());
  Level level;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<Aggregates> aggregates;
  RCP<LWGraph_kokkos> graph;

  RCP<AmalgamationFactory> amalgFact       = rcp(new AmalgamationFactory());
  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  level.Request("Graph", dropFact.get());
  level.Request(*dropFact);
  dropFact->Build(level);
  graph                    = level.Get<RCP<LWGraph_kokkos>>("Graph", dropFact.get());
  RCP<const Map> importMap = graph->GetImportMap();
  const LO numNodes        = graph->GetNodeNumVertices();
  aggregates               = rcp(new Aggregates(*graph));
  aggregates->setObjectLabel("UC");

  Kokkos::View<unsigned*, device_type> aggStat("aggStat", numNodes);
  Kokkos::deep_copy(aggStat, MueLu::READY);

  // Performing fake aggregates to generate a discontiguous aggregate
  Kokkos::View<LO**, Kokkos::LayoutLeft, device_type> vertex2AggId = aggregates->GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  Kokkos::View<LO**, Kokkos::LayoutLeft, device_type> procWinner   = aggregates->GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);

  typename Kokkos::View<LO**, Kokkos::LayoutLeft, device_type>::HostMirror vertex2AggId_h = Kokkos::create_mirror_view(vertex2AggId);
  Kokkos::deep_copy(vertex2AggId_h, vertex2AggId);
  vertex2AggId_h(0, 0) = 0;
  vertex2AggId_h(1, 0) = 1;
  vertex2AggId_h(2, 0) = 0;
  Kokkos::deep_copy(vertex2AggId, vertex2AggId_h);

  typename Kokkos::View<LO**, Kokkos::LayoutLeft, device_type>::HostMirror procWinner_h = Kokkos::create_mirror_view(procWinner);
  Kokkos::deep_copy(procWinner_h, procWinner);
  procWinner_h(0, 0) = comm->getRank();
  procWinner_h(1, 0) = comm->getRank();
  procWinner_h(2, 0) = comm->getRank();
  Kokkos::deep_copy(procWinner, procWinner_h);

  aggregates->SetNumAggregates(2);
  aggregates->AggregatesCrossProcessors(false);
  aggregates->ComputeAggregateSizes(true);

  const LO numDiscontiguous = checkAggregatesContiguous(*graph, *aggregates);
  TEST_EQUALITY(numDiscontiguous, 1);
}  // UncoupledPhase1

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, UncoupledPhase1, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  RCP<AmalgamationInfo> amalgInfo;

  RCP<Aggregates> aggregates;
  RCP<LWGraph_kokkos> graph;
  gimmeUncoupledAggregates(A, graph, aggregates, true, false, false, false);
  const GO numAggs = aggregates->GetNumAggregates();

  TEST_INEQUALITY(numAggs, 0);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);

  typename Aggregates::aggregates_sizes_type::const_type aggSizes = aggregates->ComputeAggregateSizes(true);

  LO numBadAggregates = 0;
  Kokkos::parallel_reduce(
      "Checking aggregates sizes",
      Kokkos::RangePolicy<LO, typename Aggregates::execution_space>(0, aggSizes.extent(0)),
      KOKKOS_LAMBDA(const LO aggIdx, LO& lNumBadAggregates) {
        if ((aggSizes(aggIdx) < 1) || (3 < aggSizes(aggIdx))) {
          lNumBadAggregates += 1;
        }
      },
      numBadAggregates);
  TEST_EQUALITY(numBadAggregates, 0);

  const LO numDiscontiguous = checkAggregatesContiguous(*graph, *aggregates);
  TEST_EQUALITY(numDiscontiguous, 0);
}  // UncoupledPhase1

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, UncoupledPhase2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);

  RCP<LWGraph_kokkos> graph;
  RCP<Aggregates> aggregates;
  gimmeUncoupledAggregates(A, graph, aggregates, false, true, true, false);
  GO numAggs = aggregates->GetNumAggregates();

  TEST_INEQUALITY(numAggs, 0);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);

  const LO numDiscontiguous = checkAggregatesContiguous(*graph, *aggregates);
  TEST_EQUALITY(numDiscontiguous, 0);
}  // UncoupledPhase2

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, UncoupledPhase3, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  RCP<AmalgamationInfo> amalgInfo;

  RCP<Aggregates> aggregates = gimmeUncoupledAggregates(A, false, false, false, true);
  GO numAggs                 = aggregates->GetNumAggregates();

  TEST_INEQUALITY(numAggs, 0);
  TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);

  typename Aggregates::aggregates_sizes_type::const_type aggSizes = aggregates->ComputeAggregateSizes(true);

  LO numBadAggregates = 0;
  Kokkos::parallel_reduce(
      "Checking aggregates sizes",
      Kokkos::RangePolicy<LO, typename Aggregates::execution_space>(0, aggSizes.extent(0)),
      KOKKOS_LAMBDA(const LO aggIdx, LO& lNumBadAggregates) {
        if ((aggSizes(aggIdx) < 1) || (5 < aggSizes(aggIdx))) {
          lNumBadAggregates += 1;
        }
      },
      numBadAggregates);
  TEST_EQUALITY(numBadAggregates, 0);

  // Check ComputeNodesInAggregate
  typename Aggregates::LO_view aggPtr, aggNodes, unaggregated;
  aggregates->ComputeNodesInAggregate(aggPtr, aggNodes, unaggregated);
  TEST_EQUALITY(aggPtr.extent_int(0), numAggs + 1);
  // TEST_EQUALITY(unaggregated.extent_int(0), 0); // 1 unaggregated node in the MPI_4 case

  // test aggPtr(i)+aggSizes(i)=aggPtr(i+1)
  typename Aggregates::LO_view::HostMirror aggPtr_h                 = Kokkos::create_mirror_view(aggPtr);
  typename Aggregates::aggregates_sizes_type::HostMirror aggSizes_h = Kokkos::create_mirror_view(aggSizes);
  Kokkos::deep_copy(aggPtr_h, aggPtr);
  Kokkos::deep_copy(aggSizes_h, aggSizes);
  for (LO i = 0; i < aggSizes_h.extent_int(0); ++i)
    TEST_EQUALITY(aggPtr_h(i) + aggSizes_h(i), aggPtr_h(i + 1));

}  // UncoupledPhase3

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, AllowDroppingToCreateAdditionalDirichletRows, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << "Test option that allows dropping during aggregation to create new Dirichlet rows" << std::endl;

  typedef typename Teuchos::ScalarTraits<Scalar> TST;
  using device_type = typename Aggregates::device_type;

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  A->resumeFill();

  // Create one row on every processor with small off-diagonal entries that will be dropped with
  // an appropriately chosen threshold.  Avoid domain boundaries.
  LocalOrdinal localRowToModify = 1;
  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const Scalar> values;
  A->getLocalRowView(localRowToModify, indices, values);
  Array<Scalar> newvalues(values.size(), TST::zero());
  for (int j = 0; j < indices.size(); j++) {
    if (indices[j] == localRowToModify)
      newvalues[j] = values[j];  // keep diagonal unmodified
    else
      newvalues[j] = -TST::eps();
  }
  A->replaceLocalValues(localRowToModify, indices, newvalues);
  A->fillComplete();

  // Dropping connections will not lead to the creation of new Dirichlet rows.
  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory_kokkos> dropFact = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetParameter("aggregation: dropping may create Dirichlet", Teuchos::ParameterEntry(false));
  dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(Teuchos::as<double>(-100 * TST::eps())));
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  level.Request("Graph", dropFact.get());

  // Setup aggregation factory (use default factory for graph)
  RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);
  aggFact->Build(level);
  RCP<Aggregates> aggregates                                  = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
  RCP<LWGraph_kokkos> graph                                   = level.Get<RCP<LWGraph_kokkos>>("Graph", dropFact.get());
  Kokkos::View<const bool*, device_type> dirichletBoundaryMap = graph->getLocalLWGraph().GetBoundaryNodeMap();
  int numDirichletRows                                        = 0;
  using execution_space                                       = typename LWGraph_kokkos::execution_space;
  LO numRows                                                  = graph->GetNodeNumVertices();
  Kokkos::parallel_reduce(
      "Count Dirichlet rows",
      Kokkos::RangePolicy<LO, execution_space>(0, numRows),
      KOKKOS_LAMBDA(const LO rowIdx, LO& numDirichlet) {
        if (dirichletBoundaryMap[rowIdx]) {
          ++numDirichlet;
        }
      },
      numDirichletRows);
  TEST_EQUALITY(numDirichletRows, 0);

  Array<LO> aggPtr;
  Array<LO> aggNodes;
  Array<LO> unaggregated;

  // Repeat with "aggregation: dropping may create Dirichlet" = TRUE, i.e.,
  // dropping connections may lead to the creation of new Dirichlet rows.
  // The second row should be flagged as Dirichlet because all off-diagonal entries are dropped.
  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetParameter("aggregation: dropping may create Dirichlet", Teuchos::ParameterEntry(true));
  dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(Teuchos::as<double>(-100 * TST::eps())));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);
  level.Request("Graph", dropFact.get());

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);
  aggFact->Build(level);
  aggregates           = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
  graph                = level.Get<RCP<LWGraph_kokkos>>("Graph", dropFact.get());
  dirichletBoundaryMap = graph->getLocalLWGraph().GetBoundaryNodeMap();
  numDirichletRows     = 0;
  Kokkos::parallel_reduce(
      "Count Dirichlet rows",
      Kokkos::RangePolicy<LO, execution_space>(0, numRows),
      KOKKOS_LAMBDA(const LO rowIdx, LO& numDirichlet) {
        if (dirichletBoundaryMap[rowIdx]) {
          ++numDirichlet;
        }
      },
      numDirichletRows);
  TEST_EQUALITY(numDirichletRows, 1);

}  // AllowDroppingToCreateAdditionalDirichletRows

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, GreedyDirichlet, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;
  out << "running greedy Dirichlet" << std::endl;

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  //    RCP<Matrix> A = TestHelpers::TestFactory<SC, LO, GO, NO>::Build1DPoisson(30);
  // Make a Matrix with multiple degrees of freedom per node
  GlobalOrdinal nx = 8, ny = 8;

  // Describes the initial layout of matrix rows across processors.
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", nx);
  galeriList.set("ny", ny);
  RCP<const Teuchos::Comm<int>> comm = TestHelpers_kokkos::Parameters::getDefaultComm();
  RCP<const Map> map                 = Galeri::Xpetra::CreateMap<LocalOrdinal, GlobalOrdinal, Node>(TestHelpers_kokkos::Parameters::getLib(), "Cartesian2D", comm, galeriList);

  map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, 2);  // expand map for 2 DOFs per node

  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
      Galeri::Xpetra::BuildProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, CrsMatrixWrap, MultiVector>("Elasticity2D", map, galeriList);
  RCP<Matrix> A = Pr->BuildMatrix();
  A->SetFixedBlockSize(2);

  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const Scalar> values;

  // Create a dirichlet boundary row.
  LocalOrdinal localRowToZero = 5;  // Corresponds to a Dof on local graph node 2

  A->resumeFill();
  A->getLocalRowView(localRowToZero, indices, values);
  Array<Scalar> newvalues(values.size(), TST::zero());
  for (int j = 0; j < indices.size(); j++)
    // keep diagonal
    if (indices[j] == localRowToZero) newvalues[j] = values[j];
  A->replaceLocalValues(localRowToZero, indices, newvalues);

  A->fillComplete();

  ArrayRCP<const bool> drows = Utilities::DetectDirichletRows(*A);
  TEST_EQUALITY(drows[localRowToZero], true);
  TEST_EQUALITY(drows[localRowToZero - 1], false);

  RCP<AmalgamationInfo> amalgInfo;
  Level level;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(level);
  level.Set("A", A);

  RCP<CoalesceDropFactory_kokkos> dropFact;
  RCP<AmalgamationFactory> amalgFact;
  RCP<UncoupledAggregationFactory_kokkos> aggFact;

  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetParameter("aggregation: greedy Dirichlet", Teuchos::ParameterEntry(false));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);

  level.Request("Aggregates", aggFact.get());
  level.Request("UnAmalgamationInfo", amalgFact.get());

  level.Request(*aggFact);
  aggFact->Build(level);
  RCP<Aggregates> aggregates = level.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
  Array<LO> aggPtr;
  Array<LO> aggNodes;
  Array<LO> unaggregated;

  typename Aggregates::aggregates_sizes_type::const_type aggSizes = aggregates->ComputeAggregateSizes(true);

  auto vertex2AggId = aggregates->GetVertex2AggId()->getHostLocalView(Xpetra::Access::ReadOnly);
  for (auto i = 0; i < (nx / 2 * ny / 2); i++) {
    TEST_EQUALITY(vertex2AggId(i, 0) != MUELU_UNAGGREGATED, true);  // check that all nodes are aggregated
  }

  // Repeat with greedy Dirichlet
  Level levelGreedyAndNoPreserve;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(levelGreedyAndNoPreserve);
  levelGreedyAndNoPreserve.Set("A", A);
  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetParameter("aggregation: greedy Dirichlet", Teuchos::ParameterEntry(true));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);

  levelGreedyAndNoPreserve.Request("Aggregates", aggFact.get());
  levelGreedyAndNoPreserve.Request("UnAmalgamationInfo", amalgFact.get());

  levelGreedyAndNoPreserve.Request(*aggFact);
  aggFact->Build(levelGreedyAndNoPreserve);
  aggregates = levelGreedyAndNoPreserve.Get<RCP<Aggregates>>("Aggregates", aggFact.get());
  aggFact->SetParameter("aggregation: preserve Dirichlet points", Teuchos::ParameterEntry(false));

  typename Aggregates::aggregates_sizes_type::const_type aggSizesGreedy = aggregates->ComputeAggregateSizes(true);
  // There should be the same number of aggregates in the greedy dirichlet case as in the standard case
  // though note the aggregates will be smaller in the greedy dirirchlet case.
  // This may not be true for all problem setups, but is true for the setup in this test case.
  TEST_EQUALITY(aggSizesGreedy.extent(0) == aggSizes.extent(0), true)

  vertex2AggId = aggregates->GetVertex2AggId()->getHostLocalView(Xpetra::Access::ReadOnly);
  TEST_EQUALITY(vertex2AggId(2, 0) == MUELU_UNAGGREGATED, true);  // check that the node with the Dof flagged as dirichlet is unaggregated

  // Repeat with greedy Dirichlet and preserve Dirichlet points
  Level levelGreedyAndPreserve;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(levelGreedyAndPreserve);
  levelGreedyAndPreserve.Set("A", A);

  amalgFact = rcp(new AmalgamationFactory());
  dropFact  = rcp(new CoalesceDropFactory_kokkos());
  dropFact->SetParameter("aggregation: greedy Dirichlet", Teuchos::ParameterEntry(true));
  dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

  // Setup aggregation factory (use default factory for graph)
  aggFact = rcp(new UncoupledAggregationFactory_kokkos());
  aggFact->SetFactory("Graph", dropFact);
  aggFact->SetParameter("aggregation: preserve Dirichlet points", Teuchos::ParameterEntry(true));

  levelGreedyAndPreserve.Request("Aggregates", aggFact.get());
  levelGreedyAndPreserve.Request("UnAmalgamationInfo", amalgFact.get());

  levelGreedyAndPreserve.Request(*aggFact);
  aggFact->Build(levelGreedyAndPreserve);
  aggregates = levelGreedyAndPreserve.Get<RCP<Aggregates>>("Aggregates", aggFact.get());

  typename Aggregates::aggregates_sizes_type::const_type aggSizesGreedyPreserve = aggregates->ComputeAggregateSizes(true);
  // check that dirichlet aggs are preserved
  // there should be more aggregates in the dirichlet preserved case than in the standard case.
  // This will always be true.
  TEST_EQUALITY(aggSizesGreedyPreserve.extent(0) > aggSizesGreedy.extent(0), true)

  vertex2AggId = aggregates->GetVertex2AggId()->getHostLocalView(Xpetra::Access::ReadOnly);
  for (auto i = 0; i < (nx / 2 * ny / 2); i++) {
    TEST_EQUALITY(vertex2AggId(i, 0) != MUELU_UNAGGREGATED, true);  // check that all nodes are aggregated
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, JustUncoupledAggregationFactory, SC, LO, GO, NO)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, JustUncoupledAggregation, SC, LO, GO, NO)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, JustDist2UncoupledAggregation, SC, LO, GO, NO)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, JustDist2PreserveUncoupledAggregation, SC, LO, GO, NO)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, JustDist2DeterUncoupledAggregation, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, JustOnePt2UncoupledAggregation, SC, LO, GO, NO)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, DiscontiguousAggregates, SC, LO, GO, NO)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, UncoupledPhase1, SC, LO, GO, NO)                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, UncoupledPhase2, SC, LO, GO, NO)                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, UncoupledPhase3, SC, LO, GO, NO)                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, AllowDroppingToCreateAdditionalDirichletRows, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, GreedyDirichlet, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
