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

#include "Kokkos_StaticCrsGraph.hpp"
#include "KokkosGraph_Distance2ColorHandle.hpp"
#include "KokkosGraph_Distance2Color.hpp"

#include <Xpetra_Matrix.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_AggregationPhase1Algorithm_kokkos.hpp"
#include "MueLu_AggregationPhase2aAlgorithm_kokkos.hpp"
#include "MueLu_AggregationPhase2bAlgorithm_kokkos.hpp"
#include "MueLu_AggregationPhase3Algorithm_kokkos.hpp"
#include "MueLu_AmalgamationInfo_kokkos.hpp"
#include "MueLu_AmalgamationFactory_kokkos.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_UncoupledAggregationFactory_kokkos.hpp"

//#include "MueLu_UseDefaultTypes.hpp"

namespace MueLuTests {

  // Little utility to generate uncoupled aggregates.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void gimmeUncoupledAggregates_kokkos(const Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                                       RCP<MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node> >& graph,
                                       Teuchos::RCP<MueLu::Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Node> >& aggregates,
                                       bool bPhase1 = true, bool bPhase2a = true, bool bPhase2b = true, bool bPhase3 = true) {
#   include "MueLu_UseShortNames.hpp"
    Level level;
    TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<AmalgamationFactory_kokkos> amalgFact = rcp(new AmalgamationFactory_kokkos());
    RCP<CoalesceDropFactory_kokkos> dropFact  = rcp(new CoalesceDropFactory_kokkos());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    level.Request("Graph", dropFact.get());
    level.Request(*dropFact);
    dropFact->Build(level);
    graph = level.Get<RCP<LWGraph_kokkos> >("Graph", dropFact.get());
    const LO numNodes = graph->GetNodeNumVertices();
    aggregates = rcp(new Aggregates_kokkos(*graph));
    aggregates->setObjectLabel("UC");

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
    coloringHandle->set_algorithm( KokkosGraph::COLORING_D2_SERIAL );

    //Create device views for graph rowptrs/colinds
    typename graph_t::row_map_type aRowptrs = graph->getRowPtrs();
    typename graph_t::entries_type aColinds = graph->getEntries();

    //run d2 graph coloring
    //graph is symmetric so row map/entries and col map/entries are the same
    KokkosGraph::Experimental::graph_compute_distance2_color(&kh, numNodes, numNodes,
                                                             aRowptrs, aColinds,
                                                             aRowptrs, aColinds);

    // extract the colors and store them in the aggregates
    aggregates->SetGraphColors(coloringHandle->get_vertex_colors());
    aggregates->SetGraphNumColors(static_cast<LO>(coloringHandle->get_num_colors()));

    LO numNonAggregatedNodes = 0;
    Kokkos::View<unsigned*, typename LWGraph_kokkos::memory_space> aggStat("aggStat", numNodes);
    Kokkos::deep_copy(aggStat, MueLu::READY);
    Teuchos::ParameterList params;
    params.set<int> ("aggregation: min agg size", 1);
    params.set<int> ("aggregation: max agg size", 3);
    params.set<bool>("aggregation: deterministic", false);

    if(bPhase1) {
      RCP<MueLu::AggregationAlgorithmBase_kokkos<LO,GO,NO> > phase1
        = rcp(new AggregationPhase1Algorithm_kokkos(dropFact));
      phase1->BuildAggregates(params, *graph, *aggregates, aggStat, numNonAggregatedNodes);
    }
    if(bPhase2a) {
      RCP<MueLu::AggregationAlgorithmBase_kokkos<LO,GO,NO> > phase2a
        = rcp(new AggregationPhase2aAlgorithm_kokkos(dropFact));
      phase2a->BuildAggregates(params, *graph, *aggregates, aggStat, numNonAggregatedNodes);
    }
    if(bPhase2b) {
      RCP<MueLu::AggregationAlgorithmBase_kokkos<LO,GO,NO> > phase2b
        = rcp(new AggregationPhase2bAlgorithm_kokkos(dropFact));
      phase2b->BuildAggregates(params, *graph, *aggregates, aggStat, numNonAggregatedNodes);
    }
    if(bPhase3) {
      RCP<MueLu::AggregationAlgorithmBase_kokkos<LO,GO,NO> > phase3
        = rcp(new AggregationPhase3Algorithm_kokkos(dropFact));
      phase3->BuildAggregates(params, *graph, *aggregates, aggStat, numNonAggregatedNodes);
    }
    aggregates->AggregatesCrossProcessors(false);
    aggregates->ComputeAggregateSizes(true/*forceRecompute*/);
    level.Release("Graph", dropFact.get());
  }

  // Little utility to generate uncoupled aggregates.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MueLu::Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Node> >
  gimmeUncoupledAggregates_kokkos(const Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                                  bool bPhase1 = true, bool bPhase2a = true, bool bPhase2b = true, bool bPhase3 = true) {
#   include "MueLu_UseShortNames.hpp"
    RCP<LWGraph_kokkos> graph;
    RCP<Aggregates_kokkos> aggregates;

    gimmeUncoupledAggregates_kokkos(A, graph, aggregates, bPhase1, bPhase2a, bPhase2b, bPhase3);

    return aggregates;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal checkAggregatesContiguous(MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node> graph,
                                         MueLu::Aggregates_kokkos<LocalOrdinal, GlobalOrdinal, Node> aggregates) {
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using Aggregates_kokkos = MueLu::Aggregates_kokkos<LO, GO, Node>;
    using LWGraph_kokkos    = MueLu::LWGraph_kokkos<LO, GO, Node>;
    using execution_space   = typename LWGraph_kokkos::execution_space;
    using memory_space      = typename LWGraph_kokkos::memory_space;

    const LO numNodes = graph.GetNodeNumVertices();
    auto vertex2AggId = aggregates.GetVertex2AggId()->template getLocalView<memory_space>();
    auto aggSizes     = aggregates.ComputeAggregateSizes(true);

    Kokkos::View<LO*, memory_space> discontiguousAggs("discontiguous aggregates",
                                                      aggregates.GetNumAggregates());
    Kokkos::parallel_for("Mark discontiguous aggregates",
                         Kokkos::RangePolicy<LO, execution_space>(0, numNodes),
                         KOKKOS_LAMBDA(const LO nodeIdx) {
                           const LO myAggId   = vertex2AggId(nodeIdx, 0);
                           // Check that the node is actually aggregated
                           if(myAggId == -1) {return;}
                           const LO myAggSize = aggSizes(myAggId);

                           if(myAggSize == 1) {
                             // Can't have a discontiguous singleton
                             return;
                           } else {
                             auto neighbors = graph.getNeighborVertices(nodeIdx);
                             for(LO neigh = 0; neigh < neighbors.length; ++neigh) {
                               const LO neighIdx   = neighbors(neigh);
                               const LO neighAggId = vertex2AggId(neighIdx, 0);
                               if((nodeIdx != neighIdx) && (neighAggId == myAggId)) {
                                 // This aggregate might be discontiguous
                                 // but at least not because of this node
                                 return;
                               }
                             }
                             discontiguousAggs(myAggId) = 1;
                           }
                         });

    LO numDiscontiguousAggregates = 0;
    Kokkos::parallel_reduce("Count discontiguous aggregates",
                            Kokkos::RangePolicy<LO, execution_space>(0, aggregates.GetNumAggregates()),
                            KOKKOS_LAMBDA(const LO aggIdx, LO& numDiscontiguous) {
                              if(discontiguousAggs(aggIdx) == 1) {
                                ++numDiscontiguous;
                              }
                            }, numDiscontiguousAggregates);

    return numDiscontiguousAggregates;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, JustUncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);

    RCP<Aggregates_kokkos> aggregates = gimmeUncoupledAggregates_kokkos(A);

    TEST_EQUALITY(aggregates != Teuchos::null,              true);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),  false);
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, JustDist2UncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    //TODO bmk: A lot of test code duplicated here from gimmeUncoupledAggregates
    //because it can't take a custom parameter list, add that as parameter?
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
    RCP<AmalgamationInfo_kokkos> amalgInfo;
    Level level;
    TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<AmalgamationFactory_kokkos> amalgFact = rcp(new AmalgamationFactory_kokkos());
    RCP<CoalesceDropFactory_kokkos> dropFact  = rcp(new CoalesceDropFactory_kokkos());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    // Setup aggregation factory (use default factory for graph)
    RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
    aggFact->SetFactory("Graph", dropFact);
    aggFact->SetParameter("aggregation: max agg size",           Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: min agg size",           Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(0));
    aggFact->SetParameter("aggregation: enable phase 1",         Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 2a",        Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 2b",        Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 3",         Teuchos::ParameterEntry(true));

    level.Request("Aggregates", aggFact.get());
    level.Request("UnAmalgamationInfo", amalgFact.get());

    level.Request(*aggFact);
    aggFact->Build(level);
    RCP<Aggregates_kokkos> aggregates = level.Get<RCP<Aggregates_kokkos> >("Aggregates",aggFact.get()); // fix me
    TEST_INEQUALITY(aggregates, Teuchos::null);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
    amalgInfo = level.Get<RCP<AmalgamationInfo_kokkos> >("UnAmalgamationInfo",amalgFact.get()); // fix me
    level.Release("UnAmalgamationInfo", amalgFact.get());
    level.Release("Aggregates", aggFact.get());
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, JustDist2DeterUncoupledAggregation, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
    //TODO bmk: A lot of test code duplicated here from gimmeUncoupledAggregates
    //because it can't take a custom parameter list, add that as parameter?
#   include <MueLu_UseShortNames.hpp>
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
    RCP<AmalgamationInfo_kokkos> amalgInfo;
    Level level;
    TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<AmalgamationFactory_kokkos> amalgFact = rcp(new AmalgamationFactory_kokkos());
    RCP<CoalesceDropFactory_kokkos> dropFact  = rcp(new CoalesceDropFactory_kokkos());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    // Setup aggregation factory (use default factory for graph)
    RCP<UncoupledAggregationFactory_kokkos> aggFact = rcp(new UncoupledAggregationFactory_kokkos());
    aggFact->SetFactory("Graph", dropFact);
    aggFact->SetParameter("aggregation: max agg size",           Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: min agg size",           Teuchos::ParameterEntry(3));
    aggFact->SetParameter("aggregation: max selected neighbors", Teuchos::ParameterEntry(0));
    aggFact->SetParameter("aggregation: deterministic",          Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 1",         Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 2a",        Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 2b",        Teuchos::ParameterEntry(true));
    aggFact->SetParameter("aggregation: enable phase 3",         Teuchos::ParameterEntry(true));

    level.Request("Aggregates", aggFact.get());
    level.Request("UnAmalgamationInfo", amalgFact.get());

    level.Request(*aggFact);
    aggFact->Build(level);
    RCP<Aggregates_kokkos> aggregates = level.Get<RCP<Aggregates_kokkos> >("Aggregates",aggFact.get()); // fix me
    TEST_INEQUALITY(aggregates, Teuchos::null);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);
    amalgInfo = level.Get<RCP<AmalgamationInfo_kokkos> >("UnAmalgamationInfo",amalgFact.get()); // fix me
    level.Release("UnAmalgamationInfo", amalgFact.get());
    level.Release("Aggregates", aggFact.get());
  }


// A test that creates discontiguous aggregates to make sure the detection algorithm works well
/// pretty much testing the test...
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, DiscontiguousAggregates, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    using memory_space = typename Aggregates_kokkos::device_type::memory_space;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers_kokkos::Parameters::getDefaultComm();
    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(3*comm->getSize());
    Level level;
    TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(level);
    level.Set("A", A);

    RCP<Aggregates_kokkos> aggregates;
    RCP<LWGraph_kokkos> graph;

    RCP<AmalgamationFactory_kokkos> amalgFact = rcp(new AmalgamationFactory_kokkos());
    RCP<CoalesceDropFactory_kokkos> dropFact  = rcp(new CoalesceDropFactory_kokkos());
    dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

    level.Request("Graph", dropFact.get());
    level.Request(*dropFact);
    dropFact->Build(level);
    graph = level.Get<RCP<LWGraph_kokkos> >("Graph", dropFact.get());
    RCP<const Map> importMap = graph->GetImportMap();
    const LO numNodes = graph->GetNodeNumVertices();
    aggregates = rcp(new Aggregates_kokkos(*graph));
    aggregates->setObjectLabel("UC");

    Kokkos::View<unsigned*, typename LWGraph_kokkos::memory_space> aggStat("aggStat", numNodes);
    Kokkos::deep_copy(aggStat, MueLu::READY);

    // Performing fake aggregates to generate a discontiguous aggregate
    Kokkos::View<LO**, Kokkos::LayoutLeft, memory_space> vertex2AggId = aggregates->GetVertex2AggId()->template getLocalView<memory_space>();
    Kokkos::View<LO**, Kokkos::LayoutLeft, memory_space> procWinner   = aggregates->GetProcWinner()->template getLocalView<memory_space>();

    typename Kokkos::View<LO**, Kokkos::LayoutLeft, memory_space>::HostMirror vertex2AggId_h
      = Kokkos::create_mirror_view(vertex2AggId);
    Kokkos::deep_copy(vertex2AggId_h, vertex2AggId);
    vertex2AggId_h(0, 0) = 0;
    vertex2AggId_h(1, 0) = 1;
    vertex2AggId_h(2, 0) = 0;
    Kokkos::deep_copy(vertex2AggId, vertex2AggId_h);

    typename Kokkos::View<LO**, Kokkos::LayoutLeft, memory_space>::HostMirror procWinner_h
      = Kokkos::create_mirror_view(procWinner);
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
  } //UncoupledPhase1


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, UncoupledPhase1, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<AmalgamationInfo_kokkos> amalgInfo;

    RCP<Aggregates_kokkos> aggregates;
    RCP<LWGraph_kokkos> graph;
    gimmeUncoupledAggregates_kokkos(A, graph, aggregates, true, false, false, false);
    const GO numAggs = aggregates->GetNumAggregates();

    TEST_INEQUALITY(numAggs, 0);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(), false);

    typename Aggregates_kokkos::aggregates_sizes_type::const_type aggSizes
      = aggregates->ComputeAggregateSizes(true);
    typename Aggregates_kokkos::aggregates_sizes_type::const_type::HostMirror aggSizes_h =
      Kokkos::create_mirror_view(aggSizes);
    Kokkos::deep_copy(aggSizes_h, aggSizes);

    LO numBadAggregates = 0;
    Kokkos::parallel_reduce("Checking aggregates sizes",
                            Kokkos::RangePolicy<LO, typename Aggregates_kokkos::execution_space>(0, aggSizes_h.extent(0)),
                            KOKKOS_LAMBDA(const LO aggIdx, LO& lNumBadAggregates) {
                              if ((aggSizes_h(aggIdx) < 1) || (3 < aggSizes_h(aggIdx))) {
                                lNumBadAggregates += 1;
                              }
                            }, numBadAggregates);
    TEST_EQUALITY(numBadAggregates, 0);

    const LO numDiscontiguous = checkAggregatesContiguous(*graph, *aggregates);
    TEST_EQUALITY(numDiscontiguous, 0);
  } //UncoupledPhase1

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, UncoupledPhase2, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);

    RCP<LWGraph_kokkos> graph;
    RCP<Aggregates_kokkos> aggregates;
    gimmeUncoupledAggregates_kokkos(A, graph, aggregates, false, true, true, false);
    GO numAggs = aggregates->GetNumAggregates();

    TEST_INEQUALITY(numAggs, 0);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);

    const LO numDiscontiguous = checkAggregatesContiguous(*graph, *aggregates);
    TEST_EQUALITY(numDiscontiguous, 0);
  } //UncoupledPhase2

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Aggregates_kokkos, UncoupledPhase3, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
    RCP<AmalgamationInfo_kokkos> amalgInfo;

    RCP<Aggregates_kokkos> aggregates = gimmeUncoupledAggregates_kokkos(A, false, false, false, true);
    GO numAggs = aggregates->GetNumAggregates();

    TEST_INEQUALITY(numAggs, 0);
    TEST_EQUALITY(aggregates->AggregatesCrossProcessors(),false);

    typename Aggregates_kokkos::aggregates_sizes_type::const_type aggSizes
      = aggregates->ComputeAggregateSizes(true);
    typename Aggregates_kokkos::aggregates_sizes_type::const_type::HostMirror aggSizes_h =
      Kokkos::create_mirror_view(aggSizes);
    Kokkos::deep_copy(aggSizes_h, aggSizes);

    LO numBadAggregates = 0;
    Kokkos::parallel_reduce("Checking aggregates sizes",
                            Kokkos::RangePolicy<LO, typename Aggregates_kokkos::execution_space>(0, aggSizes_h.extent(0)),
                            KOKKOS_LAMBDA(const LO aggIdx, LO& lNumBadAggregates) {
                              if ((aggSizes_h(aggIdx) < 1) || (5 < aggSizes_h(aggIdx))) {
                                lNumBadAggregates += 1;
                              }
                            }, numBadAggregates);
    TEST_EQUALITY(numBadAggregates, 0);

  } //UncoupledPhase3

#define MUELU_ETI_GROUP(SC,LO,GO,NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, JustUncoupledAggregation, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, JustDist2UncoupledAggregation, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, JustDist2DeterUncoupledAggregation, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, DiscontiguousAggregates, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, UncoupledPhase1, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, UncoupledPhase2, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Aggregates_kokkos, UncoupledPhase3, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

} // namespace MueLuTests
