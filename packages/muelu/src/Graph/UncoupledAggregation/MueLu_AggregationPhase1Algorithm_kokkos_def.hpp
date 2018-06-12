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
#ifndef MUELU_AGGREGATIONPHASE1ALGORITHM_KOKKOS_DEF_HPP
#define MUELU_AGGREGATIONPHASE1ALGORITHM_KOKKOS_DEF_HPP

#ifdef HAVE_MUELU_KOKKOS_REFACTOR

#include <queue>
#include <vector>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase1Algorithm_kokkos_decl.hpp"

#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Monitor.hpp"

#include "KokkosGraph_graph_color.hpp"
#include <Kokkos_ScatterView.hpp>

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationPhase1Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregates(const ParameterList& params, const LWGraph_kokkos& graph,
                  Aggregates_kokkos& aggregates, Kokkos::View<unsigned*, typename MueLu::
                  LWGraph_kokkos<LO,GO,Node>::local_graph_type::device_type::
                  memory_space>& aggStatView, LO& numNonAggregatedNodes,
                  Kokkos::View<LO*, typename MueLu::LWGraph_kokkos<LO, GO, Node>::
                  local_graph_type::device_type::memory_space>& colorsDevice, LO& numColors) const {
    Monitor m(*this, "BuildAggregates");

    std::string orderingStr     = params.get<std::string>("aggregation: ordering");
    int maxNeighAlreadySelected = params.get<int>        ("aggregation: max selected neighbors");
    int minNodesPerAggregate    = params.get<int>        ("aggregation: min agg size");
    int maxNodesPerAggregate    = params.get<int>        ("aggregation: max agg size");

    bool colorPermuted          = params.get<bool>       ("aggregation: permute nodes by color");

    Algorithm algorithm         = Algorithm::Serial;
    std::string algoParamName   = "aggregation: coloring algorithm";
    if(params.isParameter(algoParamName))
    {
      algorithm = algorithmFromName(params.get<std::string>("aggregation: coloring algorithm"));
    }

    TEUCHOS_TEST_FOR_EXCEPTION(maxNodesPerAggregate < minNodesPerAggregate, Exceptions::RuntimeError,
                               "MueLu::UncoupledAggregationAlgorithm::BuildAggregates: minNodesPerAggregate must be smaller or equal to MaxNodePerAggregate!");

    //Distance-2 gives less control than serial uncoupled phase 1
    //no custom row reordering because would require making deep copy of local matrix entries and permuting it
    //can only enforce max aggregate size
    if (algorithm == Algorithm::Serial)
    {
      const LO numRows = graph.GetNodeNumVertices();
      std::vector<unsigned> aggStat(numRows);
      typedef typename MueLu::LWGraph_kokkos<LO, GO, Node>::local_graph_type graph_t;
      typedef typename graph_t::device_type::memory_space memory_space;
      typename Kokkos::View<unsigned*, memory_space>::HostMirror h_aggStatView =
        Kokkos::create_mirror_view (aggStatView);
      Kokkos::deep_copy(h_aggStatView, aggStatView);
      // Note: LBV 20/03/2018
      // This last part could be done by passing a raw pointer from h_aggStatView to aggStat.
      // However I am hoping that we will depart with the serial case all together soon, so this
      // piece of code should not impact us much.
      for(LO i = 0; i < numRows; ++i) {aggStat[i] = h_aggStatView(i);}
      BuildAggregatesSerial(graph, aggregates, aggStat, numNonAggregatedNodes, minNodesPerAggregate,
                            maxNodesPerAggregate, maxNeighAlreadySelected, orderingStr);
      for(LO i = 0; i < numRows; ++i) {h_aggStatView(i) = aggStat[i];}
      Kokkos::deep_copy(aggStatView, h_aggStatView);
    }
    else if (algorithm == Algorithm::Distance2)
    {
      BuildAggregatesDistance2(graph, aggregates, aggStatView, numNonAggregatedNodes,
                               maxNodesPerAggregate, colorsDevice, numColors, colorPermuted);
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationPhase1Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregatesSerial(const LWGraph_kokkos& graph, Aggregates_kokkos& aggregates,
      std::vector<unsigned>& aggStat, LO& numNonAggregatedNodes,
      LO minNodesPerAggregate, LO maxNodesPerAggregate,
      LO maxNeighAlreadySelected, std::string& orderingStr) const
  {
    enum {
      O_NATURAL,
      O_RANDOM,
      O_GRAPH
    } ordering;

    ordering = O_NATURAL; // initialize variable (fix CID 143665)
    if (orderingStr == "natural") ordering = O_NATURAL;
    if (orderingStr == "random" ) ordering = O_RANDOM;
    if (orderingStr == "graph"  ) ordering = O_GRAPH;

    const LO  numRows = graph.GetNodeNumVertices();
    const int myRank  = graph.GetComm()->getRank();

    ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);

    LO numLocalAggregates = aggregates.GetNumAggregates();

    ArrayRCP<LO> randomVector;
    if (ordering == O_RANDOM) {
      randomVector = arcp<LO>(numRows);
      for (LO i = 0; i < numRows; i++)
        randomVector[i] = i;
      RandomReorder(randomVector);
    }

    int              aggIndex = -1;
    size_t           aggSize  =  0;
    std::vector<int> aggList(graph.getNodeMaxNumRowEntries());

    std::queue<LO>   graphOrderQueue;

    // Main loop over all local rows of graph(A)
    for (LO i = 0; i < numRows; i++) {
      // Step 1: pick the next node to aggregate
      LO rootCandidate = 0;
      if      (ordering == O_NATURAL) rootCandidate = i;
      else if (ordering == O_RANDOM)  rootCandidate = randomVector[i];
      else if (ordering == O_GRAPH) {

        if (graphOrderQueue.size() == 0) {
          // Current queue is empty for "graph" ordering, populate with one READY node
          for (LO jnode = 0; jnode < numRows; jnode++)
            if (aggStat[jnode] == READY) {
              graphOrderQueue.push(jnode);
              break;
            }
        }
        if (graphOrderQueue.size() == 0) {
          // There are no more ready nodes, end the phase
          break;
        }
        rootCandidate = graphOrderQueue.front();   // take next node from graph ordering queue
        graphOrderQueue.pop();                     // delete this node in list
      }

      if (aggStat[rootCandidate] != READY)
        continue;

      // Step 2: build tentative aggregate
      aggSize = 0;
      aggList[aggSize++] = rootCandidate;

      auto neighOfINode = graph.getNeighborVertices(rootCandidate);

      // If the number of neighbors is less than the minimum number of nodes
      // per aggregate, we know this is not going to be a valid root, and we
      // may skip it, but only for "natural" and "random" (for "graph" we still
      // need to fetch the list of local neighbors to continue)
      if ((ordering == O_NATURAL || ordering == O_RANDOM) &&
          as<int>(neighOfINode.length) < minNodesPerAggregate) {
        continue;
      }

      LO numAggregatedNeighbours = 0;

      for (int j = 0; j < as<int>(neighOfINode.length); j++) {
        LO neigh = neighOfINode(j);

        if (neigh != rootCandidate && graph.isLocalNeighborVertex(neigh)) {

          if (aggStat[neigh] == READY || aggStat[neigh] == NOTSEL) {
            // If aggregate size does not exceed max size, add node to the
            // tentative aggregate
            // NOTE: We do not exit the loop over all neighbours since we have
            // still to count all aggregated neighbour nodes for the
            // aggregation criteria
            // NOTE: We check here for the maximum aggregation size. If we
            // would do it below with all the other check too big aggregates
            // would not be accepted at all.
            if (aggSize < as<size_t>(maxNodesPerAggregate))
              aggList[aggSize++] = neigh;

          } else {
            numAggregatedNeighbours++;
          }
        }
      }

      // Step 3: check if tentative aggregate is acceptable
      if ((numAggregatedNeighbours <= maxNeighAlreadySelected) &&           // too many connections to other aggregates
          (aggSize                 >= as<size_t>(minNodesPerAggregate))) {  // too few nodes in the tentative aggregate
        // Accept new aggregate
        // rootCandidate becomes the root of the newly formed aggregate
        aggregates.SetIsRoot(rootCandidate);
        aggIndex = numLocalAggregates++;

        for (size_t k = 0; k < aggSize; k++) {
          aggStat     [aggList[k]] = AGGREGATED;
          vertex2AggId[aggList[k]] = aggIndex;
          procWinner  [aggList[k]] = myRank;
        }

        numNonAggregatedNodes -= aggSize;

      } else {
        // Aggregate is not accepted
        aggStat[rootCandidate] = NOTSEL;

        // Need this for the "graph" ordering below
        // The original candidate is always aggList[0]
        aggSize = 1;
      }

      if (ordering == O_GRAPH) {
        // Add candidates to the list of nodes
        // NOTE: the code have slightly different meanings depending on context:
        //  - if aggregate was accepted, we add neighbors of neighbors of the original candidate
        //  - if aggregate was not accepted, we add neighbors of the original candidate
        for (size_t k = 0; k < aggSize; k++) {
          auto neighOfJNode = graph.getNeighborVertices(aggList[k]);

          for (int j = 0; j < as<int>(neighOfJNode.length); j++) {
            LO neigh = neighOfJNode(j);

            if (graph.isLocalNeighborVertex(neigh) && aggStat[neigh] == READY)
              graphOrderQueue.push(neigh);
          }
        }
      }
    }

    // Reset all NOTSEL vertices to READY
    // This simplifies other algorithms
    for (LO i = 0; i < numRows; i++)
      if (aggStat[i] == NOTSEL)
        aggStat[i] = READY;

    // update aggregate object
    aggregates.SetNumAggregates(numLocalAggregates);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationPhase1Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregatesDistance2(const LWGraph_kokkos& graph, Aggregates_kokkos& aggregates,
                           Kokkos::View<unsigned*, typename MueLu::LWGraph_kokkos<LO, GO, Node>::
                           local_graph_type::device_type::memory_space>& aggStatView,
                           LO& numNonAggregatedNodes, LO maxAggSize,
                           Kokkos::View<LO*, typename MueLu::LWGraph_kokkos<LO, GO, Node>::
                           local_graph_type::device_type::memory_space>& colorsDevice,
                           LO& numColors, const bool colorPermuted) const
  {
    typedef typename MueLu::LWGraph_kokkos<LO, GO, Node>::local_graph_type graph_t;
    typedef typename graph_t::device_type::memory_space memory_space;
    typedef typename graph_t::device_type::execution_space execution_space;

    const LO  numRows = graph.GetNodeNumVertices();
    const int myRank  = graph.GetComm()->getRank();

    auto vertex2AggIdView = aggregates.GetVertex2AggId()->template getLocalView<memory_space>();
    auto procWinnerView = aggregates.GetProcWinner()    ->template getLocalView<memory_space>();
    auto colorsCardinality = aggregates.GetColorsCardinality();
    auto colorPermutation = aggregates.GetColorPermutation();

    LO numNonAggregatedNodesSerial = numNonAggregatedNodes;
    LO numLocalAggregatesSerial = aggregates.GetNumAggregates();
    LO numLocalAggregates = aggregates.GetNumAggregates();

    typename LWGraph_kokkos::local_graph_type::entries_type::non_const_type::
      HostMirror h_colors = Kokkos::create_mirror_view(colorsDevice);
    Kokkos::deep_copy(h_colors, colorsDevice);

    Kokkos::View<LO, memory_space> aggCount("aggCount");
    LO tmpNumLocalAggregates = 0;

    if(colorPermuted) {
      Kokkos::parallel_reduce("Aggregation Phase 1: initial reduction over color == 1", colorsCardinality(1),
                              KOKKOS_LAMBDA (const LO i, LO& lnumLocalAggregates) {
                                LO nodeIdx = colorPermutation(i);
                                if(aggStatView(nodeIdx) == READY) {
                                  const LO aggIdx = Kokkos::atomic_fetch_add (&aggCount(), 1);
                                  vertex2AggIdView(nodeIdx, 0) = aggIdx;
                                  aggStatView(nodeIdx) = AGGREGATED;
                                  ++lnumLocalAggregates;
                                  procWinnerView(nodeIdx, 0) = myRank;
                                }
                              }, tmpNumLocalAggregates);
    } else {
      Kokkos::parallel_reduce("Aggregation Phase 1: initial reduction over color == 1", numRows,
                              KOKKOS_LAMBDA (const LO i, LO& lnumLocalAggregates) {
                                if(colorsDevice(i) == 1 && aggStatView(i) == READY) {
                                  const LO aggIdx = Kokkos::atomic_fetch_add (&aggCount(), 1);
                                  vertex2AggIdView(i, 0) = aggIdx;
                                  aggStatView(i) = AGGREGATED;
                                  ++lnumLocalAggregates;
                                  procWinnerView(i, 0) = myRank;
                                }
                              }, tmpNumLocalAggregates);
    }
    numLocalAggregates += tmpNumLocalAggregates;

    // Compute the initial size of the aggregates.
    // Note lbv 12-21-17: I am pretty sure that the aggregates will always be of size 1
    //                    at this point so we could simplify the code below a lot if this
    //                    assumption is correct...
    Kokkos::View<LO*, memory_space> aggSizesView("aggSizes", numLocalAggregates);
    {
      // Here there is a possibility that two vertices assigned to two different threads contribute
      // to the same aggregate if somethings happened before phase 1?
      auto aggSizesScatterView = Kokkos::Experimental::create_scatter_view(aggSizesView);
      Kokkos::parallel_for("Aggregation Phase 1: compute initial aggregates size", numRows,
                           KOKKOS_LAMBDA (const LO i) {
                             auto aggSizesScatterViewAccess = aggSizesScatterView.access();
                             if(vertex2AggIdView(i, 0) >= 0)
                               aggSizesScatterViewAccess(vertex2AggIdView(i, 0)) += 1;
                           });
      Kokkos::Experimental::contribute(aggSizesView, aggSizesScatterView);
    }

    if(maxAggSize == std::numeric_limits<int>::max()) {
      // When aggSize is "unlimited" we do not need to check aggSize and atomically decrement
      // if we reach the maximum set by users.
      if(colorPermuted) {
        int colorRange = colorsCardinality(numColors) - colorsCardinality(1);
        int colorOffset = colorsCardinality(1);
        Kokkos::parallel_reduce("Aggregation Phase 1: main parallel_reduce over aggSizes", colorRange,
                                KOKKOS_LAMBDA (const size_t i, LO & lNumNonAggregatedNodes) {
                                  LO nodeIdx = colorPermutation(i + colorOffset);
                                  if(aggStatView(nodeIdx) == READY || aggStatView(nodeIdx) == NOTSEL) {
                                    // Get neighbors of vertex i and look for local, aggregated,
                                    // color 1 neighbor (valid root).
                                    auto neighbors = graph.getNeighborVertices(nodeIdx);
                                    for(LO j = 0; j < neighbors.length; ++j) {
                                      auto nei = neighbors.colidx(j);
                                      if(graph.isLocalNeighborVertex(nei) && colorsDevice(nei) == 1
                                         && aggStatView(nei) == AGGREGATED) {

                                        // This atomic guarentees that any other node trying to
                                        // join aggregate agg has the correct size.
                                        LO agg = vertex2AggIdView(nei, 0);
                                        Kokkos::atomic_add (&aggSizesView(agg), 1);
                                        //assign vertex i to aggregate with root j
                                        vertex2AggIdView(nodeIdx, 0) = agg;
                                        procWinnerView(nodeIdx, 0)   = myRank;
                                        aggStatView(nodeIdx)         = AGGREGATED;
                                        break;
                                      }
                                    }
                                  }
                                  if(aggStatView(nodeIdx) != AGGREGATED) {
                                    lNumNonAggregatedNodes++;
                                    if(aggStatView(nodeIdx) == NOTSEL) { aggStatView(nodeIdx) = READY; }
                                  }
                                }, numNonAggregatedNodes);
      } else {
        Kokkos::parallel_reduce("Aggregation Phase 1: main parallel_reduce over aggSizes", numRows,
                                KOKKOS_LAMBDA (const size_t i, LO & lNumNonAggregatedNodes) {
                                  if(colorsDevice(i) != 1
                                     && (aggStatView(i) == READY || aggStatView(i) == NOTSEL)) {
                                    // Get neighbors of vertex i and look for local, aggregated,
                                    // color 1 neighbor (valid root).
                                    auto neighbors = graph.getNeighborVertices(i);
                                    for(LO j = 0; j < neighbors.length; ++j) {
                                      auto nei = neighbors.colidx(j);
                                      if(graph.isLocalNeighborVertex(nei) && colorsDevice(nei) == 1
                                         && aggStatView(nei) == AGGREGATED) {

                                        // This atomic guarentees that any other node trying to
                                        // join aggregate agg has the correct size.
                                        LO agg = vertex2AggIdView(nei, 0);
                                        Kokkos::atomic_add (&aggSizesView(agg), 1);
                                        //assign vertex i to aggregate with root j
                                        vertex2AggIdView(i, 0) = agg;
                                        procWinnerView(i, 0)   = myRank;
                                        aggStatView(i)         = AGGREGATED;
                                        break;
                                      }
                                    }
                                  }
                                  if(aggStatView(i) != AGGREGATED) {
                                    lNumNonAggregatedNodes++;
                                    if(aggStatView(i) == NOTSEL) { aggStatView(i) = READY; }
                                  }
                                }, numNonAggregatedNodes);
      }
    } else {
      if(colorPermuted) {
        int colorRange = colorsCardinality(numColors) - colorsCardinality(1);
        int colorOffset = colorsCardinality(1);
        Kokkos::parallel_reduce("Aggregation Phase 1: main parallel_reduce over aggSizes", colorRange,
                                KOKKOS_LAMBDA (const size_t i, LO & lNumNonAggregatedNodes) {
                                  LO nodeIdx = colorPermutation(i + colorOffset);
                                  if((aggStatView(nodeIdx) == READY || aggStatView(nodeIdx) == NOTSEL)) {
                                    // Get neighbors of vertex i and look for local, aggregated,
                                    // color 1 neighbor (valid root).
                                    auto neighbors = graph.getNeighborVertices(nodeIdx);
                                    for(LO j = 0; j < neighbors.length; ++j) {
                                      auto nei = neighbors.colidx(j);
                                      if(graph.isLocalNeighborVertex(nei) && colorsDevice(nei) == 1
                                         && aggStatView(nei) == AGGREGATED) {

                                        // This atomic guarentees that any other node trying to
                                        // join aggregate agg has the correct size.
                                        LO agg = vertex2AggIdView(nei, 0);
                                        const LO aggSize =
                                          Kokkos::atomic_fetch_add (&aggSizesView(agg), 1);
                                        if(aggSize < maxAggSize) {
                                          //assign vertex i to aggregate with root j
                                          vertex2AggIdView(nodeIdx, 0) = agg;
                                          procWinnerView(nodeIdx, 0)   = myRank;
                                          aggStatView(nodeIdx)         = AGGREGATED;
                                          break;
                                        } else {
                                          // Decrement back the value of aggSizesView(agg)
                                          Kokkos::atomic_decrement(&aggSizesView(agg));
                                        }
                                      }
                                    }
                                  }
                                  if(aggStatView(nodeIdx) != AGGREGATED) {
                                    lNumNonAggregatedNodes++;
                                    if(aggStatView(nodeIdx) == NOTSEL) { aggStatView(nodeIdx) = READY; }
                                  }
                                }, numNonAggregatedNodes);
      } else {
        Kokkos::parallel_reduce("Aggregation Phase 1: main parallel_reduce over aggSizes", numRows,
                                KOKKOS_LAMBDA (const size_t i, LO & lNumNonAggregatedNodes) {
                                  if(colorsDevice(i) != 1
                                     && (aggStatView(i) == READY || aggStatView(i) == NOTSEL)) {
                                    // Get neighbors of vertex i and look for local, aggregated,
                                    // color 1 neighbor (valid root).
                                    auto neighbors = graph.getNeighborVertices(i);
                                    for(LO j = 0; j < neighbors.length; ++j) {
                                      auto nei = neighbors.colidx(j);
                                      if(graph.isLocalNeighborVertex(nei) && colorsDevice(nei) == 1
                                         && aggStatView(nei) == AGGREGATED) {

                                        // This atomic guarentees that any other node trying to
                                        // join aggregate agg has the correct size.
                                        LO agg = vertex2AggIdView(nei, 0);
                                        const LO aggSize =
                                          Kokkos::atomic_fetch_add (&aggSizesView(agg), 1);
                                        if(aggSize < maxAggSize) {
                                          //assign vertex i to aggregate with root j
                                          vertex2AggIdView(i, 0) = agg;
                                          procWinnerView(i, 0)   = myRank;
                                          aggStatView(i)         = AGGREGATED;
                                          break;
                                        } else {
                                          // Decrement back the value of aggSizesView(agg)
                                          Kokkos::atomic_decrement(&aggSizesView(agg));
                                        }
                                      }
                                    }
                                  }
                                  if(aggStatView(i) != AGGREGATED) {
                                    lNumNonAggregatedNodes++;
                                    if(aggStatView(i) == NOTSEL) { aggStatView(i) = READY; }
                                  }
                                }, numNonAggregatedNodes);
      }
    }

    // update aggregate object
    aggregates.SetNumAggregates(numLocalAggregates);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationPhase1Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::RandomReorder(ArrayRCP<LO> list) const {
    //TODO: replace int
    int n = list.size();
    for(int i = 0; i < n-1; i++)
      std::swap(list[i], list[RandomOrdinal(i,n-1)]);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  int AggregationPhase1Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::RandomOrdinal(int min, int max) const {
    return min + as<int>((max-min+1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
  }

} // end namespace

#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_AGGREGATIONPHASE1ALGORITHM_KOKKOS_DEF_HPP

