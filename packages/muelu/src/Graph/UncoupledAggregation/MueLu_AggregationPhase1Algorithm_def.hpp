// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_AGGREGATIONPHASE1ALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONPHASE1ALGORITHM_DEF_HPP_

#include <queue>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase1Algorithm_decl.hpp"

#include "MueLu_LWGraph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

#include "Kokkos_Sort.hpp"
#include <Kokkos_ScatterView.hpp>

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregatesNonKokkos(const ParameterList& params, const LWGraph& graph, Aggregates& aggregates, typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatHostType& aggStat,
                             LO& numNonAggregatedNodes) const {
  Monitor m(*this, "BuildAggregatesNonKokkos");

  std::string orderingStr     = params.get<std::string>("aggregation: ordering");
  int maxNeighAlreadySelected = params.get<int>("aggregation: max selected neighbors");
  int minNodesPerAggregate    = params.get<int>("aggregation: min agg size");
  int maxNodesPerAggregate    = params.get<int>("aggregation: max agg size");
  bool matchMLBehavior        = params.get<bool>("aggregation: match ML phase1");

  TEUCHOS_TEST_FOR_EXCEPTION(maxNodesPerAggregate < minNodesPerAggregate, Exceptions::RuntimeError,
                             "MueLu::UncoupledAggregationAlgorithm::BuildAggregatesNonKokkos: minNodesPerAggregate must be smaller or equal to MaxNodePerAggregate!");

  enum {
    O_NATURAL,
    O_RANDOM,
    O_GRAPH
  } ordering;
  ordering = O_NATURAL;  // initialize variable (fix CID 143665)
  if (orderingStr == "natural") ordering = O_NATURAL;
  if (orderingStr == "random") ordering = O_RANDOM;
  if (orderingStr == "graph") ordering = O_GRAPH;

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  LO numLocalAggregates = aggregates.GetNumAggregates();

  ArrayRCP<LO> randomVector;
  if (ordering == O_RANDOM) {
    randomVector = arcp<LO>(numRows);
    for (LO i = 0; i < numRows; i++)
      randomVector[i] = i;
    RandomReorder(randomVector);
  }

  int aggIndex   = -1;
  size_t aggSize = 0;
  std::vector<int> aggList(graph.getLocalMaxNumRowEntries());

  std::queue<LO> graphOrderQueue;

  // Main loop over all local rows of graph(A)
  for (LO i = 0; i < numRows; i++) {
    // Step 1: pick the next node to aggregate
    LO rootCandidate = 0;
    if (ordering == O_NATURAL)
      rootCandidate = i;
    else if (ordering == O_RANDOM)
      rootCandidate = randomVector[i];
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
      rootCandidate = graphOrderQueue.front();  // take next node from graph ordering queue
      graphOrderQueue.pop();                    // delete this node in list
    }

    if (aggStat[rootCandidate] != READY)
      continue;

    // Step 2: build tentative aggregate
    aggSize            = 0;
    aggList[aggSize++] = rootCandidate;

    auto neighOfINode = graph.getNeighborVertices(rootCandidate);

    // If the number of neighbors is less than the minimum number of nodes
    // per aggregate, we know this is not going to be a valid root, and we
    // may skip it, but only for "natural" and "random" (for "graph" we still
    // need to fetch the list of local neighbors to continue)
    if ((ordering == O_NATURAL || ordering == O_RANDOM) &&
        neighOfINode.length < minNodesPerAggregate) {
      continue;
    }

    LO numAggregatedNeighbours = 0;

    for (int j = 0; j < neighOfINode.length; j++) {
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

        } else if (!matchMLBehavior || aggStat[neigh] != IGNORED) {
          // NOTE: ML checks against BOUNDARY here, but boundary nodes are flagged as IGNORED by
          // the time we get to Phase 1, so we check IGNORED instead
          numAggregatedNeighbours++;
        }
      }
    }

    // Step 3: check if tentative aggregate is acceptable
    if ((numAggregatedNeighbours <= maxNeighAlreadySelected) &&  // too many connections to other aggregates
        (aggSize >= as<size_t>(minNodesPerAggregate))) {         // too few nodes in the tentative aggregate
      // Accept new aggregate
      // rootCandidate becomes the root of the newly formed aggregate
      aggregates.SetIsRoot(rootCandidate);
      aggIndex = numLocalAggregates++;

      for (size_t k = 0; k < aggSize; k++) {
        aggStat[aggList[k]]      = AGGREGATED;
        vertex2AggId[aggList[k]] = aggIndex;
        procWinner[aggList[k]]   = myRank;
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

        for (int j = 0; j < neighOfJNode.length; j++) {
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
void AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node>::RandomReorder(ArrayRCP<LO> list) const {
  // TODO: replace int
  int n = list.size();
  for (int i = 0; i < n - 1; i++)
    std::swap(list[i], list[RandomOrdinal(i, n - 1)]);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
int AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node>::RandomOrdinal(int min, int max) const {
  return min + as<int>((max - min + 1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregates(const Teuchos::ParameterList& params,
                    const LWGraph_kokkos& graph,
                    Aggregates& aggregates,
                    typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                    LO& numNonAggregatedNodes) const {
  int minNodesPerAggregate = params.get<int>("aggregation: min agg size");
  int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");

  TEUCHOS_TEST_FOR_EXCEPTION(maxNodesPerAggregate < minNodesPerAggregate,
                             Exceptions::RuntimeError,
                             "MueLu::UncoupledAggregationAlgorithm::BuildAggregates: minNodesPerAggregate must be smaller or equal to MaxNodePerAggregate!");

  // Distance-2 gives less control than serial uncoupled phase 1
  // no custom row reordering because would require making deep copy
  // of local matrix entries and permuting it can only enforce
  // max aggregate size
  {
    if (params.get<bool>("aggregation: deterministic")) {
      Monitor m(*this, "BuildAggregatesDeterministic");
      BuildAggregatesDeterministic(maxNodesPerAggregate, graph,
                                   aggregates, aggStat, numNonAggregatedNodes);
    } else {
      Monitor m(*this, "BuildAggregatesRandom");
      BuildAggregatesRandom(maxNodesPerAggregate, graph,
                            aggregates, aggStat, numNonAggregatedNodes);
    }
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregatesRandom(const LO maxAggSize,
                          const LWGraph_kokkos& graph,
                          Aggregates& aggregates,
                          typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                          LO& numNonAggregatedNodes) const {
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  // Extract data from aggregates
  auto vertex2AggId = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner   = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors       = aggregates.GetGraphColors();

  auto lclLWGraph = graph;

  LO numAggregatedNodes = 0;
  LO numLocalAggregates = aggregates.GetNumAggregates();
  Kokkos::View<LO, device_type> aggCount("aggCount");
  Kokkos::deep_copy(aggCount, numLocalAggregates);
  Kokkos::parallel_for(
      "Aggregation Phase 1: initial reduction over color == 1",
      Kokkos::RangePolicy<LO, execution_space>(0, numRows),
      KOKKOS_LAMBDA(const LO nodeIdx) {
        if (colors(nodeIdx) == 1 && aggStat(nodeIdx) == READY) {
          const LO aggIdx          = Kokkos::atomic_fetch_add(&aggCount(), 1);
          vertex2AggId(nodeIdx, 0) = aggIdx;
          aggStat(nodeIdx)         = AGGREGATED;
          procWinner(nodeIdx, 0)   = myRank;
        }
      });
  // Truely we wish to compute: numAggregatedNodes = aggCount - numLocalAggregates
  // before updating the value of numLocalAggregates.
  // But since we also do not want to create a host mirror of aggCount we do some trickery...
  numAggregatedNodes -= numLocalAggregates;
  Kokkos::deep_copy(numLocalAggregates, aggCount);
  numAggregatedNodes += numLocalAggregates;

  // Compute the initial size of the aggregates.
  // Note lbv 12-21-17: I am pretty sure that the aggregates will always be of size 1
  //                    at this point so we could simplify the code below a lot if this
  //                    assumption is correct...
  Kokkos::View<LO*, device_type> aggSizesView("aggSizes", numLocalAggregates);
  {
    // Here there is a possibility that two vertices assigned to two different threads contribute
    // to the same aggregate if somethings happened before phase 1?
    auto aggSizesScatterView = Kokkos::Experimental::create_scatter_view(aggSizesView);
    Kokkos::parallel_for(
        "Aggregation Phase 1: compute initial aggregates size",
        Kokkos::RangePolicy<LO, execution_space>(0, numRows),
        KOKKOS_LAMBDA(const LO nodeIdx) {
          auto aggSizesScatterViewAccess = aggSizesScatterView.access();
          if (vertex2AggId(nodeIdx, 0) >= 0)
            aggSizesScatterViewAccess(vertex2AggId(nodeIdx, 0)) += 1;
        });
    Kokkos::Experimental::contribute(aggSizesView, aggSizesScatterView);
  }

  LO tmpNumAggregatedNodes = 0;
  Kokkos::parallel_reduce(
      "Aggregation Phase 1: main parallel_reduce over aggSizes",
      Kokkos::RangePolicy<size_t, execution_space>(0, numRows),
      KOKKOS_LAMBDA(const size_t nodeIdx, LO& lNumAggregatedNodes) {
        if (colors(nodeIdx) != 1 && (aggStat(nodeIdx) == READY || aggStat(nodeIdx) == NOTSEL)) {
          // Get neighbors of vertex i and look for local, aggregated,
          // color 1 neighbor (valid root).
          auto neighbors = lclLWGraph.getNeighborVertices(nodeIdx);
          for (LO j = 0; j < neighbors.length; ++j) {
            auto nei = neighbors.colidx(j);
            if (lclLWGraph.isLocalNeighborVertex(nei) && colors(nei) == 1 && aggStat(nei) == AGGREGATED) {
              // This atomic guarentees that any other node trying to
              // join aggregate agg has the correct size.
              LO agg           = vertex2AggId(nei, 0);
              const LO aggSize = Kokkos::atomic_fetch_add(&aggSizesView(agg),
                                                          1);
              if (aggSize < maxAggSize) {
                // assign vertex i to aggregate with root j
                vertex2AggId(nodeIdx, 0) = agg;
                procWinner(nodeIdx, 0)   = myRank;
                aggStat(nodeIdx)         = AGGREGATED;
                ++lNumAggregatedNodes;
                break;
              } else {
                // Decrement back the value of aggSizesView(agg)
                Kokkos::atomic_decrement(&aggSizesView(agg));
              }
            }
          }
        }
        // if(aggStat(nodeIdx) != AGGREGATED) {
        //   lNumNonAggregatedNodes++;
        if (aggStat(nodeIdx) == NOTSEL) {
          aggStat(nodeIdx) = READY;
        }
        // }
      },
      tmpNumAggregatedNodes);
  numAggregatedNodes += tmpNumAggregatedNodes;
  numNonAggregatedNodes -= numAggregatedNodes;

  // update aggregate object
  aggregates.SetNumAggregates(numLocalAggregates);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationPhase1Algorithm<LocalOrdinal, GlobalOrdinal, Node>::
    BuildAggregatesDeterministic(const LO maxAggSize,
                                 const LWGraph_kokkos& graph,
                                 Aggregates& aggregates,
                                 typename AggregationAlgorithmBase<LocalOrdinal, GlobalOrdinal, Node>::AggStatType& aggStat,
                                 LO& numNonAggregatedNodes) const {
  using device_type     = typename LWGraph_kokkos::device_type;
  using execution_space = typename LWGraph_kokkos::execution_space;

  const LO numRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  auto vertex2AggId = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto procWinner   = aggregates.GetProcWinner()->getDeviceLocalView(Xpetra::Access::ReadWrite);
  auto colors       = aggregates.GetGraphColors();

  auto lclLWGraph = graph;

  LO numLocalAggregates = aggregates.GetNumAggregates();
  Kokkos::View<LO, device_type> numLocalAggregatesView("Num aggregates");
  {
    auto h_nla = Kokkos::create_mirror_view(numLocalAggregatesView);
    h_nla()    = numLocalAggregates;
    Kokkos::deep_copy(numLocalAggregatesView, h_nla);
  }

  Kokkos::View<LO*, device_type> newRoots("New root LIDs", numNonAggregatedNodes);
  Kokkos::View<LO, device_type> numNewRoots("Number of new aggregates of current color");
  auto h_numNewRoots = Kokkos::create_mirror_view(numNewRoots);

  // first loop build the set of new roots
  Kokkos::parallel_for(
      "Aggregation Phase 1: building list of new roots",
      Kokkos::RangePolicy<execution_space>(0, numRows),
      KOKKOS_LAMBDA(const LO i) {
        if (colors(i) == 1 && aggStat(i) == READY) {
          // i will become a root
          newRoots(Kokkos::atomic_fetch_add(&numNewRoots(), 1)) = i;
        }
      });
  Kokkos::deep_copy(h_numNewRoots, numNewRoots);
  // sort new roots by LID to guarantee determinism in agg IDs
  Kokkos::sort(newRoots, 0, h_numNewRoots());
  LO numAggregated = 0;
  Kokkos::parallel_reduce(
      "Aggregation Phase 1: aggregating nodes",
      Kokkos::RangePolicy<execution_space>(0, h_numNewRoots()),
      KOKKOS_LAMBDA(const LO rootIndex, LO& lnumAggregated) {
        LO root               = newRoots(rootIndex);
        LO aggID              = numLocalAggregatesView() + rootIndex;
        LO aggSize            = 1;
        vertex2AggId(root, 0) = aggID;
        procWinner(root, 0)   = myRank;
        aggStat(root)         = AGGREGATED;
        auto neighOfRoot      = lclLWGraph.getNeighborVertices(root);
        for (LO n = 0; n < neighOfRoot.length; n++) {
          LO neigh = neighOfRoot(n);
          if (lclLWGraph.isLocalNeighborVertex(neigh) && aggStat(neigh) == READY) {
            // add neigh to aggregate
            vertex2AggId(neigh, 0) = aggID;
            procWinner(neigh, 0)   = myRank;
            aggStat(neigh)         = AGGREGATED;
            aggSize++;
            if (aggSize == maxAggSize) {
              // can't add any more nodes
              break;
            }
          }
        }
        lnumAggregated += aggSize;
      },
      numAggregated);
  numNonAggregatedNodes -= numAggregated;
  // update aggregate object
  aggregates.SetNumAggregates(numLocalAggregates + h_numNewRoots());
}

}  // namespace MueLu

#endif /* MUELU_AGGREGATIONPHASE1ALGORITHM_DEF_HPP_ */
