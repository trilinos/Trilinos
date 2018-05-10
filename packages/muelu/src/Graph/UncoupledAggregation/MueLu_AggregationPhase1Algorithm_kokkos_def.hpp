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

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationPhase1Algorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregates(const ParameterList& params, const LWGraph_kokkos& graph, Aggregates_kokkos& aggregates, std::vector<unsigned>& aggStat,
                  LO& numNonAggregatedNodes) const {
    Monitor m(*this, "BuildAggregates");

    std::string orderingStr     = params.get<std::string>("aggregation: ordering");
    int maxNeighAlreadySelected = params.get<int>        ("aggregation: max selected neighbors");
    int minNodesPerAggregate    = params.get<int>        ("aggregation: min agg size");
    int maxNodesPerAggregate    = params.get<int>        ("aggregation: max agg size");

    Algorithm algorithm         = Algorithm::Serial;
    std::string algoParamName = "aggregation: phase 1 algorithm";
    if(params.isParameter(algoParamName))
    {
      algorithm = algorithmFromName(params.get<std::string>("aggregation: phase 1 algorithm"));
    }

    TEUCHOS_TEST_FOR_EXCEPTION(maxNodesPerAggregate < minNodesPerAggregate, Exceptions::RuntimeError,
                               "MueLu::UncoupledAggregationAlgorithm::BuildAggregates: minNodesPerAggregate must be smaller or equal to MaxNodePerAggregate!");

    //Distance-2 gives less control than serial uncoupled phase 1
    //no custom row reordering because would require making deep copy of local matrix entries and permuting it
    //can only enforce max aggregate size
    if(algorithm == Algorithm::Distance2)
    {
      BuildAggregatesDistance2(graph, aggregates, aggStat, numNonAggregatedNodes, maxNodesPerAggregate);
    }
    else
    {
      BuildAggregatesSerial(graph, aggregates, aggStat, numNonAggregatedNodes,
          minNodesPerAggregate, maxNodesPerAggregate, maxNeighAlreadySelected, orderingStr);
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
      std::vector<unsigned>& aggStat, LO& numNonAggregatedNodes, LO maxAggSize) const
  {
    const LO  numRows = graph.GetNodeNumVertices();
    const int myRank  = graph.GetComm()->getRank();

    ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);

    LO numLocalAggregates = aggregates.GetNumAggregates();

    //get the sparse local graph in CRS
    std::vector<LocalOrdinal> rowptrs;
    rowptrs.reserve(numRows + 1);
    std::vector<LocalOrdinal> colinds;
    colinds.reserve(graph.GetNodeNumEdges());

    rowptrs.push_back(0);
    for(LocalOrdinal row = 0; row < numRows; row++)
    {
      auto entries = graph.getNeighborVertices(row);
      for(LocalOrdinal i = 0; i < entries.length; i++)
      {
        colinds.push_back(entries.colidx(i));
      }
      rowptrs.push_back(colinds.size());
    }

    //the local CRS graph to Kokkos device views, then compute graph squared
    typedef typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type graph_t;
    typedef typename graph_t::device_type device_t;
    typedef typename device_t::memory_space memory_space;
    typedef typename device_t::execution_space execution_space;
    typedef typename graph_t::row_map_type::non_const_type rowptrs_view;
    typedef Kokkos::View<size_t*, Kokkos::HostSpace> host_rowptrs_view;
    typedef typename graph_t::entries_type::non_const_type colinds_view;
    typedef Kokkos::View<LocalOrdinal*, Kokkos::HostSpace> host_colinds_view;
    //note: just using colinds_view in place of scalar_view_t type (it won't be used at all by symbolic SPGEMM)
    typedef KokkosKernels::Experimental::KokkosKernelsHandle<
      typename rowptrs_view::const_value_type, typename colinds_view::const_value_type, typename colinds_view::const_value_type, 
      execution_space, memory_space, memory_space> KernelHandle;

    KernelHandle kh;
    //leave gc algorithm choice as the default
    kh.create_graph_coloring_handle();

    //Create device views for graph rowptrs/colinds
    rowptrs_view aRowptrs("A device rowptrs", numRows + 1);
    colinds_view aColinds("A device colinds", colinds.size());
    // Populate A in temporary host views, then copy to device
    {
      host_rowptrs_view aHostRowptrs("A host rowptrs", numRows + 1);
      for(LO i = 0; i < numRows + 1; i++)
      {
        aHostRowptrs(i) = rowptrs[i];
      }
      Kokkos::deep_copy(aRowptrs, aHostRowptrs);
      host_colinds_view aHostColinds("A host colinds", colinds.size());
      for(size_t i = 0; i < colinds.size(); i++)
      {
        aHostColinds(i) = colinds[i];
      }
      Kokkos::deep_copy(aColinds, aHostColinds);
    }
    //run d2 graph coloring
    //graph is symmetric so row map/entries and col map/entries are the same
    KokkosGraph::Experimental::d2_graph_color(&kh, numRows, numRows, aRowptrs, aColinds, aRowptrs, aColinds);

    // extract the colors
    auto coloringHandle = kh.get_graph_coloring_handle();
    auto colorsDevice = coloringHandle->get_vertex_colors();

    auto colors = Kokkos::create_mirror_view(colorsDevice);
    Kokkos::deep_copy(colors, colorsDevice);

    //clean up coloring handle
    kh.destroy_graph_coloring_handle();

    //have color 1 (first color) be the aggregate roots (add those to mapping first)
    LocalOrdinal aggCount = 0;
    for(LocalOrdinal i = 0; i < numRows; i++)
    {
      if(colors(i) == 1 && aggStat[i] == READY)
      {
        vertex2AggId[i] = aggCount++;
        aggStat[i] = AGGREGATED;
        numLocalAggregates++;
        procWinner[i] = myRank;
      }
    }
    numNonAggregatedNodes = 0;
    std::vector<LocalOrdinal> aggSizes(numLocalAggregates, 0);
    for(int i = 0; i < numRows; i++)
    {
      if(vertex2AggId[i] >= 0)
        aggSizes[vertex2AggId[i]]++;
    }
    //now assign every READY vertex to a directly connected root
    for(LocalOrdinal i = 0; i < numRows; i++)
    {
      if(colors(i) != 1 && (aggStat[i] == READY || aggStat[i] == NOTSEL))
      {
        //get neighbors of vertex i and
        //look for local, aggregated, color 1 neighbor (valid root)
        auto neighbors = graph.getNeighborVertices(i);
        for(LocalOrdinal j = 0; j < neighbors.length; j++)
        {
          auto nei = neighbors.colidx(j);
          LocalOrdinal agg = vertex2AggId[nei];
          if(graph.isLocalNeighborVertex(nei) && colors(nei) == 1 && aggStat[nei] == AGGREGATED && aggSizes[agg] < maxAggSize)
          {
            //assign vertex i to aggregate with root j
            vertex2AggId[i] = agg;
            aggSizes[agg]++;
            aggStat[i] = AGGREGATED;
            procWinner[i] = myRank;
            break;
          }
        }
      }
      if(aggStat[i] != AGGREGATED)
      {
        numNonAggregatedNodes++;
        if(aggStat[i] == NOTSEL)
          aggStat[i] = READY;
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

