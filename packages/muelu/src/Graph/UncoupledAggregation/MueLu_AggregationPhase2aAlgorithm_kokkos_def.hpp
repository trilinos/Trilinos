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
#ifndef MUELU_AGGREGATIONPHASE2AALGORITHM_KOKKOS_DEF_HPP
#define MUELU_AGGREGATIONPHASE2AALGORITHM_KOKKOS_DEF_HPP

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase2aAlgorithm_kokkos_decl.hpp"

#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Monitor.hpp"

#include "Kokkos_Sort.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationPhase2aAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregates(const ParameterList& params,
                  const LWGraph_kokkos& graph,
                  Aggregates_kokkos& aggregates,
                  Kokkos::View<unsigned*, typename LWGraph_kokkos::device_type>& aggStat,
                  LO& numNonAggregatedNodes) const {

    if(params.get<bool>("aggregation: deterministic")) {
      Monitor m(*this, "BuildAggregatesDeterministic");
      BuildAggregatesDeterministic(params, graph, aggregates, aggStat, numNonAggregatedNodes);
    } else {
      Monitor m(*this, "BuildAggregatesRandom");
      BuildAggregatesRandom(params, graph, aggregates, aggStat, numNonAggregatedNodes);
    }

  } // BuildAggregates

  template <class LO, class GO, class Node>
  void AggregationPhase2aAlgorithm_kokkos<LO, GO, Node>::
  BuildAggregatesRandom(const ParameterList& params,
                        const LWGraph_kokkos& graph,
                        Aggregates_kokkos& aggregates,
                        Kokkos::View<unsigned*, typename LWGraph_kokkos::device_type>& aggStat,
                        LO& numNonAggregatedNodes) const
  {
    const int minNodesPerAggregate = params.get<int>("aggregation: min agg size");
    const int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");
    bool matchMLbehavior = params.get<bool>("aggregation: match ML phase2a");

    const LO  numRows = graph.GetNodeNumVertices();
    const int myRank  = graph.GetComm()->getRank();

    auto vertex2AggId  = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
    auto procWinner    = aggregates.GetProcWinner()  ->getDeviceLocalView(Xpetra::Access::ReadWrite);
    auto colors        = aggregates.GetGraphColors();
    const LO numColors = aggregates.GetGraphNumColors();

    auto lclLWGraph = graph.getLocalLWGraph();

    LO numLocalNodes      = numRows;
    LO numLocalAggregated = numLocalNodes - numNonAggregatedNodes;

    const double aggFactor = 0.5;
    double       factor    = static_cast<double>(numLocalAggregated)/(numLocalNodes+1);
    factor = pow(factor, aggFactor);

    // LBV on Sept 12, 2019: this looks a little heavy handed,
    // I'm not sure a view is needed to perform atomic updates.
    // If we can avoid this and use a simple LO that would be
    // simpler for later maintenance.
    Kokkos::View<LO, device_type> numLocalAggregates("numLocalAggregates");
    typename Kokkos::View<LO, device_type>::HostMirror h_numLocalAggregates =
      Kokkos::create_mirror_view(numLocalAggregates);
    h_numLocalAggregates() = aggregates.GetNumAggregates();
    Kokkos::deep_copy(numLocalAggregates, h_numLocalAggregates);

    // Now we create new aggregates using root nodes in all colors other than the first color,
    // as the first color was already exhausted in Phase 1.
    for(int color = 2; color < numColors + 1; ++color) {
      LO tmpNumNonAggregatedNodes = 0;
      Kokkos::parallel_reduce("Aggregation Phase 2a: loop over each individual color",
                              Kokkos::RangePolicy<execution_space>(0, numRows),
                              KOKKOS_LAMBDA (const LO rootCandidate, LO& lNumNonAggregatedNodes) {
                                if(aggStat(rootCandidate) == READY &&
                                   colors(rootCandidate) == color) {

                                  LO numNeighbors = 0;
                                  LO aggSize = 0;
                                  if (matchMLbehavior) {
                                    aggSize += 1;
                                    numNeighbors +=1;
                                  }

                                  auto neighbors = lclLWGraph.getNeighborVertices(rootCandidate);

                                  // Loop over neighbors to count how many nodes could join
                                  // the new aggregate

                                  for(int j = 0; j < neighbors.length; ++j) {
                                    LO neigh = neighbors(j);
                                    if(neigh != rootCandidate) {
                                      if(lclLWGraph.isLocalNeighborVertex(neigh) &&
                                         (aggStat(neigh) == READY) &&
                                         (aggSize < maxNodesPerAggregate)) {
                                        ++aggSize;
                                      }
                                      ++numNeighbors;
                                    }
                                  }

                                  // If a sufficient number of nodes can join the new aggregate
                                  // then we actually create the aggregate.
                                  if(aggSize > minNodesPerAggregate &&
                                     (aggSize > factor*numNeighbors)) {

                                    // aggregates.SetIsRoot(rootCandidate);
                                    LO aggIndex = Kokkos::
                                      atomic_fetch_add(&numLocalAggregates(), 1);

                                    LO numAggregated = 0;

                                    if (matchMLbehavior) {
                                      // Add the root.
                                      aggStat(rootCandidate)         = AGGREGATED;
                                      vertex2AggId(rootCandidate, 0) = aggIndex;
                                      procWinner(rootCandidate, 0)   = myRank;
                                      ++numAggregated;
                                      --lNumNonAggregatedNodes;
                                    }

                                    for(int neighIdx = 0; neighIdx < neighbors.length; ++neighIdx) {
                                      LO neigh = neighbors(neighIdx);
                                      if(neigh != rootCandidate) {
                                        if(lclLWGraph.isLocalNeighborVertex(neigh) &&
                                           (aggStat(neigh) == READY) &&
                                           (numAggregated < aggSize)) {
                                          aggStat(neigh)         = AGGREGATED;
                                          vertex2AggId(neigh, 0) = aggIndex;
                                          procWinner(neigh, 0)   = myRank;

                                          ++numAggregated;
                                          --lNumNonAggregatedNodes;
                                        }
                                      }
                                    }
                                  }
                                }
                              }, tmpNumNonAggregatedNodes);
      numNonAggregatedNodes += tmpNumNonAggregatedNodes;
    }

    // update aggregate object
    Kokkos::deep_copy(h_numLocalAggregates, numLocalAggregates);
    aggregates.SetNumAggregates(h_numLocalAggregates());
  } // BuildAggregatesRandom

  template <class LO, class GO, class Node>
  void AggregationPhase2aAlgorithm_kokkos<LO, GO, Node>::
  BuildAggregatesDeterministic(const ParameterList& params,
                               const LWGraph_kokkos& graph,
                               Aggregates_kokkos& aggregates,
                               Kokkos::View<unsigned*, typename LWGraph_kokkos::device_type>& aggStat,
                               LO& numNonAggregatedNodes) const
  {
    const int minNodesPerAggregate = params.get<int>("aggregation: min agg size");
    const int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");

    const LO  numRows = graph.GetNodeNumVertices();
    const int myRank  = graph.GetComm()->getRank();

    auto vertex2AggId  = aggregates.GetVertex2AggId()->getDeviceLocalView(Xpetra::Access::ReadWrite);
    auto procWinner    = aggregates.GetProcWinner()  ->getDeviceLocalView(Xpetra::Access::ReadWrite);
    auto colors        = aggregates.GetGraphColors();
    const LO numColors = aggregates.GetGraphNumColors();

    auto lclLWGraph = graph.getLocalLWGraph();

    LO numLocalNodes      = procWinner.size();
    LO numLocalAggregated = numLocalNodes - numNonAggregatedNodes;

    const double aggFactor = 0.5;
    double       factor    = as<double>(numLocalAggregated)/(numLocalNodes+1);
    factor = pow(factor, aggFactor);

    Kokkos::View<LO, device_type> numLocalAggregates("numLocalAggregates");
    typename Kokkos::View<LO, device_type>::HostMirror h_numLocalAggregates =
      Kokkos::create_mirror_view(numLocalAggregates);
    h_numLocalAggregates() = aggregates.GetNumAggregates();
    Kokkos::deep_copy(numLocalAggregates, h_numLocalAggregates);

    // Now we create new aggregates using root nodes in all colors other than the first color,
    // as the first color was already exhausted in Phase 1.
    //
    // In the deterministic version, exactly the same set of aggregates will be created
    // (as the nondeterministic version)
    // because no vertex V can be a neighbor of two vertices of the same color, so two root
    // candidates can't fight over V
    //
    // But, the precise values in vertex2AggId need to match exactly, so just sort the new
    // roots of each color before assigning aggregate IDs

    //numNonAggregatedNodes is the best available upper bound for the number of aggregates
    //which may be created in this phase, so use it for the size of newRoots
    Kokkos::View<LO*, device_type> newRoots("New root LIDs", numNonAggregatedNodes);
    Kokkos::View<LO, device_type> numNewRoots("Number of new aggregates of current color");
    auto h_numNewRoots = Kokkos::create_mirror_view(numNewRoots);
    for(int color = 1; color < numColors + 1; ++color) {
      h_numNewRoots() = 0;
      Kokkos::deep_copy(numNewRoots, h_numNewRoots);
      Kokkos::parallel_for("Aggregation Phase 2a: determining new roots of current color",
                           Kokkos::RangePolicy<execution_space>(0, numRows),
                           KOKKOS_LAMBDA(const LO rootCandidate) {
                             if(aggStat(rootCandidate) == READY &&
                                colors(rootCandidate) == color) {
                               LO aggSize = 0;
                               auto neighbors = lclLWGraph.getNeighborVertices(rootCandidate);
                               // Loop over neighbors to count how many nodes could join
                               // the new aggregate
                               LO numNeighbors = 0;
                               for(int j = 0; j < neighbors.length; ++j) {
                                 LO neigh = neighbors(j);
                                 if(neigh != rootCandidate)
                                   {
                                     if(lclLWGraph.isLocalNeighborVertex(neigh) &&
                                        aggStat(neigh) == READY &&
                                        aggSize < maxNodesPerAggregate)
                                       {
                                         ++aggSize;
                                       }
                                     ++numNeighbors;
                                   }
                               }
                               // If a sufficient number of nodes can join the new aggregate
                               // then we mark rootCandidate as a future root.
                               if(aggSize > minNodesPerAggregate && aggSize > factor*numNeighbors) {
                                 LO newRootIndex = Kokkos::atomic_fetch_add(&numNewRoots(), 1);
                                 newRoots(newRootIndex) = rootCandidate;
                               }
                             }
                           });
      Kokkos::deep_copy(h_numNewRoots, numNewRoots);

      if(h_numNewRoots() > 0) {
        //sort the new root indices
        Kokkos::sort(newRoots, 0, h_numNewRoots());
        //now, loop over all new roots again and actually create the aggregates
        LO tmpNumNonAggregatedNodes = 0;
        //First, just find the set of color vertices which will become aggregate roots
        Kokkos::parallel_reduce("Aggregation Phase 2a: create new aggregates",
                                Kokkos::RangePolicy<execution_space>(0, h_numNewRoots()),
                                KOKKOS_LAMBDA (const LO newRootIndex, LO& lNumNonAggregatedNodes) {
                                  LO root = newRoots(newRootIndex);
                                  LO newAggID = numLocalAggregates() + newRootIndex;
                                  auto neighbors = lclLWGraph.getNeighborVertices(root);
                                  // Loop over neighbors and add them to new aggregate
                                  aggStat(root)      = AGGREGATED;
                                  vertex2AggId(root, 0) = newAggID;
                                  LO aggSize = 1;
                                  for(int j = 0; j < neighbors.length; ++j) {
                                    LO neigh = neighbors(j);
                                    if(neigh != root) {
                                      if(lclLWGraph.isLocalNeighborVertex(neigh) &&
                                         aggStat(neigh) == READY &&
                                         aggSize < maxNodesPerAggregate) {
                                        aggStat(neigh)      = AGGREGATED;
                                        vertex2AggId(neigh, 0) = newAggID;
                                        procWinner(neigh, 0)   = myRank;
                                        aggSize++;
                                      }
                                    }
                                  }
                                  lNumNonAggregatedNodes -= aggSize;
                                }, tmpNumNonAggregatedNodes);
        numNonAggregatedNodes += tmpNumNonAggregatedNodes;
        h_numLocalAggregates() += h_numNewRoots();
        Kokkos::deep_copy(numLocalAggregates, h_numLocalAggregates);
      }
    }
    aggregates.SetNumAggregates(h_numLocalAggregates());
  }

} // end namespace

#endif // MUELU_AGGREGATIONPHASE2AALGORITHM_KOKKOS_DEF_HPP
