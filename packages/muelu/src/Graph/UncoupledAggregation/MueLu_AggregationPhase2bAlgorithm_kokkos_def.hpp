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
#ifndef MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_DEF_HPP
#define MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_DEF_HPP

#ifdef HAVE_MUELU_KOKKOS_REFACTOR

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase2bAlgorithm_kokkos_decl.hpp"

#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  // Try to stick unaggregated nodes into a neighboring aggregate if they are
  // not already too big
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationPhase2bAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregates(const ParameterList& params, const LWGraph_kokkos& graph,
                  Aggregates_kokkos& aggregates, Kokkos::View<unsigned*, typename MueLu::
                  LWGraph_kokkos<LO,GO,Node>::local_graph_type::device_type::
                  memory_space>& aggStatView, LO& numNonAggregatedNodes, Kokkos::View<LO*,
                  typename MueLu::LWGraph_kokkos<LO, GO, Node>::local_graph_type::device_type::
                  memory_space>& colorsDevice, LO& numColors) const {
    Monitor m(*this, "BuildAggregates");

    typedef typename MueLu::LWGraph_kokkos<LO, GO, Node>::local_graph_type graph_t;
    typedef typename graph_t::device_type::memory_space memory_space;
    typedef typename graph_t::device_type::execution_space execution_space;
    typedef typename graph_t::device_type::execution_space::scratch_memory_space scratch_space;
    typedef Kokkos::View<LO*,
                         scratch_space,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> > ScratchViewType;

    const LO  numRows = graph.GetNodeNumVertices();
    int myRank  = graph.GetComm()->getRank();

    bool colorPermuted          = params.get<bool>       ("aggregation: permute nodes by color");

    auto vertex2AggIdView = aggregates.GetVertex2AggId()->template getLocalView<memory_space>();
    auto procWinnerView   = aggregates.GetProcWinner()  ->template getLocalView<memory_space>();
    auto colorsCardinality = aggregates.GetColorsCardinality();
    auto colorPermutation = aggregates.GetColorPermutation();

    LO numLocalAggregates = aggregates.GetNumAggregates();

    const int defaultConnectWeight = 100;
    const int penaltyConnectWeight = 10;

    // This actually corresponds to the maximum number of entries per row in the matrix.
    const size_t maxNumNeighbors = graph.getNodeMaxNumRowEntries();
    int scratch_size = ScratchViewType::shmem_size( 3*maxNumNeighbors );

    Kokkos::View<int*, memory_space> connectWeightView("connectWeight", numRows);
    Kokkos::View<int*, memory_space> aggPenaltiesView ("aggPenalties",  numRows);

    Kokkos::parallel_for("Aggregation Phase 2b: Initialize connectWeightView", numRows,
                         KOKKOS_LAMBDA (const LO i) {
                           connectWeightView(i) = defaultConnectWeight;
                         });

    // taw: by running the aggregation routine more than once there is a chance that also
    // non-aggregated nodes with a node distance of two are added to existing aggregates.
    // Assuming that the aggregate size is 3 in each direction running the algorithm only twice
    // should be sufficient.
    // lbv: If the prior phase of aggregation where run without specifying an aggregate size,
    // the distance 2 coloring and phase 1 aggregation actually guarantee that only one iteration
    // is needed to reach distance 2 neighbors.
    int maxIters = 2;
    int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");
    if(maxNodesPerAggregate == std::numeric_limits<int>::max()) {maxIters = 1;}
    std::cout << "maxIters set to: " << maxIters << std::endl;
    for (int iter = 0; iter < maxIters; ++iter) {
      // total work = numberOfTeams * teamSize
      typedef typename Kokkos::TeamPolicy<execution_space>::member_type  member_type;
      LO tmpNumNonAggregatedNodes = 0;
      if(colorPermuted) {
        int colorRange = colorsCardinality(numColors) - colorsCardinality(1);
        int colorOffset = colorsCardinality(1);
        Kokkos::TeamPolicy<execution_space> outerPolicy(colorRange, Kokkos::AUTO);
        Kokkos::parallel_reduce("Aggregation Phase 2b: aggregates expansion",
                                outerPolicy.set_scratch_size( 0, Kokkos::PerTeam( scratch_size ) ),
                                KOKKOS_LAMBDA (const member_type &teamMember,
                                               LO& lNumNonAggregatedNodes) {

                                  // Retrieve the id of the vertex we are currently working on and
                                  // allocate view locally so that threads do not trash the weigth
                                  // when working on the same aggregate.
                                  const int vertexIdx = colorPermutation(teamMember.league_rank()
                                                                         + colorOffset);
                                  int numAggregatedNeighbors = 0;
                                  ScratchViewType aggregatedNeighbors(teamMember.team_scratch( 0 ),
                                                                      maxNumNeighbors);
                                  ScratchViewType vertex2AggLIDView(teamMember.team_scratch( 0 ),
                                                                    maxNumNeighbors);
                                  ScratchViewType aggWeightView(teamMember.team_scratch( 0 ),
                                                                maxNumNeighbors);

                                  if (aggStatView(vertexIdx) == READY) {

                                    // neighOfINode should become scratch and shared among
                                    // threads in the team...
                                    auto neighOfINode = graph.getNeighborVertices(vertexIdx);

                                    // create a mapping from neighbor "lid" to aggregate "lid"
                                    Kokkos::single( Kokkos::PerTeam( teamMember ), [&] () {
                                        int aggLIDCount = 0;
                                        for (int j = 0; j < neighOfINode.length; ++j) {
                                          LO neigh = neighOfINode(j);
                                          if( graph.isLocalNeighborVertex(neigh) &&
                                              (aggStatView(neigh) == AGGREGATED) ) {
                                            aggregatedNeighbors(numAggregatedNeighbors) = j;

                                            bool useNewLID = true;
                                            for(int k = 0; k < numAggregatedNeighbors; ++k) {
                                              LO lowerNeigh = neighOfINode(aggregatedNeighbors(k));
                                              if(vertex2AggIdView(neigh, 0)
                                                 == vertex2AggIdView(lowerNeigh, 0)) {
                                                vertex2AggLIDView(numAggregatedNeighbors) = 
                                                  vertex2AggLIDView(k);
                                                useNewLID = false;
                                              }
                                            }
                                            if(useNewLID) {
                                              vertex2AggLIDView(numAggregatedNeighbors) = aggLIDCount;
                                              ++aggLIDCount;
                                            }

                                            ++numAggregatedNeighbors;
                                          }
                                        }
                                      });

                                    for (int j = 0; j < numAggregatedNeighbors; j++) {
                                      LO localNeigh = aggregatedNeighbors(j);
                                      LO neigh = neighOfINode(j);

                                      aggWeightView(vertex2AggLIDView(j)) =
                                        aggWeightView(vertex2AggLIDView(j))
					+ connectWeightView(neigh);
                                    }

                                    int bestScore   = -100000;
                                    int bestAggId   = -1;
                                    int bestConnect = -1;

                                    for (int j = 0; j < numAggregatedNeighbors; j++) {
                                      LO localNeigh = aggregatedNeighbors(j);
                                      LO neigh = neighOfINode(localNeigh);
                                      int aggId = vertex2AggIdView(neigh, 0);
                                      int score = aggWeightView(vertex2AggLIDView(j))
                                        - aggPenaltiesView(aggId);

                                      if (score > bestScore) {
                                        bestAggId   = aggId;
                                        bestScore   = score;
                                        bestConnect = connectWeightView(neigh);

                                      } else if (aggId == bestAggId
                                                 && connectWeightView(neigh) > bestConnect) {
                                        bestConnect = connectWeightView(neigh);
                                      }

                                      // Reset the weights for the next loop
                                      // LBV: this looks a little suspicious, it would probably
                                      // need to be taken out of this inner for loop...
                                      aggWeightView(vertex2AggLIDView(j)) = 0;
                                    }

                                    // Do the actual aggregate update with a single thread!
                                    Kokkos::single( Kokkos::PerTeam( teamMember ), [&] () {
                                        if (bestScore >= 0) {
                                          aggStatView     (vertexIdx)    = AGGREGATED;
                                          vertex2AggIdView(vertexIdx, 0) = bestAggId;
                                          procWinnerView  (vertexIdx, 0) = myRank;

                                          lNumNonAggregatedNodes--;

                                          // This does not protect bestAggId's aggPenalties from being
                                          // fetched by another thread before this update happens, it just
                                          // guarantees that the update is performed correctly...
                                          Kokkos::atomic_add(&aggPenaltiesView(bestAggId), 1);
                                          connectWeightView(vertexIdx) = bestConnect
                                            - penaltyConnectWeight;
                                        }
                                      });
                                  }
                                }, tmpNumNonAggregatedNodes);
      } else {
        Kokkos::TeamPolicy<execution_space> outerPolicy(numRows, Kokkos::AUTO);
        Kokkos::parallel_reduce("Aggregation Phase 2b: aggregates expansion",
                                outerPolicy.set_scratch_size( 0, Kokkos::PerTeam( scratch_size ) ),
                                KOKKOS_LAMBDA (const member_type &teamMember,
                                               LO& lNumNonAggregatedNodes) {

                                  // Retrieve the id of the vertex we are currently working on and
                                  // allocate view locally so that threads do not trash the weigth
                                  // when working on the same aggregate.
                                  const int vertexIdx = teamMember.league_rank();
                                  int numAggregatedNeighbors = 0;
                                  ScratchViewType aggregatedNeighbors(teamMember.team_scratch( 0 ),
                                                                      maxNumNeighbors);
                                  ScratchViewType vertex2AggLIDView(teamMember.team_scratch( 0 ),
                                                                    maxNumNeighbors);
                                  ScratchViewType aggWeightView(teamMember.team_scratch( 0 ),
                                                                maxNumNeighbors);

                                  if (aggStatView(vertexIdx) == READY) {

                                    // neighOfINode should become scratch and shared among
                                    // threads in the team...
                                    auto neighOfINode = graph.getNeighborVertices(vertexIdx);

                                    // create a mapping from neighbor "lid" to aggregate "lid"
                                    Kokkos::single( Kokkos::PerTeam( teamMember ), [&] () {
                                        int aggLIDCount = 0;
                                        for (int j = 0; j < neighOfINode.length; ++j) {
                                          LO neigh = neighOfINode(j);
                                          if( graph.isLocalNeighborVertex(neigh) &&
                                              (aggStatView(neigh) == AGGREGATED) ) {
                                            aggregatedNeighbors(numAggregatedNeighbors) = j;

                                            bool useNewLID = true;
                                            for(int k = 0; k < numAggregatedNeighbors; ++k) {
                                              LO lowerNeigh = neighOfINode(aggregatedNeighbors(k));
                                              if(vertex2AggIdView(neigh, 0)
                                                 == vertex2AggIdView(lowerNeigh, 0)) {
                                                vertex2AggLIDView(numAggregatedNeighbors) = 
                                                  vertex2AggLIDView(k);
                                                useNewLID = false;
                                              }
                                            }
                                            if(useNewLID) {
                                              vertex2AggLIDView(numAggregatedNeighbors) = aggLIDCount;
                                              ++aggLIDCount;
                                            }

                                            ++numAggregatedNeighbors;
                                          }
                                        }
                                      });

                                    for (int j = 0; j < numAggregatedNeighbors; j++) {
                                      LO localNeigh = aggregatedNeighbors(j);
                                      LO neigh = neighOfINode(j);

                                      aggWeightView(vertex2AggLIDView(j)) =
                                        aggWeightView(vertex2AggLIDView(j))
					+ connectWeightView(neigh);
                                    }

                                    int bestScore   = -100000;
                                    int bestAggId   = -1;
                                    int bestConnect = -1;

                                    for (int j = 0; j < numAggregatedNeighbors; j++) {
                                      LO localNeigh = aggregatedNeighbors(j);
                                      LO neigh = neighOfINode(localNeigh);
                                      int aggId = vertex2AggIdView(neigh, 0);
                                      int score = aggWeightView(vertex2AggLIDView(j))
                                        - aggPenaltiesView(aggId);

                                      if (score > bestScore) {
                                        bestAggId   = aggId;
                                        bestScore   = score;
                                        bestConnect = connectWeightView(neigh);

                                      } else if (aggId == bestAggId
                                                 && connectWeightView(neigh) > bestConnect) {
                                        bestConnect = connectWeightView(neigh);
                                      }

                                      // Reset the weights for the next loop
                                      // LBV: this looks a little suspicious, it would probably
                                      // need to be taken out of this inner for loop...
                                      aggWeightView(vertex2AggLIDView(j)) = 0;
                                    }

                                    // Do the actual aggregate update with a single thread!
                                    Kokkos::single( Kokkos::PerTeam( teamMember ), [&] () {
                                        if (bestScore >= 0) {
                                          aggStatView     (vertexIdx)    = AGGREGATED;
                                          vertex2AggIdView(vertexIdx, 0) = bestAggId;
                                          procWinnerView  (vertexIdx, 0) = myRank;

                                          lNumNonAggregatedNodes--;

                                          // This does not protect bestAggId's aggPenalties from being
                                          // fetched by another thread before this update happens, it just
                                          // guarantees that the update is performed correctly...
                                          Kokkos::atomic_add(&aggPenaltiesView(bestAggId), 1);
                                          connectWeightView(vertexIdx) = bestConnect
                                            - penaltyConnectWeight;
                                        }
                                      });
                                  }
                                }, tmpNumNonAggregatedNodes);
      }
      numNonAggregatedNodes += tmpNumNonAggregatedNodes;
    } // loop over k

  }

} // end namespace

#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_DEF_HPP
