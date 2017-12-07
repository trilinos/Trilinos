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
    const int myRank  = graph.GetComm()->getRank();

    auto vertex2AggIdView = aggregates.GetVertex2AggId()->template getLocalView<memory_space>();
    auto procWinnerView   = aggregates.GetProcWinner()  ->template getLocalView<memory_space>();

    LO numLocalAggregates = aggregates.GetNumAggregates();

    const int defaultConnectWeight = 100;
    const int penaltyConnectWeight = 10;

    const size_t maxNumNeighbors = graph.getNodeMaxNumRowEntries();
    int scratch_size = ScratchViewType::shmem_size( maxNumNeighbors );

    Kokkos::View<int*, memory_space> connectWeightView("connectWeight", numRows);
    Kokkos::View<int*, memory_space> aggPenaltiesView ("aggPenalties",  numRows);

    typename Kokkos::View<int*, memory_space>::HostMirror h_connectWeightView =
      Kokkos::create_mirror_view (connectWeightView);
    Kokkos::parallel_for("Aggregation Phase 2b: Initialize connectWeightView", numRows,
                         KOKKOS_LAMBDA (const LO i) {
                           h_connectWeightView(i) = defaultConnectWeight;
                         });
    Kokkos::deep_copy(connectWeightView,    h_connectWeightView);

    // We do this cycle twice.
    // I don't know why, but ML does it too
    // taw: by running the aggregation routine more than once there is a chance that also
    // non-aggregated nodes with a node distance of two are added to existing aggregates.
    // Assuming that the aggregate size is 3 in each direction running the algorithm only twice
    // should be sufficient.
    for (int k = 0; k < 2; k++) {
      // total work = numberOfTeams * teamSize
      Kokkos::TeamPolicy<execution_space> outerPolicy(numRows, Kokkos::AUTO);
      typedef typename Kokkos::TeamPolicy<execution_space>::member_type  member_type;
      // Kokkos::RangePolicy<execution_space> numRowsPolicy(0, numRows);
      LO tmpNumNonAggregatedNodes = 0;
      Kokkos::parallel_reduce("Aggregation Phase 2b: aggregates expansion",
                              // numRowsPolicy,
                              outerPolicy.set_scratch_size( 0, Kokkos::PerTeam( scratch_size ) ),
                              KOKKOS_LAMBDA (const member_type &teamMember,
                                             LO& lNumNonAggregatedNodes) {

                                // Retrieve the id of the vertex we are currently working on and
                                // allocate view locally so that threads do not trash the weigth
                                // when working on the same aggregate.
                                const int vertexIdx = teamMember.league_rank();
                                ScratchViewType aggWeightView(teamMember.team_scratch( 0 ),
                                                              maxNumNeighbors);

                                if (aggStatView(vertexIdx) == READY) {

                                  auto neighOfINode = graph.getNeighborVertices(vertexIdx);

                                  for (int j = 0; j < as<int>(neighOfINode.length); j++) {
                                    LO neigh = neighOfINode(j);

                                    // We don't check (neigh != i), as it is covered by checking
                                    // (aggStat[neigh] == AGGREGATED)
                                    if ( graph.isLocalNeighborVertex(neigh)
                                         && (aggStatView(neigh) == AGGREGATED) )
                                      aggWeightView(j) = aggWeightView(j)
                                        + connectWeightView(neigh);
                                  }

                                  int bestScore   = -100000;
                                  int bestAggId   = -1;
                                  int bestConnect = -1;

                                  for (int j = 0; j < as<int>(neighOfINode.length); j++) {
                                    LO neigh = neighOfINode(j);

                                    if ( graph.isLocalNeighborVertex(neigh)
                                         && (aggStatView(neigh) == AGGREGATED) ) {
                                      int aggId = vertex2AggIdView(neigh, 0);
                                      int score = aggWeightView(j) - aggPenaltiesView(aggId);

                                      if (score > bestScore) {
                                        bestAggId   = aggId;
                                        bestScore   = score;
                                        bestConnect = connectWeightView(neigh);

                                      } else if (aggId == bestAggId
                                                 && connectWeightView(neigh) > bestConnect) {
                                        bestConnect = connectWeightView(neigh);
                                      }

                                      // Reset the weights for the next loop
                                      aggWeightView(j) = 0;
                                    }
                                  }

                                  if (bestScore >= 0) {
                                    aggStatView   (vertexIdx)    = AGGREGATED;
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
                                }
                              }, tmpNumNonAggregatedNodes);
      numNonAggregatedNodes += tmpNumNonAggregatedNodes;
    } // loop over k

  }

} // end namespace

#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_AGGREGATIONPHASE2BALGORITHM_KOKKOS_DEF_HPP
