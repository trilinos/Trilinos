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

#ifdef HAVE_MUELU_KOKKOS_REFACTOR

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationPhase2aAlgorithm_kokkos_decl.hpp"

#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationPhase2aAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregates(const ParameterList& params, const LWGraph_kokkos& graph,
                  Aggregates_kokkos& aggregates, Kokkos::View<unsigned*, typename MueLu::
                  LWGraph_kokkos<LO,GO,Node>::local_graph_type::device_type::
                  memory_space>& aggStatView, LO& numNonAggregatedNodes,
                  Kokkos::View<LO*, typename MueLu::LWGraph_kokkos<LO, GO, Node>::
                  local_graph_type::device_type::memory_space>& colorsDevice, LO& numColors) const {
    Monitor m(*this, "BuildAggregates");

    typedef typename MueLu::LWGraph_kokkos<LO, GO, Node>::local_graph_type graph_t;
    typedef typename graph_t::device_type::memory_space memory_space;

    int minNodesPerAggregate = params.get<int>("aggregation: min agg size");
    int maxNodesPerAggregate = params.get<int>("aggregation: max agg size");

    const LO  numRows = graph.GetNodeNumVertices();
    const int myRank  = graph.GetComm()->getRank();

    auto vertex2AggId = aggregates.GetVertex2AggId()->template getLocalView<memory_space>();
    auto procWinner   = aggregates.GetProcWinner()  ->template getLocalView<memory_space>();

    LO numLocalNodes      = numRows;
    LO numLocalAggregated = numLocalNodes - numNonAggregatedNodes;

    const double aggFactor = 0.5;
    double       factor    = as<double>(numLocalAggregated)/(numLocalNodes+1);
    factor = pow(factor, aggFactor);

    Kokkos::View<LO, memory_space> numLocalAggregates("numLocalAggregates");
    typename Kokkos::View<LO, memory_space>::HostMirror h_numLocalAggregates =
      Kokkos::create_mirror_view(numLocalAggregates);
    h_numLocalAggregates() = aggregates.GetNumAggregates();
    Kokkos::deep_copy(numLocalAggregates, h_numLocalAggregates);

    // Now we create new aggregates using root nodes in all colors other than the first color,
    // as the first color was already exhausted in Phase 1.
    for(int color = 2; color < numColors + 1; ++color) {

      LO tmpNumNonAggregatedNodes = 0;
      Kokkos::parallel_reduce("Aggregation Phase 2a: loop over each individual color", numRows,
                              KOKKOS_LAMBDA (const LO rootCandidate, LO& lNumNonAggregatedNodes) {
                                if(aggStatView(rootCandidate) == READY &&
                                   colorsDevice(rootCandidate) == color) {

                                  LO aggSize = 0;
                                  auto neighbors = graph.getNeighborVertices(rootCandidate);

                                  // Loop over neighbors to count how many nodes could join
                                  // the new aggregate
                                  LO numNeighbors = 0;
                                  for(int j = 0; j < neighbors.length; ++j) {
                                    LO neigh = neighbors(j);
                                    if(neigh != rootCandidate) {
                                      if(graph.isLocalNeighborVertex(neigh) &&
                                         aggStatView(neigh) == READY &&
                                         aggSize < maxNodesPerAggregate) {
                                        // aggList(aggSize) = neigh;
                                        ++aggSize;
                                      }
                                      ++numNeighbors;
                                    }
                                  }

                                  // If a sufficient number of nodes can join the new aggregate
                                  // then we actually create the aggregate.
                                  if(aggSize > minNodesPerAggregate &&
                                     aggSize > factor*numNeighbors) {

                                    // aggregates.SetIsRoot(rootCandidate);
                                    LO aggIndex = Kokkos::
                                      atomic_fetch_add(&numLocalAggregates(), 1);

                                    for(int j = 0; j < neighbors.length; ++j) {
                                      LO neigh = neighbors(j);
                                      if(neigh != rootCandidate) {
                                        if(graph.isLocalNeighborVertex(neigh) &&
                                           aggStatView(neigh) == READY &&
                                           aggSize < maxNodesPerAggregate) {
                                          aggStatView(neigh)   = AGGREGATED;
                                          vertex2AggId(neigh, 0) = aggIndex;
                                          procWinner(neigh, 0)   = myRank;
                                        }
                                      }
                                    }
                                    lNumNonAggregatedNodes -= aggSize;
                                  }
                                }
                              }, tmpNumNonAggregatedNodes);
      numNonAggregatedNodes += tmpNumNonAggregatedNodes;
    }

    // update aggregate object
    Kokkos::deep_copy(h_numLocalAggregates, numLocalAggregates);
    aggregates.SetNumAggregates(h_numLocalAggregates());
  }

} // end namespace

#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_AGGREGATIONPHASE2AALGORITHM_KOKKOS_DEF_HPP
