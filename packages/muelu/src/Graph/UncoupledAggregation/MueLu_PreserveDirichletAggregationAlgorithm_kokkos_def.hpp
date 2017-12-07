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
#ifndef MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_KOKKOS_DEF_HPP
#define MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_KOKKOS_DEF_HPP

#ifdef HAVE_MUELU_KOKKOS_REFACTOR

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_PreserveDirichletAggregationAlgorithm_kokkos.hpp"

#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_Aggregates_kokkos.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void PreserveDirichletAggregationAlgorithm_kokkos<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregates(Teuchos::ParameterList const & params, LWGraph_kokkos const & graph,
                  Aggregates_kokkos & aggregates, Kokkos::View<unsigned*, typename MueLu::
                  LWGraph_kokkos<LO,GO,Node>::local_graph_type::device_type::
                  memory_space>& aggStatView, LO& numNonAggregatedNodes, Kokkos::View<LO*,
                  typename MueLu::LWGraph_kokkos<LO, GO, Node>::local_graph_type::device_type::
                  memory_space>& colorsDevice, LO& numColors) const {
    Monitor m(*this, "BuildAggregates");

    bool preserve = params.get<bool>("aggregation: preserve Dirichlet points");

    const LO  numRows = graph.GetNodeNumVertices();
    const int myRank  = graph.GetComm()->getRank();

    typedef typename MueLu::LWGraph_kokkos<LO, GO, Node>::local_graph_type graph_t;
    typedef typename graph_t::device_type::memory_space memory_space;

    LO numLocalAggregates = aggregates.GetNumAggregates();

    auto vertex2AggIdView = aggregates.GetVertex2AggId()->template getLocalView<memory_space>();
    auto procWinnerView = aggregates.GetProcWinner()    ->template getLocalView<memory_space>();

    LO numAggregatedNodes = 0;
    if(preserve) {
      // Note LBV 03/20/2018
      // This does not look very good to me as the level of contention due to the atomic operation
      // can potentially be high. It should not be to bad as long as the ratio
      // (BC nodes)/(total nodes) is not too large...
      Kokkos::View<LO, memory_space> countBC("countBC");
      Kokkos::parallel_for("Aggregation preserve boundary nodes (preserve)", numRows,
                           KOKKOS_LAMBDA(const LO i) {
                             if(aggStatView(i) == BOUNDARY) {
                               // update the count atomically across threads for correctness
                               const LO idx = Kokkos::atomic_fetch_add (&countBC(), 1);
                               aggStatView(i) = IGNORED;
                               vertex2AggIdView(i, 0) = idx + numLocalAggregates;
                               procWinnerView(i, 0)   = myRank;
                             // aggregates.SetIsRoot(i);
                             }
                           });

      // Extract the number of aggregated nodes, which is also the number of new aggregates.
      Kokkos::deep_copy(numAggregatedNodes, countBC);
      numLocalAggregates += numAggregatedNodes;
    } else {
      // This could be handled in a single loop as it is done in serial but I like the idea
      // of removing the if statement from GPU kernels and the reduction below could be running
      // faster than the above code that includes atomics...
      Kokkos::parallel_reduce("Aggregation preserve boundary nodes (no-preserve)", numRows,
                              KOKKOS_LAMBDA(const LO i, LO & aggregatedNodesAccumulator) {
                                if(aggStatView(i) == BOUNDARY) {
                                  aggStatView(i) = IGNORED;
                                  aggregatedNodesAccumulator++;
                                }
                              }, numAggregatedNodes);
    }
    // update the number of non aggregated nodes.
    numNonAggregatedNodes -= numAggregatedNodes;

    // update aggregate object
    aggregates.SetNumAggregates(numLocalAggregates);
  }

} // end namespace

#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_KOKKOS_DEF_HPP
