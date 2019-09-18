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
  BuildAggregates(Teuchos::ParameterList const & params,
                  LWGraph_kokkos const & graph,
                  Aggregates_kokkos & aggregates,
                  Kokkos::View<unsigned*, typename LWGraph_kokkos::memory_space>& aggstat,
                  LO& numNonAggregatedNodes) const {
    Monitor m(*this, "BuildAggregates");

    typename Kokkos::View<unsigned*, memory_space>::HostMirror aggstatHost
      = Kokkos::create_mirror(aggstat);
    Kokkos::deep_copy(aggstatHost, aggstat);
    std::vector<unsigned> aggStat;
    aggStat.resize(aggstatHost.extent(0));
    for(size_t idx = 0; idx < aggstatHost.extent(0); ++idx) {
      aggStat[idx] = aggstatHost(idx);
    }

    bool preserve = params.get<bool>("aggregation: preserve Dirichlet points");

    const LO  numRows = graph.GetNodeNumVertices();
    const int myRank  = graph.GetComm()->getRank();

    ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);

    LO numLocalAggregates = aggregates.GetNumAggregates();

    for (LO i = 0; i < numRows; i++)
      if (aggStat[i] == BOUNDARY) {
        aggStat[i] = IGNORED;
        numNonAggregatedNodes--;

        if (preserve) {
          aggregates.SetIsRoot(i);

          vertex2AggId[i] = numLocalAggregates++;
          procWinner  [i] = myRank;
        }
      }

    for(size_t idx = 0; idx < aggstatHost.extent(0); ++idx) {
      aggstatHost(idx) = aggStat[idx];
    }
    Kokkos::deep_copy(aggstat, aggstatHost);

    // update aggregate object
    aggregates.SetNumAggregates(numLocalAggregates);
  }

} // end namespace

#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_PRESERVEDIRICHLETAGGREGATIONALGORITHM_KOKKOS_DEF_HPP
