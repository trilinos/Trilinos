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
#ifndef MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DEF_HPP_
#define MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DEF_HPP_


#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_AggregationStructuredAlgorithm_decl.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_IndexManager.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationStructuredAlgorithm<LocalOrdinal, GlobalOrdinal, Node>::
  BuildAggregates(const ParameterList& params, const GraphBase& graph, Aggregates& aggregates, std::vector<unsigned>& aggStat,
                  LO& numNonAggregatedNodes) const {
    Monitor m(*this, "BuildAggregates");

    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_STRUCTUREDALGORITHM_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    std::string coupling = params.get<std::string>("aggregation: coupling");
    const bool coupled = (coupling == "coupled" ? true : false);
    ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
    ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()  ->getDataNonConst(0);
    RCP<IndexManager> geoData = aggregates.GetIndexManager();
    Array<LO>  ghostedCoarseNodeCoarseLIDs;
    Array<int> ghostedCoarseNodeCoarsePIDs;

    *out << "Extract data for ghosted nodes" << std::endl;
    geoData->getGhostedNodesData(graph.GetDomainMap(), ghostedCoarseNodeCoarseLIDs,
                                 ghostedCoarseNodeCoarsePIDs);

    LO iGhosted, jGhosted, kGhosted, iCoarse, jCoarse, kCoarse, iRem, jRem, kRem;
    LO ghostedCoarseNodeCoarseLID, aggId, rate;
    *out << "Loop over fine nodes and assign them to an aggregate and a rank" << std::endl;
    for(LO nodeIdx = 0; nodeIdx < geoData->getNumLocalFineNodes(); ++nodeIdx) {
      // Compute coarse ID associated with fine LID
      geoData->getFineNodeGhostedTuple(nodeIdx, iGhosted, jGhosted, kGhosted);
      iCoarse = iGhosted / geoData->getCoarseningRate(0);
      iRem    = iGhosted % geoData->getCoarseningRate(0);
      if(iGhosted - geoData->getOffset(0)
         < geoData->getLocalFineNodesInDir(0) - geoData->getCoarseningEndRate(0)) {
        rate = geoData->getCoarseningRate(0);
      } else {
        rate = geoData->getCoarseningEndRate(0);
      }
      if(iRem > (rate / 2)) { ++iCoarse; }
      if(coupled && (geoData->getStartGhostedCoarseNode(0)*geoData->getCoarseningRate(0)
                     > geoData->getStartIndex(0))) { --iCoarse; }
      jCoarse = jGhosted / geoData->getCoarseningRate(1);
      jRem    = jGhosted % geoData->getCoarseningRate(1);
      if(jGhosted - geoData->getOffset(1)
         < geoData->getLocalFineNodesInDir(1) - geoData->getCoarseningEndRate(1)) {
        rate = geoData->getCoarseningRate(1);
      } else {
        rate = geoData->getCoarseningEndRate(1);
      }
      if(jRem > (rate / 2)) { ++jCoarse; }
      if(coupled && (geoData->getStartGhostedCoarseNode(1)*geoData->getCoarseningRate(1)
                     > geoData->getStartIndex(1))) { --jCoarse; }
      kCoarse = kGhosted / geoData->getCoarseningRate(2);
      kRem    = kGhosted % geoData->getCoarseningRate(2);
      if(kGhosted - geoData->getOffset(2)
         < geoData->getLocalFineNodesInDir(2) - geoData->getCoarseningEndRate(2)) {
        rate = geoData->getCoarseningRate(2);
      } else {
        rate = geoData->getCoarseningEndRate(2);
      }
      if(kRem > (rate / 2)) { ++kCoarse; }
      if(coupled && (geoData->getStartGhostedCoarseNode(2)*geoData->getCoarseningRate(2)
                     > geoData->getStartIndex(2))) { --kCoarse; }
      geoData->getCoarseNodeGhostedLID(iCoarse, jCoarse, kCoarse, ghostedCoarseNodeCoarseLID);

      aggId                 = ghostedCoarseNodeCoarseLIDs[ghostedCoarseNodeCoarseLID];
      vertex2AggId[nodeIdx] = aggId;
      procWinner[nodeIdx]   = ghostedCoarseNodeCoarsePIDs[ghostedCoarseNodeCoarseLID];
      aggStat[nodeIdx]      = AGGREGATED;
      --numNonAggregatedNodes;
    } // Loop over fine points

  }


} // end namespace


#endif /* MUELU_AGGREGATIONSTRUCTUREDALGORITHM_DEF_HPP_ */
