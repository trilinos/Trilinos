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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_MaxLinkAggregationAlgorithm_def.hpp
 *
 *  Created on: Sep 18, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_MAXLINKAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_MAXLINKAGGREGATIONALGORITHM_DEF_HPP_

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_MaxLinkAggregationAlgorithm.hpp"

#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
//#include "MueLu_LinkedList.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
MaxLinkAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MaxLinkAggregationAlgorithm(RCP<const FactoryBase> const &graphFact)
{
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
LocalOrdinal MaxLinkAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildAggregates(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat) const {
  Monitor m(*this, "BuildAggregates");

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  const LocalOrdinal nRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  // loop over all local rows
  for (LocalOrdinal iNode=0; iNode<Teuchos::as<LocalOrdinal>(nRows); iNode++) {
    if(aggStat[iNode] != NodeStats::AGGREGATED ) {
      LocalOrdinal selected_aggregate = -1;

      std::map<LocalOrdinal,LocalOrdinal> aggid2cntconnections;

      Teuchos::ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(iNode);
      for(typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it!=neighOfINode.end(); ++it) {
        LocalOrdinal index = *it;
        if(graph.isLocalNeighborVertex(index) && aggStat[index] == NodeStats::AGGREGATED) {
          LocalOrdinal aggid = vertex2AggId[index]; // get (local) aggregate id
          if(aggid2cntconnections.count(aggid)>0) aggid2cntconnections[aggid] += 1;
          else aggid2cntconnections[aggid] = 1;
        } // if node is aggregated but not a one-pt node

      } // loop over all columns

      // find aggregate id with most connections
      LocalOrdinal maxcnt = 0;
      for(typename std::map<LocalOrdinal,LocalOrdinal>::const_iterator it=aggid2cntconnections.begin(); it!=aggid2cntconnections.end();++it) {
        if(maxcnt < it->second) {
          maxcnt = it->second;
          selected_aggregate = it->first;
        }
      }

      // add node iNode to aggregate
      if(selected_aggregate != -1) {
        aggStat[iNode] = NodeStats::AGGREGATED;
        vertex2AggId[iNode] = selected_aggregate;
        procWinner[iNode] = myRank; //graph.GetComm()->getRank();
      }
    } // end if aggState == NOTSEL...
  } // end loop over all local rows

  // print aggregation information
  this->PrintAggregationInformation("Phase 2 (max_link, extend aggregates):", graph, aggregates, aggStat);

  // collect some local information
  LO nLocalAggregated    = 0;
  LO nLocalNotAggregated = 0;
  for (LO i = 0; i < nRows; i++) {
    if (aggStat[i] == NodeStats::AGGREGATED) nLocalAggregated++;
    else nLocalNotAggregated++;
  }

  return nLocalNotAggregated;
}

} // end namespace


#endif /* MUELU_MAXLINKAGGREGATIONALGORITHM_DEF_HPP_ */
