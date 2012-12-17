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
 * MueLu_UncoupledAggregationAlgorithm_def.hpp
 *
 *  Created on: Sep 17, 2012
 *      Author: Tobias Wiesner
 */

#ifndef MUELU_UNCOUPLEDAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_UNCOUPLEDAGGREGATIONALGORITHM_DEF_HPP_

#include <queue>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_UncoupledAggregationAlgorithm.hpp"

#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
UncoupledAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UncoupledAggregationAlgorithm(RCP<const FactoryBase> const &graphFact)
{
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
LocalOrdinal UncoupledAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildAggregates(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<unsigned int> & aggStat) const {
  Monitor m(*this, "Coarsen Uncoupled (UncoupledAggregationAlgorithm)");

  std::string orderingType;
  switch (this->GetOrdering()) {
  case NATURAL:
    orderingType = "Natural";
    break;
  case RANDOM:
    orderingType = "Random";
    break;
  case GRAPH:
    orderingType = "Graph";
    break;
  default:
    break;
  }

  const LocalOrdinal nRows = graph.GetNodeNumVertices();
  const int myRank = graph.GetComm()->getRank();

  // vertex ids for output
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);

  // some internal variables
  LocalOrdinal nLocalAggregates = aggregates.GetNumAggregates();    // number of local aggregates on current proc
  std::queue<LocalOrdinal> graph_ordering_inodes; // inodes for graph ordering
  LocalOrdinal iNode2 = 0;        // local iteration variable
  LocalOrdinal iNode1  = 0;        // current node
  Teuchos::ArrayRCP<LO> randomVector;

  if ( this->GetOrdering() == RANDOM ) {
    randomVector = Teuchos::arcp<LO>(nRows); //size_t or int ?-> to be propagated
    for (LocalOrdinal i = 0; i < nRows; ++i) randomVector[i] = i;
    RandomReorder(randomVector);
  }

  // main loop over all local rows of grpah(A)
  while (iNode2 < nRows) {

    // pick the next node to aggregate
    if      (this->GetOrdering() == NATURAL) iNode1 = iNode2++;
    else if (this->GetOrdering()== RANDOM ) iNode1 = randomVector[iNode2++];
    else if (this->GetOrdering() == GRAPH) {
      // if there are no nodes for graph ordering scheme
      if(graph_ordering_inodes.size() == 0) {
        // add exactly one ready node for graph ordering aggregates
        for(LocalOrdinal jnode=0; jnode<nRows; jnode++) {
          if(aggStat[jnode] == NodeStats::READY) {
            graph_ordering_inodes.push(jnode);
            break;
          }
        }
      }
      if(graph_ordering_inodes.size()==0) break; // there's no ready node any more -> end phase 1
      iNode1 = graph_ordering_inodes.front(); // take next node from graph ordering queue
      graph_ordering_inodes.pop();           // delete this node in list
    }

    // consider iNode1 only if it is not aggregated yet
    if(aggStat[iNode1] == NodeStats::READY ) {
      // build new aggregate
      Aggregate ag;
      ag.list.push_back(iNode1);

      // extract column information from graph for current row on current proc
      Teuchos::ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(iNode1);

      LocalOrdinal cnt_neighbours = 0;  // number of possible neighbour nodes for current new aggregate

      // build tentative new aggregate
      for( typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
        // note: this is uncoupled coarsening
        // only column indices of current proc are allowed
        if(graph.isLocalNeighborVertex(*it)) {
          // check status of current neighbor node
          if(aggStat[*it] == NodeStats::READY ||
              aggStat[*it] == NodeStats::NOTSEL) {
            ag.list.push_back(*it); // add neighbor node to current aggregate
          }
          else
            cnt_neighbours++; // increment number of neighbour nodes that are already aggregated
        } // end if: current column index belongs to this proc
      } // end for: loop over all columns in current row

      // if there are too many neighbours aggregated or the number of nodes
      // in the new aggregate is too few, don't do this one.

      // check if aggregate ag is acceptable
      if((cnt_neighbours > this->GetMaxNeighAlreadySelected() ) || // aggregate has too many nodes that are already aggregated
          (ag.list.size() < (unsigned int) this->GetMinNodesPerAggregate())) {      // not enough nodes in new aggregate
        // failed to build a new aggregate
        ag.list.clear();
        aggStat[iNode1] = NodeStats::NOTSEL;
        if(this->GetOrdering() == GRAPH) {
          // even though the aggregate around iNode1 is not perfect, we try the ndoes where iNode1 is connected to
          // loop over all column indices
          for (typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
            if(graph.isLocalNeighborVertex(*it) &&
                aggStat[*it] == NodeStats::READY )
              graph_ordering_inodes.push(*it);
          }
        }
      } else {
        // accept new aggregate
        aggregates.SetIsRoot(iNode1);    // mark iNode1 as root node for new aggregate 'ag'
        ag.index = nLocalAggregates++;       // aggregate accepted, increment aggregate counter
        vertex2AggId[iNode1] = ag.index;
        procWinner[iNode1] = myRank; //graph.GetComm()->getRank();
        /*std::cout << "build new aggregate of size " << ag.list.size() << " nodes" << std::endl;
        std::cout << "nodes: ";
        for (unsigned int k=0; k<ag.list.size(); k++)
          std::cout << ag.list[k] << " ";
        std::cout << std::endl;*/

        for (unsigned int k=0; k<ag.list.size(); k++) {
          aggStat[ag.list[k]] = NodeStats::AGGREGATED;  // mark node as aggregated
          vertex2AggId[ag.list[k]] = ag.index;  // fill vertex2AggId and procWinner structure with information
          procWinner[ag.list[k]] = myRank;
          if(this->GetOrdering() == GRAPH) {
            Teuchos::ArrayView<const LocalOrdinal> neighOfJNode = graph.getNeighborVertices(ag.list[k]);
            for(typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfJNode.begin(); it!=neighOfJNode.end(); ++it) {
              if(graph.isLocalNeighborVertex(*it) &&  // TODO check me (index < nRows)
                  aggStat[*it] == NodeStats::READY) // jnode not aggregated and not selected
                graph_ordering_inodes.push(*it);
            }
          } // end GRAPH
        } // end loop over all nodes in aggregate
      } // end if accept aggs or decline aggs
    } // end if "node not aggregated yet
  } // end while

  // update aggregate object
  aggregates.SetNumAggregates(nLocalAggregates);

  // clean up
  if(graph_ordering_inodes.size()>0){
    for(unsigned int k=0; k<graph_ordering_inodes.size(); k++)
      graph_ordering_inodes.pop();
  }

  // print aggregation information
  this->PrintAggregationInformation("UncoupledAggregationAlgorithm:", graph, aggregates, aggStat);

  // collect some local information
  LO nLocalAggregated    = 0;
  LO nLocalNotAggregated = 0;
  for (LO i = 0; i < nRows; i++) {
    if (aggStat[i] == NodeStats::AGGREGATED) nLocalAggregated++;
    else nLocalNotAggregated++;
  }

  return nLocalNotAggregated;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void UncoupledAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomReorder(Teuchos::ArrayRCP<LO> list) const {
  //TODO: replace int
  int n = list.size();
  for(int i=0; i<n-1; i++) {
    std::swap(list[i], list[RandomOrdinal(i,n-1)]);
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
int UncoupledAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomOrdinal(int min, int max) const {
  return min + static_cast<int>((max-min+1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
}

} // end namespace


#endif /* MUELU_UNCOUPLEDAGGREGATIONALGORITHM_DEF_HPP_ */
