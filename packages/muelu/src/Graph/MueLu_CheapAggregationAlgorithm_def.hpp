/*
 * MueLu_CheapAggregationAlgorithm_def.hpp
 *
 *  Created on: Jul 25, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CHEAPAGGREGATIONALGORITHM_DEF_HPP_
#define MUELU_CHEAPAGGREGATIONALGORITHM_DEF_HPP_

#include <queue>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_CheapAggregationAlgorithm_decl.hpp"

#include "MueLu_Graph.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_LinkedList.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CheapAggregationAlgorithm(RCP<FactoryBase> const &graphFact)
    : ordering_(NATURAL), minNodesPerAggregate_(1), maxNeighAlreadySelected_(0)
  { }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  LocalOrdinal CheapAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoarsenUncoupled(Graph const & graph, Aggregates & aggregates) const {
    Monitor m(*this, "Coarsen Uncoupled (CheapAggregationAlgorithm)");

    RCP<const Map> gBoundaryNodeMap = graph.GetBoundaryNodeMap();

    std::string orderingType;
    switch (ordering_) {
    case NATURAL:
      orderingType = "Natural";
      break;
/*    case RANDOM:
      orderingType = "Random";
      break;*/
    case GRAPH:
      orderingType = "Graph";
      break;
    default:
      break;
    }


    // aggregate state
    Teuchos::ArrayRCP<NodeState> aggStat;
    const LocalOrdinal nRows = graph.GetNodeNumVertices();
    if(nRows > 0) aggStat = Teuchos::arcp<NodeState>(nRows);
    for(LocalOrdinal i=0; i<nRows; ++i) {
      aggStat[i] = READY;
      if(gBoundaryNodeMap->isNodeGlobalElement(graph.GetDomainMap()->getGlobalElement(i))) aggStat[i] = BDRY;
    }


    // vertex ids for output
    Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);

    // some internal variables
    LocalOrdinal nAggregates = 0;   // number of local aggregates on current proc
    std::queue<LocalOrdinal> graph_ordering_inodes; // inodes for graph ordering
    LocalOrdinal iNode2 = 0;        // local iteration variable
    LocalOrdinal iNode  = 0;        // current node

    // main loop over all local rows of grpah(A)
    while (iNode2 < nRows) {

      // pick the next node to aggregate
      if      (ordering_ == NATURAL) iNode = iNode2++;
      /*else if (ordering_ == RANDOM )*/  // NOT SUPPORTED YET
      else if (ordering_ == GRAPH) {
        // if there are no nodes for graph ordering scheme
        if(graph_ordering_inodes.size() == 0) {
          // add exactly one ready node for graph ordering aggregates
          for(LocalOrdinal jnode=0; jnode<nRows; jnode++) {
            if(aggStat[jnode] == READY) {
              graph_ordering_inodes.push(jnode);
              break;
            }
          }
        }
        if(graph_ordering_inodes.size()==0) break; // there's no ready node any more -> end phase 1
        iNode = graph_ordering_inodes.front(); // take next node from graph ordering queue
        graph_ordering_inodes.pop();           // delete this node in list
      }

      // consider iNode only if it is in READY state
      if(aggStat[iNode] == READY) {
        // build new aggregate
        Aggregate ag;
        ag.list.push_back(iNode);

        // extract column information from graph for current row on current proc
        Teuchos::ArrayView<const LocalOrdinal> neighOfINode = graph.getNeighborVertices(iNode);
        //typename Teuchos::ArrayView<const LocalOrdinal>::size_type nnz = neighOfINode.size();

        LocalOrdinal cnt_neighbours = 0;  // number of possible neighbour nodes for current new aggregate

        // build tentative new aggregate
        for( typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
            LocalOrdinal index = *it;

            // note: this is uncoupled coarsening
            // only column indices of current proc are allowed
            if(graph.isLocalNeighborVertex(index)) {
              // check status of current neighbor node
              if(aggStat[index]==READY ||
                 aggStat[index]==NOTSEL ) {
                ag.list.push_back(index); // add neighbor node to current aggregate
              }
              else if(aggStat[index]!=BDRY){
                cnt_neighbours++;
              }
            } // end if: current column index belongs to this proc
        } // end for: loop over all columns in current row

        // if there are too many neighbours aggregated or the number of nodes
        // in the new aggregate is too few, don't do this one.

        // check if aggregate ag is acceptable
        if((cnt_neighbours > GetMaxNeighAlreadySelected() ) || // aggregate is too close to bdry nodes
          (ag.list.size() < (unsigned int) GetMinNodesPerAggregate())) {      // not enough nodes in new aggregate
          // failed to build a new aggregate
          ag.list.clear();
          aggStat[iNode] = NOTSEL; // node not valid to be supernode, mark it as NOTSEL
          if(ordering_ == GRAPH) {
            // even though the aggregate around iNode is not perfect, we try the ndoes where iNode is connected to
            // loop over all column indices
            for (typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              LocalOrdinal index = *it;
              if(graph.isLocalNeighborVertex(index) && aggStat[index] == READY) // if node connected to iNode is not aggregated
                graph_ordering_inodes.push(index);
            }
          }
        } else {
          // accept new aggregate
          aggregates.SetIsRoot(iNode);    // mark iNode as root node for new aggregate 'ag'
          ag.index = nAggregates++;       // aggregate accepted, increment aggregate counter
          vertex2AggId[iNode] = ag.index;
          for (unsigned int k=0; k<ag.list.size(); k++) {
            aggStat[ag.list[k]] = SELECTED;
            vertex2AggId[ag.list[k]] = ag.index;
            if(ordering_ == GRAPH) {
              Teuchos::ArrayView<const LocalOrdinal> neighOfJNode = graph.getNeighborVertices(ag.list[k]);
              for(typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator it = neighOfJNode.begin(); it!=neighOfJNode.end(); ++it) {
                LocalOrdinal index = *it;
                if(graph.isLocalNeighborVertex(index) && aggStat[index] == READY)
                  graph_ordering_inodes.push(index);
              }
            } // end GRAPH
          } // end loop over all nodes in aggregate
        } // end if accept aggs or decline aggs


      } // end if aggStat[iNode] == READY

    } // end while

    // update aggregate object
    aggregates.SetNumAggregates(nAggregates);

    // clean up
    if(graph_ordering_inodes.size()>0){
      for(unsigned int k=0; k<graph_ordering_inodes.size(); k++)
        graph_ordering_inodes.pop();
    }

    // verbose
    {
      const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();

      if(IsPrint(Warnings0)) {
        GO localReady  = 0;
        GO globalReady = 0;

        // compute number of local nodes with status "ready"
        for(LO i=0; i<nRows; ++i)
          if(aggStat[i] == READY) localReady++;
        // compute global number of nodes with status "ready"
        sumAll(comm, localReady, globalReady);
        if(globalReady > 0)
          GetOStream(Warnings0,0) << "Aggregation (UC): Phase 1 (WARNING) " << globalReady << " unaggregated nodes left (status READY)" << std::endl;
      }

      if(IsPrint(Statistics1)) {
        LO localSelected  = 0;
        GO globalSelected = 0;
        GO globalNRows    = 0;
        for(LO i=0; i<nRows; ++i)
          if(aggStat[i] == SELECTED) localSelected++;
        sumAll(comm, localSelected, globalSelected);
        sumAll(comm, nRows, globalNRows);
        GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 1: Nodes aggregated = " << globalSelected << " out of " << globalNRows << " nodes" << std::endl;
        GO nAggregatesGlobal = 0;
        sumAll(comm, nAggregates, nAggregatesGlobal);
        GetOStream(Statistics1, 0) << "Aggregation (UC): Phase 1: Total aggregates = " << nAggregatesGlobal << std::endl;
      }
    }

    // collect some local information
    LO nLocalSelected    = 0;
    LO nLocalBdry        = 0;
    LO nLocalNotSelected = 0;
    LO nLocalReady       = 0;
    for (LO iNode = 0; iNode < nRows; iNode++) {
      if      (aggStat[iNode] == SELECTED) nLocalSelected++;
      else if (aggStat[iNode] == BDRY)     nLocalBdry++;
      else if (aggStat[iNode] == NOTSEL)   nLocalNotSelected++;
      else if (aggStat[iNode] == READY)    nLocalReady++;
    }

    return nLocalReady + nLocalNotSelected;
  }
} // end namespace


#endif /* MUELU_CHEAPAGGREGATIONALGORITHM_DEF_HPP_ */
