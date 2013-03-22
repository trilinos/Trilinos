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
#ifndef MUELU_LOCALAGGREGATIONALGORITHM_DEF_HPP
#define MUELU_LOCALAGGREGATIONALGORITHM_DEF_HPP

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Xpetra_Vector.hpp>

#include "MueLu_LocalAggregationAlgorithm_decl.hpp"

#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_LinkedList.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::LocalAggregationAlgorithm()
    : ordering_(NATURAL), minNodesPerAggregate_(1), maxNeighAlreadySelected_(0)
  { }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoarsenUncoupled(GraphBase const & graph, Aggregates & aggregates) const {
    Monitor m(*this, "Coarsen Uncoupled");

    std::string orderingType;
    switch(ordering_) {
      case NATURAL:
        orderingType="Natural";
        break;
      case RANDOM:
        orderingType="Random";
        break;
      case GRAPH:
        orderingType="Graph";
        break;
      default:
        break;
    }
    GetOStream(Runtime1, 0) << "Ordering:                  " << orderingType << std::endl;
    GetOStream(Runtime1, 0) << "Min nodes per aggregate:   " << minNodesPerAggregate_ << std::endl;
    GetOStream(Runtime1, 0) << "Max nbrs already selected: " << maxNeighAlreadySelected_ << std::endl;

    /* Create Aggregation object */
    my_size_t nAggregates = 0;

    /* ============================================================= */
    /* aggStat indicates whether this node has been aggreated, and   */
    /* vertex2AggId stores the aggregate number where this node has  */
    /* been aggregated into.                                         */
    /* ============================================================= */

    Teuchos::ArrayRCP<NodeState> aggStat;
    const my_size_t nRows = graph.GetNodeNumVertices();
    if (nRows > 0) aggStat = Teuchos::arcp<NodeState>(nRows);
    for ( my_size_t i = 0; i < nRows; ++i ) aggStat[i] = READY;

    /* ============================================================= */
    /* Phase 1  :                                                    */
    /*    for all nodes, form a new aggregate with its neighbors     */
    /*    if the number of its neighbors having been aggregated does */
    /*    not exceed a given threshold                               */
    /*    (GetMaxNeighAlreadySelected() = 0 ===> Vanek's scheme) */
    /* ============================================================= */

    /* some general variable declarations */
    Teuchos::ArrayRCP<LO> randomVector;
    RCP<MueLu::LinkedList> nodeList; /* list storing the next node to pick as a root point for ordering_ == GRAPH */
    MueLu_SuperNode  *aggHead=NULL, *aggCurrent=NULL, *supernode=NULL;
    /**/

    if ( ordering_ == RANDOM )       /* random ordering */
      {
        //TODO: could be stored in a class that respect interface of LinkedList

        randomVector = Teuchos::arcp<LO>(nRows); //size_t or int ?-> to be propagated
        for (my_size_t i = 0; i < nRows; ++i) randomVector[i] = i;
        RandomReorder(randomVector);
      }
    else if ( ordering_ == GRAPH )  /* graph ordering */
      {
        nodeList = rcp(new MueLu::LinkedList());
        nodeList->Add(0);
      }

    /* main loop */
    {
      LO iNode  = 0;
      LO iNode2 = 0;

      Teuchos::ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0); // output only: contents ignored

      while (iNode2 < nRows)
        {

          /*------------------------------------------------------ */
          /* pick the next node to aggregate                       */
          /*------------------------------------------------------ */

          if      ( ordering_ == NATURAL ) iNode = iNode2++;
          else if ( ordering_ == RANDOM )  iNode = randomVector[iNode2++];
          else if ( ordering_ == GRAPH )
            {
              if ( nodeList->IsEmpty() )
                {
                  for ( int jNode = 0; jNode < nRows; ++jNode )
                    {
                      if ( aggStat[jNode] == READY )
                        {
                          nodeList->Add(jNode); //TODO optim: not necessary to create a node. Can just set iNode value and skip the end
                          break;
                        }
                    }
                }
              if ( nodeList->IsEmpty() ) break; /* end of the while loop */ //TODO: coding style :(

              iNode = nodeList->Pop();
            }
          else {
            throw(Exceptions::RuntimeError("CoarsenUncoupled: bad aggregation ordering option"));
          }

          /*------------------------------------------------------ */
          /* consider further only if the node is in READY mode    */
          /*------------------------------------------------------ */

          if ( aggStat[iNode] == READY )
            {
              // neighOfINode is the neighbor node list of node 'iNode'.
              Teuchos::ArrayView<const LO> neighOfINode = graph.getNeighborVertices(iNode);
              typename Teuchos::ArrayView<const LO>::size_type length = neighOfINode.size();

              supernode = new MueLu_SuperNode;
              try {
                supernode->list = Teuchos::arcp<int>(length+1);
              } catch (std::bad_alloc&) {
                TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::LocalAggregationAlgorithm::CoarsenUncoupled(): Error: couldn't allocate memory for supernode! length=" + Teuchos::toString(length));
              }

              supernode->maxLength = length;
              supernode->length = 1;
              supernode->list[0] = iNode;

              int selectFlag = 1;
              {
                /*--------------------------------------------------- */
                /* count the no. of neighbors having been aggregated  */
                /*--------------------------------------------------- */

                int count = 0;
                for (typename Teuchos::ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
                  {
                    int index = *it;
                    if ( index < nRows )
                      {
                        if ( aggStat[index] == READY ||
                             aggStat[index] == NOTSEL )
                          supernode->list[supernode->length++] = index;
                        else count++;

                      }
                  }

                /*--------------------------------------------------- */
                /* if there are too many neighbors aggregated or the  */
                /* number of nodes in the new aggregate is too few,   */
                /* don't do this one                                  */
                /*--------------------------------------------------- */

                if ( count > GetMaxNeighAlreadySelected() ) selectFlag = 0;
              }

              // Note: the supernode length is actually 1 more than the
              //       number of nodes in the candidate aggregate. The
              //       root is counted twice. I'm not sure if this is
              //       a bug or a feature ... so I'll leave it and change
              //       < to <= in the if just below.

              if (selectFlag != 1 ||
                  supernode->length <= GetMinNodesPerAggregate())
                {
                  aggStat[iNode] = NOTSEL;
                  delete supernode;
                  if ( ordering_ == GRAPH ) /* if graph ordering */
                    {
                      for (typename Teuchos::ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
                        {
                          int index = *it;
                          if  ( index < nRows && aggStat[index] == READY )
                            {
                              nodeList->Add(index);
                            }
                        }
                    }
                }
              else
                {
                  aggregates.SetIsRoot(iNode);
                  for ( int j = 0; j < supernode->length; ++j )
                    {
                      int jNode = supernode->list[j];
                      aggStat[jNode] = SELECTED;
                      vertex2AggId[jNode] = nAggregates;
                      if ( ordering_ == GRAPH ) /* if graph ordering */
                        {

                          Teuchos::ArrayView<const LO> neighOfJNode = graph.getNeighborVertices(jNode);

                          for (typename Teuchos::ArrayView<const LO>::const_iterator it = neighOfJNode.begin(); it != neighOfJNode.end(); ++it)
                            {
                              int index = *it;
                              if ( index < nRows && aggStat[index] == READY )
                                {
                                  nodeList->Add(index);
                                }
                            }
                        }
                    }
                  supernode->next = NULL;
                  supernode->index = nAggregates;
                  if ( nAggregates == 0 )
                    {
                      aggHead = supernode;
                      aggCurrent = supernode;
                    }
                  else
                    {
                      aggCurrent->next = supernode;
                      aggCurrent = supernode;
                    }
                  nAggregates++;
                  // unused aggCntArray[nAggregates] = supernode->length;
                }
            }
        } // end of 'for'

      // views on distributed vectors are freed here.

    } // end of 'main loop'

    nodeList = Teuchos::null;

    /* Update aggregate object */
    aggregates.SetNumAggregates(nAggregates);

    /* Verbose */
    {
      const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();

      if (IsPrint(Warnings0)) {
        GO localReady=0, globalReady;

        // Compute 'localReady'
        for ( my_size_t i = 0; i < nRows; ++i )
          if (aggStat[i] == READY) localReady++;

        // Compute 'globalReady'
        sumAll(comm, localReady, globalReady);

        if(globalReady > 0)
          GetOStream(Warnings0, 0) << "Warning: " << globalReady << " READY nodes left" << std::endl;
      }

      if (IsPrint(Statistics1)) {
        // Compute 'localSelected'
        LO localSelected=0;
        for ( my_size_t i = 0; i < nRows; ++i )
          if ( aggStat[i] == SELECTED ) localSelected++;

        // Compute 'globalSelected'
        GO globalSelected; sumAll(comm, (GO)localSelected, globalSelected);

        // Compute 'globalNRows'
        GO globalNRows; sumAll(comm, (GO)nRows, globalNRows);

        GetOStream(Statistics1, 0) << "Nodes aggregated = " << globalSelected << " (" << globalNRows << ")" << std::endl;
      }

      if (IsPrint(Statistics1)) {
        GO nAggregatesGlobal; sumAll(comm, (GO)nAggregates, nAggregatesGlobal);
        GetOStream(Statistics1, 0) << "Total aggregates = " << nAggregatesGlobal << std::endl;
      }

    } // verbose

      /* ------------------------------------------------------------- */
      /* clean up                                                      */
      /* ------------------------------------------------------------- */

    aggCurrent = aggHead;
    while ( aggCurrent != NULL )
      {
        supernode = aggCurrent;
        aggCurrent = aggCurrent->next;
        delete supernode;
      }

  } // CoarsenUncoupled

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomReorder(Teuchos::ArrayRCP<LO> list) const {
    //TODO: replace int
    int n = list.size();
    for(int i=0; i<n-1; i++) {
      std::swap(list[i], list[RandomOrdinal(i,n-1)]);
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  int LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomOrdinal(int min, int max) const {
    return min + static_cast<int>((max-min+1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
  }

} // namespace MueLu

#endif // MUELU_LOCALAGGREGATIONALGORITHM_DEF_HPP
