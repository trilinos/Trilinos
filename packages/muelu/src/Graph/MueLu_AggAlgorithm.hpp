/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to aggregate vertices in a graph                                */
/* ************************************************************************* */
/* Inital author : Ray Tuminaro                                              */
/* Date          : Jan, 2010                                                 */
/* ************************************************************************* */

#define CLEAN_DEBUG

#include <assert.h>
#include <stdio.h>

#include <Cthulhu_VectorFactory.hpp>
//#include <Cthulhu_ConfigDefs.hpp> // CombineMode

#include "MueLu_AggregationOptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Graph.hpp"

// MPI helper
#define sumAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));

using namespace std;
using namespace MueLu;
using Teuchos::ArrayView;
typedef ArrayView<const int>::const_iterator iter;

int MueLu_PrintLevel() { return 7; }    /* Normally this should be some general*/
                                        /* attribute the indicates the level   */
                                        /* verbosity.                          */

/* ************************************************************************* */
/* linked list structures from ML for holding free node information          */
/* ------------------------------------------------------------------------- */
typedef struct MueLu_Node_Struct
{
  int    nodeId;
  struct MueLu_Node_Struct *next;
} MueLu_Node;
/* ************************************************************************* */
/* definition of the structure from ML for holding aggregate information     */
/* ------------------------------------------------------------------------- */
typedef struct MueLu_SuperNode_Struct
{
  int    length;
  int    maxLength;
  int    index;
  Teuchos::ArrayRCP<int> list;
  struct MueLu_SuperNode_Struct *next;
} MueLu_SuperNode;

/* In the algorithm, aggStat[]=READY/NOTSEL/SELECTED indicates whether a node has been aggreated. */
enum NodeState {
  READY   = -11,   /* indicates that a node is available to be */
                   /* selected as a root node of an aggregate  */

  NOTSEL  = -12,   /* indicates that a node has been rejected  */
                   /* as a root node. This could perhaps be    */
                   /* because if this node had been selected a */
                   /* small aggregate would have resulted.     */

  SELECTED = -13   /* indicates that a node has been assigned  */
                   /* to an aggregate.                         */
};

int MueLu_RandomReorder(Teuchos::ArrayRCP<int> randomVector, const Map &map);

/* ************************************************************************ */
/* Coarsen Graph by aggregation. In particular, each processor works        */
/* independently on its subgraph (ignoring parts of the graph on other      */
/* processors). The basic idea is to grab an unaggregated node which is not */
/* Adjacent to an already aggregated node. A new aggregate is formed by     */
/* taking this root node and all of its immediate neighbors. The order in   */
/* which nodes are considered as candidate root nodes is an input option.   */
/* On termination, a set of aggregates is created on each processor. There  */
/* are also a set of nodes which remain unaggregated. Normally, these       */
/* unaggregated nodes should be Adjacented to at least one aggregated node  */
/* (otherwise it would have been chosen as a root node).                    */
/*                                                                          */
/* New aggregates are specified by setting vertex2AggId. In particular,     */
/* vertex2AggId[k] = j >= 0 indicates that the kth node resides in the      */
/* jth aggregate. vertex2AggId[k] == MUELU_UNAGGREGATED indicates that the */
/* kth node is  unaggregated.                                               */
/*                                                                          */
/* NOTE: This function does not set procWinner[]'s. The main trickness      */
/*       with setting procWinner[]'s is that one needs to make sure that the*/
/*       ghost ids are properly set on all processors. Since there is no    */
/*       arbitration, this can have easily been done with an import. Instead*/
/*       procWinner[] will be set in MueLu_AggregateLeftOvers().           */
/*       MueLu_AggregateLeftOvers() should also function properly when     */
/*       procWinner[] is set during Phase 1.                                */
/* ------------------------------------------------------------------------ */

RCP<Aggregates<int,int> > MueLu_Aggregate_CoarsenUncoupled(const AggregationOptions & aggOptions, const Graph<int,int> & graph)
{
  /* Create Aggregation object */
  const std::string name = "Uncoupled";
  int nAggregates = 0;
  RCP<Aggregates<int,int> > aggregates = Teuchos::rcp(new Aggregates<int,int>(graph, name));

  /* ============================================================= */
  /* aggStat indicates whether this node has been aggreated, and   */
  /* vertex2AggId stores the aggregate number where this node has  */
  /* been aggregated into.                                         */
  /* ============================================================= */

  Teuchos::ArrayRCP<NodeState> aggStat;
  const int nRows = graph.GetNodeNumVertices();
  if (nRows > 0) aggStat = Teuchos::arcp<NodeState>(nRows);
  for ( int i = 0; i < nRows; ++i ) aggStat[i] = READY;

  /* unused */
  // Teuchos::ArrayRCP<int> aggCntArray = Teuchos::arcp<int>(nRows+1);
  // for ( int i = 0; i <= nRows; ++i ) aggCntArray[i] = 0;

  /* ============================================================= */
  /* Phase 1  :                                                    */
  /*    for all nodes, form a new aggregate with its neighbors     */
  /*    if the number of its neighbors having been aggregated does */
  /*    not exceed a given threshold                               */
  /*    (aggOptions.GetMaxNeighAlreadySelected() = 0 ===> Vanek's scheme) */
  /* ============================================================= */

  /* some general variable declarations */   
  const int ordering = aggOptions.GetOrdering();
  Teuchos::ArrayRCP<int> randomVector;
  MueLu_Node       *nodeHead=NULL, *nodeTail=NULL, *newNode=NULL;
  MueLu_SuperNode  *aggHead=NULL, *aggCurrent=NULL, *supernode=NULL;
  /**/

  if ( ordering == 1 )       /* random ordering */
    {
      randomVector = Teuchos::arcp<int>(nRows);
      for (int i = 0; i < nRows; ++i) randomVector[i] = i;
      MueLu_RandomReorder(randomVector, *graph.GetDomainMap());
    } 
  else if ( ordering == 2 )  /* graph ordering */
    {
      newNode = new MueLu_Node;      
      newNode->nodeId = 0;
      nodeHead = newNode;
      nodeTail = newNode;
      newNode->next = NULL;
    }

  /* main loop */
  {
    int iNode  = 0;
    int iNode2 = 0;
    
    Teuchos::ArrayRCP<int> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0); // output only: contents ignored
    
    while (iNode2 < nRows)
      {
        /*------------------------------------------------------ */
        /* pick the next node to aggregate                       */
        /*------------------------------------------------------ */

        if      ( ordering == 0 ) iNode = iNode2++;
        else if ( ordering == 1 ) iNode = randomVector[iNode2++];
        else if ( ordering == 2 ) 
          {
            if ( nodeHead == NULL ) 
              {
                for ( int jNode = 0; jNode < nRows; ++jNode ) 
                  {
                    if ( aggStat[jNode] == READY )
                      { 
                        newNode = new MueLu_Node;
                        newNode->nodeId = jNode;
                        nodeHead = newNode;
                        nodeTail = newNode;
                        newNode->next = NULL;
                        break;
                      }
                  }
              }
            if ( nodeHead == NULL ) break;
            newNode = nodeHead;
            iNode = newNode->nodeId;
            nodeHead = newNode->next;
            delete newNode;
          }

        /*------------------------------------------------------ */
        /* consider further only if the node is in READY mode    */
        /*------------------------------------------------------ */

        if ( aggStat[iNode] == READY ) 
          {
            // neighOfINode is the neighbor node list of node 'iNode'.
            ArrayView<const int> neighOfINode = graph.getNeighborVertices(iNode);
            int length = neighOfINode.size();
          
            supernode = new MueLu_SuperNode;
            try {
              supernode->list = Teuchos::arcp<int>(length+1);
            } catch (std::bad_alloc&) {
              printf("Error:couldn't allocate memory for supernode! %d\n", length);
              exit(1);
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
              for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
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
            
              if ( count > aggOptions.GetMaxNeighAlreadySelected() ) selectFlag = 0;
            }

            // Note: the supernode length is actually 1 more than the 
            //       number of nodes in the candidate aggregate. The 
            //       root is counted twice. I'm not sure if this is 
            //       a bug or a feature ... so I'll leave it and change
            //       < to <= in the if just below.

            if (selectFlag != 1 || 
                supernode->length <= aggOptions.GetMinNodesPerAggregate()) 
              {
                aggStat[iNode] = NOTSEL;
                delete supernode;
                if ( ordering == 2 ) /* if graph ordering */
                  {
                    for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
                      {
                        int index = *it;
                        if ( aggStat[index] == READY )
                          { 
                            newNode = new MueLu_Node;
                            newNode->nodeId = index;
                            newNode->next = NULL;
                            if ( nodeHead == NULL )
                              {
                                nodeHead = newNode;
                                nodeTail = newNode;
                              } else {
                              nodeTail->next = newNode;
                              nodeTail = newNode;
                            }
                          } 
                      } 
                  } 
              } 
            else 
              {
                aggregates->SetIsRoot(iNode);
                for ( int j = 0; j < supernode->length; ++j ) 
                  {
                    int jNode = supernode->list[j];
                    aggStat[jNode] = SELECTED;
                    vertex2AggId[jNode] = nAggregates;
                    if ( ordering == 2 ) /* if graph ordering */
                      {
                        for (iter it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
                          {
                            int index = *it;
                            if ( aggStat[index] == READY )
                              { 
                                newNode = new MueLu_Node;
                                newNode->nodeId = index;
                                newNode->next = NULL;
                                if ( nodeHead == NULL )
                                  {
                                    nodeHead = newNode;
                                    nodeTail = newNode;
                                  } else {
                                  nodeTail->next = newNode;
                                  nodeTail = newNode;
                                }
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

  if ( ordering == 2 ) 
    {
      while ( nodeHead != NULL )
        {
          newNode = nodeHead;
          nodeHead = newNode->next;
          delete newNode;
        }
    }
  
  /* Update aggregate object */  
  aggregates->SetNumAggregates(nAggregates);

  /* Verbose */
  // TODO: replace AllReduce by Reduce to proc 0
  int myPid = graph.GetComm()->getRank();
  if ( myPid == 0 && aggOptions.GetPrintFlag() < MueLu_PrintLevel()) {
    const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();

    {
      int localReady=0, globalReady;
      
      // Compute 'localReady'
      for ( int i = 0; i < nRows; ++i ) 
        if ( aggStat[i] == READY ) localReady++;
      
      // Compute 'globalReady'
      sumAll(comm, localReady, globalReady);
      
      if (globalReady > 0)
        printf("Aggregation(UC) : Phase 1 (WARNING) - %d READY nodes left\n",globalReady);
    }
    
    {
      int localSelected=0, globalSelected;
      int globalNRows;
      
      // Compute 'localSelected'
      for ( int i = 0; i < nRows; ++i ) 
          if ( aggStat[i] == SELECTED ) localSelected++;
      
      // Compute 'globalSelected' and 'globalNRows'
      sumAll(comm, localSelected, globalSelected);
      sumAll(comm, nRows, globalNRows);
      
      printf("Aggregation(UC) : Phase 1 - nodes aggregated = %d (%d)\n",globalSelected, globalNRows);
    }
    
    {
      int nAggregatesGlobal; 
      sumAll(comm, nAggregates, nAggregatesGlobal);
      printf("Aggregation(UC) : Phase 1 - total aggregates = %d \n",nAggregatesGlobal);
    }
    
  } // if myPid == 0 ...
  
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

  return aggregates;
}

// Utility to take a list of integers (which should be the same 
// length as the number of local ids in Map) and reorder them randomly.
// Input,
//     list[]      A bunch of integers
// Output,    
//     list[]      Same integers as on input but in a different order
//                 that is determined randomly.
//
int MueLu_RandomReorder(Teuchos::ArrayRCP<int> list, const Map &map)
{

  TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "RandomReorder: TODO");

//   Epetra_Vector     RandVec(map);
//   Epetra_IntVector iRandVec(map);

//   double *ptr;
//   int  *iptr;

//   RandVec.Random(); RandVec.ExtractView(&ptr); iRandVec.ExtractView(&iptr);
//   for (int i=0; i <  map.NumMyElements(); ++i) iptr[i] = (int) (10000.*ptr[i]);
//   Epetra_Util::Sort(true,RandVec.getMap().NumMyElements(), iptr, 0,NULL,1,&list);

  return 0; 
}

// JG TODO: rename variables:
//  Adjacent-> adjacent
//  homogenization of variables names :
//  - colj and j
//  - i and iNode
//  - k->kNode
//  - ...
