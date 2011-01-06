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
#include <math.h>
#include <vector>
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Time.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <Cthulhu_VectorFactory.hpp>
#include <Cthulhu_Import.hpp>
//#include <Cthulhu_ConfigDefs.hpp> // CombineMode
#include <Cthulhu_EpetraImport.hpp> //tmp

#include "MueLu_AggregationOptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Graph.hpp"

// MPI helper
#define sumAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

using namespace std;
using namespace MueLu;
using Teuchos::ArrayView;

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
  int    *list;
  struct MueLu_SuperNode_Struct *next;
} MueLu_SuperNode;

//extern int MueLu_RandomReorder(int *randomVector, const Map &map);
int MueLu_RandomReorder(int *randomVector, const Epetra_BlockMap &map); //RRTODO

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
/* jth aggregate. vertex2AggId[k] == MUELOO_UNAGGREGATED indicates that the */
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
  int     i, j, k, m, iNode = 0, jNode, length, nRows;
  int     selectFlag, nAggregates, index, myPid, iNode2;
  Teuchos::ArrayRCP<int> vertex2AggId; //, *itmpArray = NULL,
  int     count;
  int     *aggStat = NULL, ordering;
  double  printFlag;
  int     *randomVector = NULL, *aggCntArray = NULL;
  int     minNodesPerAggregate, maxNeighSelected;
  unsigned int nBytes;
  MueLu_Node       *nodeHead=NULL, *nodeTail=NULL, *newNode=NULL;
  MueLu_SuperNode  *aggHead=NULL, *aggCurrent=NULL, *supernode=NULL;

  std::string name = "Uncoupled";
  RCP<Aggregates<int,int> > aggregates = Teuchos::rcp(new Aggregates<int,int>(graph, name));

  vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);

  /* ============================================================= */
  /* get the machine information and matrix references             */
  /* ============================================================= */

  myPid                   = graph.GetComm()->getRank();
  minNodesPerAggregate    = aggOptions.GetMinNodesPerAggregate();
  maxNeighSelected        = aggOptions.GetMaxNeighAlreadySelected();
  ordering                = aggOptions.GetOrdering();
  printFlag               = aggOptions.GetPrintFlag();
  nRows                   = graph.GetNodeNumVertices();

  /* ============================================================= */
  /* aggStat indicates whether this node has been aggreated, and */
  /* vertex2AggId stores the aggregate number where this node has    */
  /* been aggregated into.                                         */
  /* ============================================================= */

  nBytes = nRows * sizeof( int );
  if (nBytes > 0) aggStat = (int *) malloc(nBytes);
  for ( i = 0; i < nRows; i++ ) aggStat[i] = MUELOO_AGGR_READY;

  /* ============================================================= */
  /* Set up the data structures for aggregation                    */
  /* ============================================================= */

  nAggregates = 0;
  aggHead = NULL;
  nBytes = (nRows+1)*sizeof(int);
  aggCntArray = (int *) malloc(nBytes);
  for ( i = 0; i <= nRows; i++ ) aggCntArray[i] = 0;

  /* ============================================================= */
  /* Phase 1  :                                                    */
  /*    for all nodes, form a new aggregate with its neighbors     */
  /*    if the number of its neighbors having been aggregated does */
  /*    not exceed a given threshold                               */
  /*    (maxNeighSelected = 0 ===> Vanek's scheme)               */
  /* ============================================================= */

  if ( ordering == 1 )       /* random ordering */
    {
      nBytes = nRows * sizeof(int);
      randomVector = (int *) malloc(nBytes);
      for (i = 0; i < nRows; i++) randomVector[i] = i;
//RRTODO      MueLu_RandomReorder(randomVector, *graph.GetDomainMap()); //RTODO
    } 
  else if ( ordering == 2 )  /* graph ordering */
    {
      newNode = (MueLu_Node *) malloc(sizeof(MueLu_Node));      
      newNode->nodeId = 0;
      nodeHead = newNode;
      nodeTail = newNode;
      newNode->next = NULL;
    }
   
  iNode2 = 0;
  while ( iNode2 < nRows)
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
              for ( jNode = 0; jNode < nRows; jNode++ ) 
                {
                  if ( aggStat[jNode] == MUELOO_AGGR_READY )
                    { 
                      newNode = (MueLu_Node *) malloc(sizeof(MueLu_Node));      
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
          free(newNode);
        }

      /*------------------------------------------------------ */
      /* consider further only if the node is in READY mode    */
      /*------------------------------------------------------ */

      if ( aggStat[iNode] == MUELOO_AGGR_READY ) 
        {
          // neighOfINode is the neighbor node list of node 'iNode'.
          ArrayView<const int> neighOfINode = graph.getNeighborVertices(iNode);
          length = neighOfINode.size();
          
          supernode = (MueLu_SuperNode *) malloc(sizeof(MueLu_SuperNode));      
          supernode->list = (int*) malloc((length+1)*sizeof(int));

          if ((supernode->list) == NULL) 
            {
              printf("Error:couldn't allocate memory for supernode! %d\n",
                     length);
              exit(1);
            }

          supernode->maxLength = length;
          supernode->length = 1;
          supernode->list[0] = iNode;
          selectFlag = 1;

          /*--------------------------------------------------- */
          /* count the no. of neighbors having been aggregated  */
          /*--------------------------------------------------- */

          count = 0;
          for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
            {
              index = *it;
              if ( index < nRows ) 
                {
                  if ( aggStat[index] == MUELOO_AGGR_READY || 
                       aggStat[index] == MUELOO_AGGR_NOTSEL ) 
                    supernode->list[supernode->length++] = index;
                  else count++;

                }
            }

          /*--------------------------------------------------- */
          /* if there are too many neighbors aggregated or the  */
          /* number of nodes in the new aggregate is too few,   */
          /* don't do this one                                  */
          /*--------------------------------------------------- */

          if ( count > maxNeighSelected ) selectFlag = 0;

          // Note: the supernode length is actually 1 more than the 
          //       number of nodes in the candidate aggregate. The 
          //       root is counted twice. I'm not sure if this is 
          //       a bug or a feature ... so I'll leave it and change
          //       < to <= in the if just below.

          if (selectFlag != 1 || 
              supernode->length <= minNodesPerAggregate) 
            {
              aggStat[iNode] = MUELOO_AGGR_NOTSEL;
              free( supernode->list );
              free( supernode );
              if ( ordering == 2 ) /* if graph ordering */
                {
                  for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
                    {
                      index = *it;
                      if ( aggStat[index] == MUELOO_AGGR_READY )
                        { 
                          newNode = (MueLu_Node *) malloc(sizeof(MueLu_Node));      
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
              for ( j = 0; j < supernode->length; j++ ) 
                {
                  jNode = supernode->list[j];
                  aggStat[jNode] = MUELOO_AGGR_SELECTED;
                  vertex2AggId[jNode] = nAggregates;
                  if ( ordering == 2 ) /* if graph ordering */
                    {
                      for (ArrayView<const int>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it)
                        {
                          index = *it;
                          if ( aggStat[index] == MUELOO_AGGR_READY )
                            { 
                              newNode = (MueLu_Node *) malloc(sizeof(MueLu_Node));      
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
              aggCntArray[nAggregates++] = supernode->length;
            }
        }
    }
  if ( ordering == 1 ) free(randomVector);
  else if ( ordering == 2 ) 
    {
      while ( nodeHead != NULL )
        {
          newNode = nodeHead;
          nodeHead = newNode->next;
          free( newNode );
        }
    }

  // TODO: replace AllReduce by Reduce to proc 0

  // Compute 'm'
  m = 0;
  for ( i = 0; i < nRows; i++ ) 
    if ( aggStat[i] == MUELOO_AGGR_READY ) m++;

  sumAll(graph.GetComm(), m, k);

  if ( k > 0 && myPid == 0 && printFlag  < MueLu_PrintLevel())
    printf("Aggregation(UC) : Phase 1 (WARNING) - %d READY nodes left\n",k);

  // Compute 'm'
  m = 0;
  for ( i = 0; i < nRows; i++ ) 
    if ( aggStat[i] == MUELOO_AGGR_SELECTED ) m++;

  sumAll(graph.GetComm(), m, k);
  sumAll(graph.GetComm(), nRows, m);
  sumAll(graph.GetComm(), nAggregates, j);

  aggregates->SetNumAggregates(nAggregates);
  if ( myPid == 0 && printFlag  < MueLu_PrintLevel()) 
    {
      printf("Aggregation(UC) : Phase 1 - nodes aggregated = %d (%d)\n",k,m);
      printf("Aggregation(UC) : Phase 1 - total aggregates = %d \n",j);
    }

  /* ------------------------------------------------------------- */
  /* clean up                                                      */
  /* ------------------------------------------------------------- */

  if (aggStat != NULL) free(aggStat);
  free(aggCntArray);
  aggCurrent = aggHead;
  while ( aggCurrent != NULL ) 
    {
      supernode = aggCurrent;
      aggCurrent = aggCurrent->next;
      if ( supernode->maxLength > 0 ) free( supernode->list );
      free( supernode );
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
//int MueLu_RandomReorder(int *list, const Map &map) //RRTODO
int MueLu_RandomReorder(int *list, const Epetra_BlockMap &map)
{

  TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::RuntimeError, "RandomReorder: TODO");

//   Epetra_Vector     RandVec(map);
//   Epetra_IntVector iRandVec(map);

//   double *ptr;
//   int  *iptr;

//   RandVec.Random(); RandVec.ExtractView(&ptr); iRandVec.ExtractView(&iptr);
//   for (int i=0; i <  map.NumMyElements(); i++) iptr[i] = (int) (10000.*ptr[i]);
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
