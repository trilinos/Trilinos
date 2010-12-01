/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to aggregate vertices in a graph                                */
/* ************************************************************************* */
/* Author        : Ray Tuminaro                                              */
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


#include "MueLu_AggOptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Graph.hpp"

using namespace std;
using namespace MueLu;

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

extern int MueLu_RandomReorder(int *randomVector, const Epetra_BlockMap &map);


extern Aggregates *MueLu_Aggregate_CoarsenUncoupled(MueLu_AggOptions *aggregateOptions,
                        MueLu_Graph *graph);

extern int MueLu_AggregateLeftOvers(MueLu_AggOptions *aggregateOptions, 
                                  Aggregates *aggregates,
				  const char *label, MueLu_Graph *graph);

extern int MueLu_NonUnique2NonUnique(const Epetra_Vector &source, 
         Epetra_Vector &dest, const Epetra_Map &uniqueMap, 
         const Epetra_Import &unique2NonUniqueWidget, 
         const Epetra_CombineMode what);

extern int MueLu_NonUnique2NonUnique(const Epetra_IntVector &source, 
         Epetra_IntVector &dest, const Epetra_Map &uniqueMap, 
         const Epetra_Import &unique2NonUniqueWidget, 
         const Epetra_CombineMode what);

int MueLu_ArbitrateAndCommunicate(Epetra_Vector &weight, 
                                  Epetra_IntVector &procWinner,
                                  Epetra_IntVector *companion, const Epetra_Map &uniqueMap, 
                                  const Epetra_Import &unique2NonUniqueWidget, const bool perturb);

extern int MueLu_RootCandidates(int nVertices, int *vertex2AggId, MueLu_Graph *graph,
            int *Candidates, int &NCandidates, int &NCandidatesGlobal);

extern int MueLu_RemoveSmallAggs(Aggregates *aggregates, int min_size,
    Epetra_Vector &weights, const Epetra_Map &uniqueMap, 
    const Epetra_Import &unique2NonUniqueWidget);

extern int MueLu_ComputeAggSizes(Aggregates *aggregates, int *AggSizes);




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

Aggregates *MueLu_Aggregate_CoarsenUncoupled(MueLu_AggOptions 
        *aggregateOptions, MueLu_Graph *graph)
{
   int     i, j, k, m, kk, iNode = 0, jNode, length, nRows;
   int     selectFlag, nAggregates, index, myPid, iNode2;
   int     *vertex2AggId = NULL; //, *itmpArray = NULL,
   int     count;
   int     *aggStat = NULL, ordering;
   double  printFlag;
   int     *randomVector = NULL, *aggCntArray = NULL;
   int     minNodesPerAggregate, maxNeighSelected;
   unsigned int nBytes;
   MueLu_Node       *nodeHead=NULL, *nodeTail=NULL, *newNode=NULL;
   MueLu_SuperNode  *aggHead=NULL, *aggCurrent=NULL, *supernode=NULL;
   Aggregates *aggregates=NULL;

   std::string name = "Uncoupled";

   aggregates = new Aggregates(graph, name.c_str());
   aggregates->GetVertex2AggId()->ExtractView(&vertex2AggId);
   
   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   myPid                   = graph->eGraph->Comm().MyPID();
   minNodesPerAggregate    = aggregateOptions->minNodesPerAggregate;
   maxNeighSelected        = aggregateOptions->maxNeighAlreadySelected;
   ordering                = aggregateOptions->ordering;
   printFlag               = aggregateOptions->printFlag;
   nRows                   = graph->nVertices;

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
      MueLu_RandomReorder(randomVector, graph->eGraph->DomainMap());
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
         length = graph->vertexNeighborsPtr[iNode+1] - 
                  graph->vertexNeighborsPtr[iNode] + 1;
         supernode = (MueLu_SuperNode *) malloc(sizeof(MueLu_SuperNode));      
         supernode->list = (int*) malloc(length*sizeof(int));

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
         for (jNode=graph->vertexNeighborsPtr[iNode];jNode<graph->vertexNeighborsPtr[iNode+1];jNode++) 
         {
            index = graph->vertexNeighbors[jNode];
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
               for (jNode=graph->vertexNeighborsPtr[iNode];jNode<graph->vertexNeighborsPtr[iNode+1];jNode++) 
               {
                  index = graph->vertexNeighbors[jNode];
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
                    for (kk=graph->vertexNeighborsPtr[jNode];kk<graph->vertexNeighborsPtr[jNode+1];kk++) 
                  {
                     if ( aggStat[(graph->vertexNeighbors)[kk]] == MUELOO_AGGR_READY )
                     { 
                        newNode = (MueLu_Node *) malloc(sizeof(MueLu_Node));      
                        newNode->nodeId = graph->vertexNeighbors[kk];
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

   m = 0;
   for ( i = 0; i < nRows; i++ ) 
      if ( aggStat[i] == MUELOO_AGGR_READY ) m++;

   graph->eGraph->Comm().SumAll(&m,&k,1);
   if ( k > 0 && myPid == 0 && printFlag  < MueLu_PrintLevel())
      printf("Aggregation(UC) : Phase 1 (WARNING) - %d READY nodes left\n",k);
   m = 0;
   for ( i = 0; i < nRows; i++ ) 
      if ( aggStat[i] == MUELOO_AGGR_SELECTED ) m++;

   graph->eGraph->Comm().SumAll(&m,&k,1);
   graph->eGraph->Comm().SumAll(&nRows,&m,1);
   graph->eGraph->Comm().SumAll(&nAggregates,&j,1);
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
// Take a partially aggregated graph and complete the aggregation. This is
// typically needed to take care of vertices that are left over after
// creating a bunch of ideal aggregates (root plus immediate neighbors).
//
// On input, the structure Aggregates describes already aggregated vertices.
// The field procWinners[] indicates the processor owning the aggregate to
// which a vertex is "definitively" assigned. If on entry 
// procWinners[i] == MUELOO_UNASSIGNED, MueLu_ArbitrateAndCommunicate() 
// will arbitrate and decide which processor's aggregate really has
// the vertex. If only one processor claims ownership (as in
// the Uncoupled case), no real arbitration is needed. Otherwise,
// random arbitration is done.
//
// This cleanup has many phases:
//   
//   Phase 1b: Invoke MueLu_ArbitrateAndCommunicate() to ensure that
//             all processors have the same view of aggregated vertices
//             (e.g., to which aggregate they have been assigend and
//             which processor owns that aggregate).
//   Phase 2:  Check for vertices (local or nonlocal) which are Adjacent
//             to root nodes. Tentatively assign these to the aggregate
//             associated with the root. Arbitrate any cases where 
//             several processors claim the same vertex for one of 
//             their aggregates via MueLu_ArbitrateAndCommunicate().
//   Phase 3:  Try to create new aggregates if it looks like there are
//             root node candidates which have many unaggregated neighbors.
//             This decision to make a new aggregate is based on only local
//             information. However, the new aggregate will be tentatively
//             assigned any unaggregated ghost vertices. Arbitration is again
//             done by MueLu_ArbitrateAndCommunicate() where local vertices
//             use a weight[] = 2 and ghost vertices have weight[] = 1.
//             The basic idea is that after arbitration, each aggregate
//             is guaranteed to keep all local vertices assigned in
//             this phase. Thus, by basing the aggregation creation logic 
//             on local information, we are guarantee to have a sufficiently
//             large aggregation. The only local vertices that might be
//             assigned to another processor's aggregates are unclaimed
//             during this phase of the aggregation.
//   Phase 5:  Sweep new points into existing aggregates. Each processor tries
//             to assign any (whether it is a ghost or local) unaggregated
//             vertex that it has to an aggregate that it owns. In other words,
//             processor p attempts to assign vertex v to aggregate y where
//             y is owned by p but v may be a ghost vertex (and so really 
//             assigned to another processor). Deciding which aggregate
//             a vertex is assigned to is done by scoring. Roughly, we want 
//
//                  a) larger scores for y if v is is close (graph distance)
//                     to y's root.
//                  b) larger scores for y if v has direct connections to 
//                     several different vertices already assigned to y.
//                  c) lower scores for y if several vertices have already
//                     been swept into y during this phase.
//
//             Some care must be taken for vertices that are shared (either
//             local vertices that are sent to other processors or ghost
//             vertices) in that each processor sharing the vertex
//             will attempt to assign it to different aggregates. 
//             MueLu_ArbitrateAndCommunicate() is again used for arbitration
//             with the score being given as the weight. 
//
//             The main tricky thing occurs when v is tentatively added to y.
//             When assigning vprime to y, the assumed connection with v should
//             not be the sole basis of this decisioin if there is some chance
//             that v might be lost in arbitration. This could actually lead to
//             vprime being disconnected from the rest of the aggregate.  This
//             is by building a list of shared ids and checking that there is
//             at least one vertex in the scoring that 
//             is either not shared or has already been definitively
//             assigned to this processor's aggregate (i.e. have been assigned
//             to a local aggregate and have been through arbitration).
//             
//             Scoring is done by first giving a mark to vertices that have been
//             already been assigned to aggregates. This mark essentially
//             reflects the distance of this point to the root. Specifically,
//
//               mark(v) <-- MUELOO_DISTONE_VERTEX_WEIGHT if v assigned to 
//                                                        aggregate prior
//                                                        to this phase.
//
//               mark(v) <-- max(mark(vk))/2              otherwise
//
//             where max(mark(vk)) considers all vertices definitively
//             assigned to y that have direct connections to v.
//
//             Finally,
//               score(vtilde,y)<--sum(mark(vkhat)) - AggregateIncrementPenalty
//
//             where vtilde is an unaggregated vertex being considered for
//             assignment in aggregate y and vkhat are all vertices in y
//             with a direct connection to vtilde. AggregateIncrementPenalty
//             is equal to 
//                 max (INCR_SCALING*NNewVtx,
//                      sum(mark(vkhat))*(1-MUELOO_PENALTYFACTOR))
//             where NNewVtx is the number of phase 5 vertices already
//             assigned to y.
//
//             One last wrinkle, is that we have wrapped the whole 
//             scoring/assigning of vertices around a big loop that
//             looks something like
//                for ( Threshold = big; Threshold >= 0; Reduce(Threshold)){
//                         .
//                         .
//                         .
//                     MueLu_ArbitrateAndCommunicate() i
//                }
//
//             New vertices are swept into aggregates only if their best
//             score is >= a Threshold.  This encourages only very good
//             vertices to be assigned first followed by somewhat less
//             well connected ones in later iterations of the loop.
//             It also helps minimize the number of exclusions that would
//             occur to address the issue mentioned above where we don't want
//             to make assignment decisions based on connections to vertices
//             that might be later lost in arbitration.
//   Phase 6:  Aggregate remaining vertices and avoid small aggregates (e.g.,
//             singletons) at all costs. Typically, most everything should
//             be aggregated by Phase's 1-5.  One way that we could still have 
//             unaggegated vertices is if processor p was never assigned a
//             root node (perhaps because the number of local vertices on p
//             is less than minNodesPerAggregate) and additionally p has
//             local ids which are not shared with any other processors (so
//             that no other processor's aggregate can claim these vertices).
//             
//             Phase 6 looks at the first unassigned vertex and all of its
//             local unassigned neighbors and makes a new aggregate. If this
//             aggregate has at least minNodesPerAggregate vertices, 
//             we continue this process of creating new aggregates by 
//             examining other unassigned vertices. If the new aggregate
//             is too small, we try add the next unassigned vertex
//             and its neighbors to the same newly created aggregate. 
//             Once again, we check the size of this new aggregate to
//             decide whether other unassigned vertices should be added
//             to this aggregate or used to create a new aggregate. 
//             If the last newly created aggregate (of course there may be just
//             one newly created aggregate) is too small, we then see if
//             there is at least one other aggregate that this processor owns.
//             If so, we merge the small newly created aggregate with aggregate
//             0. If not, we accept the fact that a small aggregate has been
//             created.

//  
// One final note about the use of MueLu_ArbitrateAndCommunicate(). No
// arbitration occurs (which means the procWinner[] is not set as well) for a
// global shared id if and only if all weights on all processors corresponding
// to this id is zero. Thus, the general idea is that any global id where we
// want arbitration to occur should have at least one associated weight on 
// one processor which is nonzero. Any global id where we don't want any
// arbitration should have all weights set to 0.
//
// Note: procWinners is also set to MyPid() by MueLu_ArbitrateAndCommunicate()
// for any nonshared gid's with a nonzero weight.
//

int MueLu_AggregateLeftOvers(MueLu_AggOptions *aggregateOptions, 
                                  Aggregates *aggregates,
				  const char *label,
                                  MueLu_Graph *graph)
{
  int      Nphase1_agg, phase_one_aggregated, i, j, k, kk, nonaggd_neighbors;
  int      AdjacentAgg, total_aggs, *agg_incremented = NULL;
  int      *Mark = NULL, *SumOfMarks = NULL;
  int      best_score, score, best_agg, BestMark, myPid;
  int      nAggregates, *vertex2AggId, nVertices, exp_nRows;
  int      *rowi_col = NULL, rowi_N,total_vertices;
  int      Nleftover = 0, Nsingle = 0, Adjacent;
  double   printFlag, factor = 1., penalty;

  nVertices    = graph->nVertices;
  exp_nRows    = graph->nVertices + graph->nGhost;
  myPid        = graph->eGraph->Comm().MyPID();
  printFlag    = aggregateOptions->printFlag;
  nAggregates  = aggregates->GetNumAggregates();

  int minNodesPerAggregate = aggregateOptions->minNodesPerAggregate;
  Nphase1_agg = nAggregates;

  const Epetra_Map &nonUniqueMap=(Epetra_Map &) aggregates->GetVertex2AggId()->Map();
  const Epetra_Map &uniqueMap=(Epetra_Map &) graph->eGraph->DomainMap();
  const Epetra_Import unique2NonUniqueWidget(nonUniqueMap, uniqueMap);

  // Pull stuff out of epetra vectors

  Epetra_IntVector &Vtx2AggId= (Epetra_IntVector &) *(aggregates->GetVertex2AggId());
  aggregates->GetVertex2AggId()->ExtractView(&vertex2AggId);
  Epetra_IntVector &procWinner = *(aggregates->GetProcWinner());


  Epetra_Vector weights(nonUniqueMap);

  // Aggregated vertices not "definitively" assigned to processors are
  // arbitrated by MueLu_ArbitrateAndCommunicate(). There is some
  // additional logic to prevent losing root nodes in arbitration.
  
  weights.PutScalar(0.);
  for (int i=0;i<nonUniqueMap.NumMyElements();i++) {
     if (procWinner[i] == MUELOO_UNASSIGNED) {
        if (vertex2AggId[i] != MUELOO_UNAGGREGATED) {
           weights[i] = 1.;
           if (aggregates->IsRoot(i)) weights[i] = 2.;
         }
     }
  }
  MueLu_ArbitrateAndCommunicate(weights, procWinner, &Vtx2AggId, 
                uniqueMap, unique2NonUniqueWidget, true);
  weights.PutScalar(0.);//All tentatively assigned vertices are now definitive

  // Tentatively assign any vertex (ghost or local) which neighbors a root
  // to the aggregate associated with the root.

   for (int i =0; i < nVertices; i++) {
     if ( aggregates->IsRoot(i) && (procWinner[i] == myPid) ){
         rowi_col = &(graph->vertexNeighbors[graph->vertexNeighborsPtr[i]]);
         rowi_N   = graph->vertexNeighborsPtr[i+1] - graph->vertexNeighborsPtr[i];
         for (j = 0; j < rowi_N; j++) {
            int colj = rowi_col[j];
            if (vertex2AggId[colj] == MUELOO_UNAGGREGATED) {
               weights[colj]= 1.;
               vertex2AggId[colj] = vertex2AggId[i];
            }
         }
       }
   }
  MueLu_ArbitrateAndCommunicate(weights, procWinner,&Vtx2AggId, 
                uniqueMap, unique2NonUniqueWidget, true);
  weights.PutScalar(0.);//All tentatively assigned vertices are now definitive

  // Record the number of aggregated vertices

  phase_one_aggregated = 0;
  for (i = 0; i < nVertices; i++) {
     if (vertex2AggId[i] != MUELOO_UNAGGREGATED)
         phase_one_aggregated++;
  }
  graph->eGraph->Comm().SumAll(&phase_one_aggregated,&phase_one_aggregated,1);
  graph->eGraph->Comm().SumAll(&nVertices,&total_vertices,1);


   /* Among unaggregated points, see if we can make a reasonable size    */
   /* aggregate out of it. We do this by looking at neighbors and seeing */
   /* how many are unaggregated and on my processor. Loosely,            */
   /* base the number of new aggregates created on the percentage of     */
   /* unaggregated nodes.                                                */

   factor = ((double) phase_one_aggregated)/((double)(total_vertices + 1));
   factor = pow(factor, aggregateOptions->phase3AggCreation);

   for (i = 0; i < nVertices; i++) {
     if (vertex2AggId[i] == MUELOO_UNAGGREGATED) 
     {
       rowi_col = &(graph->vertexNeighbors[graph->vertexNeighborsPtr[i]]);
       rowi_N   = graph->vertexNeighborsPtr[i+1] - graph->vertexNeighborsPtr[i];
       nonaggd_neighbors = 0;
       for (j = 0; j < rowi_N; j++) {
         int colj = rowi_col[j];
         if (vertex2AggId[colj] == MUELOO_UNAGGREGATED && colj < nVertices)
           nonaggd_neighbors++;
       }
       if (  (nonaggd_neighbors > minNodesPerAggregate) &&
          (((double) nonaggd_neighbors)/((double) rowi_N) > factor))
       {
         vertex2AggId[i] = (nAggregates)++;
         for (j = 0; j < rowi_N; j++) {
           int colj = rowi_col[j];
           if (vertex2AggId[colj]==MUELOO_UNAGGREGATED) {
             vertex2AggId[colj] = vertex2AggId[i];
             if (colj < nVertices) weights[colj] = 2.;
             else                  weights[colj] = 1.;
           }
         }
         aggregates->SetIsRoot(i);
         weights[i] = 2.;
       }
     }
   } /*for (i = 0; i < nVertices; i++)*/

  MueLu_ArbitrateAndCommunicate(weights, procWinner,&Vtx2AggId, 
                uniqueMap, unique2NonUniqueWidget, true);
  weights.PutScalar(0.);//All tentatively assigned vertices are now definitive


   if ( printFlag < MueLu_PrintLevel()) {

     graph->eGraph->Comm().SumAll(&Nphase1_agg,&total_aggs,1);
     if (myPid == 0) {
       printf("Aggregation(%s) : Phase 1 - nodes aggregated = %d \n",label,
             phase_one_aggregated);
       printf("Aggregation(%s) : Phase 1 - total aggregates = %d\n",label, total_aggs);
     }
     i = nAggregates - Nphase1_agg;
     graph->eGraph->Comm().SumAll(&i,&i,1);
     if ( myPid == 0 ) {
       printf("Aggregation(%s) : Phase 3 - additional aggregates = %d\n",label, i);
     }
   }

   // Determine vertices that are not shared by setting Temp to all ones
   // and doing MueLu_NonUnique2NonUnique(..., Add). This sums values of all
   // local copies associated with each Gid. Thus, sums > 1 are shared.

   Epetra_Vector temp(nonUniqueMap);
   Epetra_Vector tempOutput(nonUniqueMap);

   temp.PutScalar(1.);  
   tempOutput.PutScalar(0.); 
   MueLu_NonUnique2NonUnique(temp, tempOutput, uniqueMap, 
                             unique2NonUniqueWidget, Add);
   
   vector<bool> gidNotShared(exp_nRows);
   for (int i = 0; i < exp_nRows; i++) {
      if (tempOutput[i] > 1.) gidNotShared[i] = false; 
      else  gidNotShared[i] = true; 
   }

   // Phase 4. 

   double nAggregatesTarget;
   int    nAggregatesGlobal, minNAggs, maxNAggs;

   nAggregatesTarget = ((double) uniqueMap.NumGlobalElements())*
                       (((double) uniqueMap.NumGlobalElements())/
                       ((double)graph->eGraph->NumGlobalNonzeros()));

   graph->eGraph->Comm().SumAll(&nAggregates,&nAggregatesGlobal,1);
   graph->eGraph->Comm().MinAll(&nAggregates,&minNAggs,1);
   graph->eGraph->Comm().MaxAll(&nAggregates,&maxNAggs,1);

   //
   // Only do this phase if things look really bad. THIS
   // CODE IS PRETTY EXPERIMENTAL 
   //
#define MUELOO_PHASE4BUCKETS 6
   if ((nAggregatesGlobal < graph->eGraph->Comm().NumProc()) &&
        (2.5*nAggregatesGlobal < nAggregatesTarget) &&
        (minNAggs ==0) && (maxNAggs <= 1)) {

      Epetra_Util util;
      util.SetSeed( (unsigned int) myPid*2 + (int) (11*rand()));
      k = (int)ceil( (10.*myPid)/graph->eGraph->Comm().NumProc());
      for (i = 0; i < k+7; i++) util.SetSeed( (unsigned int) util.RandomInt() );
      temp.SetSeed( (unsigned int) util.RandomInt() );
      temp.Random(); 

      // build a list of candidate root nodes (vertices not adjacent
      // to aggregated vertices)

      int nCandidates = 0, nCandidatesGlobal;
      int *candidates = new int[nVertices+1];

      double priorThreshold = 0.;
     for (int kkk = 0 ; kkk < MUELOO_PHASE4BUCKETS; kkk++) {
       MueLu_RootCandidates(nVertices, vertex2AggId, graph,
            candidates, nCandidates, nCandidatesGlobal);


       double nTargetNewGuys =  nAggregatesTarget - nAggregatesGlobal;
       double threshold      =  priorThreshold + (1. - priorThreshold)*
                                nTargetNewGuys/(nCandidatesGlobal + .001);
   

       threshold = (threshold*(kkk+1.))/((double) MUELOO_PHASE4BUCKETS);
       priorThreshold = threshold;

       for (k = 0; k < nCandidates ; k++ ) {
         i = candidates[k];
         if ((vertex2AggId[i] == MUELOO_UNAGGREGATED) && 
             (fabs(temp[i])  < threshold)) {
                 // Note: priorThreshold <= fabs(temp[i]) <= 1
            rowi_col = &(graph->vertexNeighbors[graph->vertexNeighborsPtr[i]]);
            rowi_N   = graph->vertexNeighborsPtr[i+1] - graph->vertexNeighborsPtr[i];

            if (rowi_N >= minNodesPerAggregate) {
               int count = 0;
               for (j = 0; j < rowi_N; j++) {
                  Adjacent    = rowi_col[j];
                  // This might not be true if someone close to i
                  // is chosen as a root via fabs(temp[]) < Threshold
                  if (vertex2AggId[Adjacent] == MUELOO_UNAGGREGATED){
                     count++;
                     vertex2AggId[Adjacent] = nAggregates;
                     weights[Adjacent] = 1.;
                  }
               }
               if (count >= minNodesPerAggregate) {
                  vertex2AggId[i] = nAggregates++;
                  weights[i] = 2.;
                  aggregates->SetIsRoot(i);
               }
               else { // undo things
                  for (j = 0; j < rowi_N; j++) {
                     Adjacent    = rowi_col[j];
                     if (vertex2AggId[Adjacent] == nAggregates){
                        vertex2AggId[Adjacent] = MUELOO_UNAGGREGATED;
                        weights[Adjacent] = 0.;
                     }
                  }
               }

            }
         }
       }
       MueLu_ArbitrateAndCommunicate(weights, procWinner,&Vtx2AggId, 
                uniqueMap, unique2NonUniqueWidget, true);
       weights.PutScalar(0.);//All tentatively assigned vertices are now definitive
       graph->eGraph->Comm().SumAll(&nAggregates,&nAggregatesGlobal,1);

     // check that there are no aggregates sizes below minNodesPerAggregate
     
       aggregates->SetNumAggregates(nAggregates);

      MueLu_RemoveSmallAggs(aggregates, minNodesPerAggregate, 
                             weights, uniqueMap, unique2NonUniqueWidget);

      nAggregates = aggregates->GetNumAggregates();
     }   // one possibility
   }

   // Initialize things for Phase 5. This includes building the transpose
   // of the matrix ONLY for transposed rows that correspond to unaggregted
   // ghost vertices. Further, the transpose is only a local transpose. 
   // Nonzero edges which exist on other processors are not represented.

   Mark = (int *) malloc(sizeof(int)* (exp_nRows+1));
   agg_incremented = (int *) malloc(sizeof(int)* (nAggregates+1)); 
   SumOfMarks = (int *) malloc(sizeof(int)*(nAggregates+1));

   for (i = 0; i < exp_nRows; i++)   Mark[i] = MUELOO_DISTONE_VERTEX_WEIGHT;
   for (i = 0; i < nAggregates; i++) agg_incremented[i] = 0;
   for (i = 0; i < nAggregates; i++) SumOfMarks[i] = 0;

   // Grab the transpose matrix graph for unaggregated ghost vertices.
   //     a) count the number of nonzeros per row in the transpose

   vector<int> RowPtr(exp_nRows+1-nVertices);
   for (i = nVertices; i < exp_nRows ;  i++) RowPtr[i-nVertices] = 0;
   for (i = 0; i < nVertices;  i++) {
      rowi_col = &(graph->vertexNeighbors[graph->vertexNeighborsPtr[i]]);
      rowi_N   = graph->vertexNeighborsPtr[i+1] - graph->vertexNeighborsPtr[i];
      for (int k = 0; k < rowi_N; k++) {
         j = rowi_col[k];
         if ( (j >= nVertices) && (vertex2AggId[j] == MUELOO_UNAGGREGATED)){
            RowPtr[j-nVertices]++;
         }
      }
   }

   //     b) Convert RowPtr[i] to point to 1st first nnz spot in row i.

   int iSum = 0, iTemp;
   for (i = nVertices; i < exp_nRows ;  i++) {
      iTemp = RowPtr[i-nVertices];
      RowPtr[i-nVertices] = iSum;
      iSum += iTemp;
   }
   RowPtr[exp_nRows-nVertices] = iSum;
   vector<int> cols(iSum+1);
   
   //     c) Traverse matrix and insert entries in proper location.
   for (i = 0; i < nVertices;  i++) {
      rowi_col = &(graph->vertexNeighbors[graph->vertexNeighborsPtr[i]]);
      rowi_N   = graph->vertexNeighborsPtr[i+1] - graph->vertexNeighborsPtr[i];
      for (int k = 0; k < rowi_N; k++) {
         j = rowi_col[k];
         if ( (j >= nVertices) && (vertex2AggId[j] == MUELOO_UNAGGREGATED)){
            cols[RowPtr[j-nVertices]++] = i;
         }
      }
   }
   //     d) RowPtr[i] points to beginning of row i+1 so shift by one location.

   for (i = exp_nRows; i > nVertices;  i--)
      RowPtr[i-nVertices] = RowPtr[i-1-nVertices];
   RowPtr[0] = 0;
  
   int bestScoreCutoff;
   int thresholds[10] = {300,200,100,50,25,13,7,4,2,0};

   // Stick unaggregated vertices into existing aggregates as described above. 
   
   bool cannotLoseAllFriends; // Used to address possible loss of vertices in 
                              // arbitration of shared nodes discussed above.
   for (kk = 0; kk < 10; kk += 2) {
     bestScoreCutoff = thresholds[kk];
     for (i = 0; i < exp_nRows; i++) {
       if (vertex2AggId[i] == MUELOO_UNAGGREGATED) {

         // Grab neighboring vertices which is either in graph for local ids
         // or sits in transposed fragment just constructed above for ghosts.
         if (i < nVertices) {
            rowi_col = &(graph->vertexNeighbors[graph->vertexNeighborsPtr[i]]);
            rowi_N   = graph->vertexNeighborsPtr[i+1] - graph->vertexNeighborsPtr[i];
         }
         else {
            rowi_col = &(cols[RowPtr[i-nVertices]]);
            rowi_N   = RowPtr[i+1-nVertices] - RowPtr[i-nVertices];
         }
         for (j = 0; j < rowi_N; j++) {
           Adjacent    = rowi_col[j];
           AdjacentAgg = vertex2AggId[Adjacent];

          //Adjacent is aggregated and either I own the aggregate
          // or I could own the aggregate after arbitration.
          if ((AdjacentAgg != MUELOO_UNAGGREGATED) && 
              ((procWinner[Adjacent] == myPid) ||     
               (procWinner[Adjacent] == MUELOO_UNASSIGNED))){
            SumOfMarks[AdjacentAgg] += Mark[Adjacent];
          }
         }
         best_score = MUELOO_NOSCORE;
         for (j = 0; j < rowi_N; j++) {
           Adjacent    = rowi_col[j];
           AdjacentAgg = vertex2AggId[Adjacent];
           //Adjacent is unaggregated, has some value and no
           //other processor has definitively claimed him
           if ((AdjacentAgg != MUELOO_UNAGGREGATED) && 
               (SumOfMarks[AdjacentAgg] != 0) &&
               ((procWinner[Adjacent] == myPid) ||
                (procWinner[Adjacent] == MUELOO_UNASSIGNED ))) {

             // first figure out the penalty associated with
             // AdjacentAgg having already been incremented 
             // during this phase, then compute score.

             penalty = (double) (INCR_SCALING*agg_incremented[AdjacentAgg]);
             if (penalty > MUELOO_PENALTYFACTOR*((double)SumOfMarks[AdjacentAgg]))
               penalty = MUELOO_PENALTYFACTOR*((double)SumOfMarks[AdjacentAgg]);
             score = SumOfMarks[AdjacentAgg]- ((int) floor(penalty));

             if (score > best_score) { 
               best_agg             = AdjacentAgg; 
               best_score           = score;
               BestMark             = Mark[Adjacent];
               cannotLoseAllFriends = false;
 
               // This address issue mentioned above by checking whether
               // Adjacent could be lost in arbitration. weight==0 means that
               // Adjacent was not set during this loop of Phase 5 (and so it
               // has already undergone arbitration). GidNotShared == true 
               // obviously implies that Adjacent cannot be lost to arbitration
               if ((weights[Adjacent]== 0.) || (gidNotShared[Adjacent] == true))
                  cannotLoseAllFriends = true;
             }
             // Another vertex within current best aggregate found.
             // We should have (best_score == score). We need to see
             // if we can improve BestMark and cannotLoseAllFriends.
             else if (best_agg == AdjacentAgg) {
               if ((weights[Adjacent]== 0.) || (gidNotShared[Adjacent] == true))
                  cannotLoseAllFriends = true;
               if (Mark[Adjacent] > BestMark) BestMark = Mark[Adjacent];
             }
           }
         }
         // Clean up
         for (j = 0; j < rowi_N; j++) {
           AdjacentAgg = vertex2AggId[rowi_col[j]];
           if (AdjacentAgg >= 0) SumOfMarks[AdjacentAgg] = 0;
         }
         // Tentatively assign vertex to best_agg. 
         if ( (best_score >= bestScoreCutoff) && (cannotLoseAllFriends)) { 
           vertex2AggId[i] = best_agg;
           weights[i] = best_score;
           agg_incremented[best_agg]++;
           Mark[i] = (int) ceil(   ((double) BestMark)/2.);
         }
       }
     }
     MueLu_ArbitrateAndCommunicate(weights, procWinner,&Vtx2AggId, 
                uniqueMap, unique2NonUniqueWidget, true);
     weights.PutScalar(0.);//All tentatively assigned vertices are now definitive
   }

   // Phase 6: Aggregate remain unaggregated vertices and try at all costs
   //          to avoid small aggregates.
   //          One case where we can find ourselves in this situation
   //          is if all vertices vk adjacent to v have already been
   //          put in other processor's aggregates and v does not have
   //          a direct connection to a local vertex in any of these
   //          aggregates.

   int count = 0;
   for (i = 0; i < nVertices; i++) { 
     if ((vertex2AggId[i] == MUELOO_UNAGGREGATED) ) {
       Nleftover++;
       rowi_col = &(graph->vertexNeighbors[graph->vertexNeighborsPtr[i]]);
       rowi_N   = graph->vertexNeighborsPtr[i+1] - graph->vertexNeighborsPtr[i];

	 // We don't want too small of an aggregate. So lets see if there is an
         // unaggregated neighbor that we can also put with this vertex

         vertex2AggId[i] = nAggregates;
         weights[i] = 1.;
         if (count == 0) aggregates->SetIsRoot(i);
          count++;
          for (j = 0; j < rowi_N; j++) {
            if ((rowi_col[j] != i)&&(vertex2AggId[rowi_col[j]] == MUELOO_UNAGGREGATED)&&
                (rowi_col[j] < nVertices)) {
               vertex2AggId[rowi_col[j]] = nAggregates;
               weights[rowi_col[j]] = 1.;
               count++;
            }
          }
          if ( count >= minNodesPerAggregate) {
             nAggregates++; 
             count = 0;
          }
      }
    }
    // We have something which is under minNodesPerAggregate when 
    if (count != 0) {
       // Can stick small aggregate with 0th aggregate?
       if (nAggregates > 0) {
          for (i = 0; i < nVertices; i++) {
             if ((vertex2AggId[i] == nAggregates) && (procWinner[i] == myPid)){
                vertex2AggId[i] = 0;
                aggregates->SetIsRoot(i,false);
            }
         }
      }
      else {
         Nsingle++;
         nAggregates++;
      }
   }
   MueLu_ArbitrateAndCommunicate(weights, procWinner,&Vtx2AggId, 
                 uniqueMap, unique2NonUniqueWidget, false);


   if (printFlag < MueLu_PrintLevel()) {
     graph->eGraph->Comm().SumAll(&Nsingle,&Nsingle,1);
     graph->eGraph->Comm().SumAll(&Nleftover,&Nleftover,1);
     graph->eGraph->Comm().SumAll(&nAggregates,&total_aggs,1);
     if ( myPid == 0 ) {
       printf("Aggregation(%s) : Phase 3 - total aggregates = %d\n",label, total_aggs);
       printf("Aggregation(%s) : Phase 6 - leftovers = %d and singletons = %d\n",label ,Nleftover, Nsingle);
     }
   }
  if (agg_incremented != NULL) free(agg_incremented);
  if (Mark != NULL) free(Mark);
  if (SumOfMarks != NULL) free(SumOfMarks);
  aggregates->SetNumAggregates(nAggregates);

  return 0;
}

// Utility to take a list of integers (which should be the same 
// length as the number of local ids in Map) and reorder them randomly.
// Input,
//     list[]      A bunch of integers
// Output,    
//     list[]      Same integers as on input but in a different order
//                 that is determined randomly.
//
int MueLu_RandomReorder(int *list, const Epetra_BlockMap &map)
{
   Epetra_Vector     RandVec(map);
   Epetra_IntVector iRandVec(map);

   double *ptr;
   int  *iptr;

   RandVec.Random(); RandVec.ExtractView(&ptr); iRandVec.ExtractView(&iptr);
   for (int i=0; i <  map.NumMyElements(); i++) iptr[i] = (int) (10000.*ptr[i]);
   Epetra_Util::Sort(true,RandVec.Map().NumMyElements(), iptr, 0,NULL,1,&list);
   return 0; 
}

// Redistribute data in source to dest where both source and dest might have 
// multiple copies of the same global id across many processors. The source
// may not have the same value for all of these multiple copies, but on 
// termination dest will have a unique value for each global id.  When multiple
// copies exist in source, 'what' determines how they are combined to make a 
// unique value in dest (see Epetra_CombineMode).
//
//  Input:
//     source                   Vector where multiple copies of some GlobalIds
//                              might exist and might have different values.
//
//     dest                     Allocated but contents ignored.
//
//     uniqueMap                A subset of source.Map() where each GlobalId 
//                              has only one unique copy on one processor.
//                              Normally, source.Map() would have both locals
//                              and ghost elements while uniqueMap would just
//                              have the locals. It should be possible to
//                              remove this or make it an optional argument
//                              and use some existing Epetra capability to 
//                              make a uniqueMap.
//
//     unique2NonUniqueWidget   This corresponds precisely to 
//                                   Epetra_Import unique2NonUniqueWidget(
//                                           source.Map(), uniqueMap);
//                              This could also be eliminated and created
//                              here, but for efficiency user's must pass in.
//
//     what                     Determines how multiple copies of the same
//                              GlobalId are combined (see Epetra_CombineMode).
//
//  Output:
//
//     dest                     dest has redistributed data from source where
//                              'what' determines how multiple copies of source
//                              values associated with the same GlobalId are
//                              combined into a unique value on all processors.
//
int MueLu_NonUnique2NonUnique(const Epetra_Vector &source, 
     Epetra_Vector &dest, const Epetra_Map &uniqueMap, 
     const Epetra_Import &unique2NonUniqueWidget, const Epetra_CombineMode what)
{
  Epetra_Vector temp(uniqueMap);
  temp.Export(source, unique2NonUniqueWidget, what);
  dest.Import(temp,   unique2NonUniqueWidget, Insert);

  return 0;
}

//
// For each GlobalId associated with weight.Map():
//
//      1) find the maximum absolute value of weight[] distributed across all
//         processors and assign this to all local elements of weight[] (across 
//         processors) associated with the GlobalId.
//      2) set procWinner[] to the MyPid() that had the largest element.
//         procWinner[] is still set if only one processor owns a GlobalId. 
//
//         The ONLY CASE when procWinner[i] is NOT set corresponds to when
//         all local weights associated with a GlobalId are zero. This allows
//         one to effectively skip the maximum/winner calculation for a subset
//         of GlobalId's.  This might occur when a processor has already
//         claimed ownership for a GlobalId and so all local copies have
//         the same value. We want to skip the maximum calculation with 
//         tiebreaking to avoid another processor claiming ownership.
//
//      3) optionally, set companion[] (across all relevant processors) to the
//         local companion value associated with the procWinner[] processor.
//
//  Input:
//     weight                   Vector of weights. ASSUMED TO BE nonnegative.
//
//     procWinner               Allocated but contents ignored.
//
//     companion                Either NULL or allocated but contents ignored.
//                              If NULL, step 3 above is skipped.
//
//     uniqueMap                A subset of weight.Map() where each GlobalId 
//                              has only one unique copy on one processor.
//                              Normally, weight.Map() would have both locals
//                              and ghost elements while uniqueMap would just
//                              have the locals. It should be possible to
//                              remove this or make it an optional argument
//                              and use some existing Epetra capability to 
//                              make a uniqueMap.
//
//     unique2NonUniqueWidget   This corresponds precisely to 
//                                   Epetra_Import unique2NonUniqueWidget(
//                                           weight.Map(), uniqueMap);
//                              This could also be eliminated and created
//                              here, but for efficiency user's must pass in.
//
//     perturb                  Optional arguments that is either true or 
//                              false (default: true). weight is perturbed
//                              and the perturbed values are used in step 1)
//                              above. Returned values reflect the perturbed
//                              data. This option avoids having lots of
//                              tiebreaks where the large MyPid() always wins.
//
//  Output:
//     weight                   weight[k] <-- Max(weight[k1],...,weight[kn])
//                              where weight[kj] live on different processors
//                              but have the same GlobalId as weight[k] on
//                              this processor.
//
//     procWinner               procWinner[k] <-- MyPid associated with the
//                              kj yielding the max in 
//                                    Max(weight[k1],...,weight[kn]) .
//                              See weight Output comments.
//                              NOTE: If all input weight[kj]'s are zero,
//                                    then procWinner[k] is left untouched.
//
//     companion                If not null, 
//                                 companion[k] <-- companion[kj] where
//                              companion[kj] lives on processor procWinner[k].
//                              and corresponds to the same GlobalId as k.
//                              NOTE: If for a particlar GlobalId, no processor
//                                    has a value of procWinner that matches
//                                    its MyPid, the corresponding companion
//                                    is not altered.
//
int MueLu_ArbitrateAndCommunicate(Epetra_Vector &weight, 
               Epetra_IntVector &procWinner,
               Epetra_IntVector *companion, const Epetra_Map &uniqueMap, 
               const Epetra_Import &unique2NonUniqueWidget, const bool perturb)
{
   int MyPid = weight.Comm().MyPID();

   if (perturb) {
      Epetra_Vector perturbWt(weight.Map());

      double largestGlobalWeight;
      weight.MaxValue(&largestGlobalWeight);

      Epetra_Util util;
      util.SetSeed( (unsigned int) MyPid*2 + (int) (11*rand()));
      for (int i = 0; i < 10; i++) util.SetSeed( (unsigned int) util.RandomInt() );

      perturbWt.SetSeed( (unsigned int) util.RandomInt() );
      perturbWt.Random(); 

      for (int i=0; i < weight.Map().NumMyElements(); i++) {
         if (weight[i] == 0.) perturbWt[i] = 0.;
         else (perturbWt)[i] = weight[i] + 1.0e-7*largestGlobalWeight*fabs(perturbWt[i]);
      }
      for (int i=0; i < weight.Map().NumMyElements(); i++) weight[i]= perturbWt[i]; 
   }

   // Communicate weights and store results in PostComm (which will be copied
   // back into weights later. When multiple processors have different weights
   // for the same GID, we take the largest weight. After this fragment every
   // processor should have the same value for PostComm[] even when multiple
   // copies of the same Gid are involved.

   Epetra_Vector postComm(weight.Map());
   postComm.PutScalar(0.0);

   MueLu_NonUnique2NonUnique(weight, postComm, uniqueMap, unique2NonUniqueWidget, AbsMax);

   
   // Let every processor know who is the procWinner. For nonunique
   // copies of the same Gid, this corresponds to the procesosr with
   // the highest Wt[]. When several processors have the same positive value
   // for weight[] (which is also the maximum value), the highest proc id
   // is declared the procWinner.
   //
   // Note:This is accomplished by filling a vector with MyPid+1 if weight[k] is
   //      nonzero and PostComm[k]==weight[k]. NonUnique2NonUnique(...,AbsMax)
   //      is invoked to let everyone know the procWinner.
   //      One is then subtracted so that procWinner[i] indicates the 
   //      Pid of the winning processor.
   //      When all weight's for a GID are zero, the associated procWinner's
   //      are left untouched.

   Epetra_Vector candidateWinners(weight.Map());
   candidateWinners.PutScalar(0.0);

   for (int i=0; i < weight.Map().NumMyElements(); i++) {
       if (postComm[i] == weight[i]) candidateWinners[i] = (double) MyPid+1;
   }

   for (int i=0; i < weight.Map().NumMyElements(); i++) weight[i]=postComm[i]; 

   MueLu_NonUnique2NonUnique(candidateWinners, postComm, uniqueMap, unique2NonUniqueWidget, AbsMax);

   // Note: 
   //                      associated CandidateWinners[]
   //    weight[i]!=0  ==> on some proc is equal to its ==> postComm[i]!=0
   //                      MyPid+1.
   //          
   for (int i=0; i < weight.Map().NumMyElements(); i++)  {
      if ( weight[i] != 0.) procWinner[i] = ((int) (postComm[i])) - 1;
   }

   if (companion != NULL) {
      // Now build a new Map, WinnerMap which just consists of procWinners. 
      // This is done by extracting the Gids for Wt, and shoving
      // the subset that correspond to procWinners in MyWinners.
      // WinnerMap is then constructed using MyWinners.
   
      int numMyWinners = 0;
      for (int i = 0; i < weight.Map().NumMyElements(); i++) {
         if (procWinner[i] == MyPid) numMyWinners++;
      }
   
      int *myGids    = new int[weight.Map().NumMyElements()+1];
      int *myWinners = new int[numMyWinners+1];
   
      weight.Map().MyGlobalElements(myGids);
   
      numMyWinners = 0;
      for (int i = 0; i < weight.Map().NumMyElements(); i++) {
         if (procWinner[i] == MyPid)
            myWinners[numMyWinners++] = myGids[i];
      }
      Epetra_Map winnerMap(-1, numMyWinners, myWinners, 0, weight.Comm());
   
      // Pull the Winners out of companion
      //     JustWinners <-- companion[Winners];
   
      Epetra_IntVector justWinners(winnerMap);
      Epetra_Import winnerImport(winnerMap,weight.Map());
      justWinners.Import(*companion, winnerImport, Insert);
   
      // Put the JustWinner values back into companion so that
      // all nonunique copies of the same Gid have the procWinner's
      // version of the companion.
   
      Epetra_Import pushWinners(weight.Map(), winnerMap);
      companion->Import(justWinners, pushWinners, Insert);

      delete [] myWinners;
      delete [] myGids;
   }

   return 0; //TODO
}

// build a list of candidate root nodes (vertices not adjacent to already
// aggregated vertices)
int MueLu_RootCandidates(int nVertices, int *vertex2AggId, MueLu_Graph *graph,
            int *candidates, int &nCandidates, int &nCandidatesGlobal) {

   int rowi_N, *rowi_col, adjacent;
   bool noAggdNeighbors;

   nCandidates = 0;
 
   for (int i = 0; i < nVertices; i++ ) {
      if (vertex2AggId[i] == MUELOO_UNAGGREGATED) {
         noAggdNeighbors = true;
         rowi_col = &(graph->vertexNeighbors[graph->vertexNeighborsPtr[i]]);
         rowi_N   = graph->vertexNeighborsPtr[i+1] - graph->vertexNeighborsPtr[i];

         for (int j = 0; j < rowi_N; j++) {
            adjacent    = rowi_col[j];
            if (vertex2AggId[adjacent] != MUELOO_UNAGGREGATED) 
               noAggdNeighbors = false;
         }
         if (noAggdNeighbors == true) candidates[nCandidates++] = i;
      }
   }
   graph->eGraph->Comm().SumAll(&nCandidates,&nCandidatesGlobal,1);

   return 0;
}

// Compute sizes of all the aggregates.
int MueLu_ComputeAggSizes(Aggregates *aggregates, int *aggSizes)
{
  int *vertex2AggId;
  int nAggregates = aggregates->GetNumAggregates();
  Epetra_IntVector &procWinner = *(aggregates->GetProcWinner());
  int N = procWinner.Map().NumMyElements();
  int myPid = procWinner.Comm().MyPID();

  aggregates->GetVertex2AggId()->ExtractView(&vertex2AggId);
  for (int i = 0; i < nAggregates; i++) aggSizes[i] = 0;
  for (int k = 0; k < N; k++ ) {
     if (procWinner[k] == myPid) aggSizes[vertex2AggId[k]]++;
  }

  return 0; //TODO
}

int MueLu_RemoveSmallAggs(Aggregates *aggregates, int min_size,
    Epetra_Vector &weights, const Epetra_Map &uniqueMap, 
    const Epetra_Import &unique2NonUniqueWidget) 
{
  int nAggregates = aggregates->GetNumAggregates();
  int *AggInfo = new int[nAggregates+1];
  int *vertex2AggId;

  Epetra_IntVector &procWinner = *(aggregates->GetProcWinner());
  int N = procWinner.Map().NumMyElements();
  int myPid = procWinner.Comm().MyPID();

  Epetra_IntVector &Vtx2AggId= (Epetra_IntVector &) *(aggregates->GetVertex2AggId());
  aggregates->GetVertex2AggId()->ExtractView(&vertex2AggId);

  MueLu_ComputeAggSizes(aggregates, AggInfo);

  // Make a list of all aggregates indicating New AggId
  // Use AggInfo array for this.

  int NewNAggs = 0; 
  for (int i = 0; i < nAggregates; i++) {
     if ( AggInfo[i] < min_size) { 
            AggInfo[i] =  MUELOO_UNAGGREGATED;
     }
     else AggInfo[i] = NewNAggs++;
  }

  for (int k = 0; k < N; k++ ) {
     if (procWinner[k] == myPid) {
        if (vertex2AggId[k] !=  MUELOO_UNAGGREGATED) {
           vertex2AggId[k] = AggInfo[vertex2AggId[k]];
           weights[k] = 1.;
        }
        if (vertex2AggId[k] ==  MUELOO_UNAGGREGATED) 
          aggregates->SetIsRoot(k,false);
     }
  }
  nAggregates = NewNAggs;

  MueLu_ArbitrateAndCommunicate(weights, procWinner,&Vtx2AggId, 
                uniqueMap, unique2NonUniqueWidget, true);
  weights.PutScalar(0.);//All tentatively assigned vertices are now definitive

  // procWinner is not set correctly for aggregates which have 
  // been eliminated
  for (int i = 0; i < N; i++) {
     if (vertex2AggId[i] == MUELOO_UNAGGREGATED) 
            procWinner[i] = MUELOO_UNASSIGNED;
  }
  aggregates->SetNumAggregates(nAggregates);

  return 0; //TODO
}
