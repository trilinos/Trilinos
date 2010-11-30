/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */
/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to aggregate vertices in a graph                                */
/* ************************************************************************* */
/* Author        : Ray Tuminaro
/* Date          : Jan, 2010                                              */
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
#include "MueLu_Aggregate.hpp"
#include "MueLu_Graph.hpp"


using namespace std;

int MueLu_PrintLevel() { return 7; }   /* Normally this should be some general*/
                                        /* attribute the indicates the level   */
                                        /* verbosity.                          */

/* ************************************************************************* */
/* linked list structures from ML for holding free node information          */
/* ------------------------------------------------------------------------- */
typedef struct MueLu_Node_Struct
{
   int    node_id;
   struct MueLu_Node_Struct *next;
} MueLu_Node;
/* ************************************************************************* */
/* definition of the structure from ML for holding aggregate information     */
/* ------------------------------------------------------------------------- */
typedef struct MueLu_SuperNode_Struct
{
   int    length;
   int    maxlength;
   int    index;
   int    *list;
   struct MueLu_SuperNode_Struct *next;
} MueLu_SuperNode;

extern int MueLu_RandomReorder(int *randomVector, const Epetra_BlockMap &Map);


extern MueLu_Aggregate *MueLu_Aggregate_CoarsenUncoupled(MueLu_AggOptions *AggregateOptions,
                        MueLu_Graph *Graph);

extern int MueLu_AggregateLeftOvers(MueLu_AggOptions *AggregateOptions, 
                                  MueLu_Aggregate *Aggregates,
				  const char *label, MueLu_Graph *Graph);

extern int MueLu_NonUnique2NonUnique(const Epetra_Vector &source, 
         Epetra_Vector &dest, const Epetra_Map &UniqueMap, 
         const Epetra_Import &Unique2NonUniqueWidget, 
         const Epetra_CombineMode what);

extern int MueLu_NonUnique2NonUnique(const Epetra_IntVector &source, 
         Epetra_IntVector &dest, const Epetra_Map &UniqueMap, 
         const Epetra_Import &Unique2NonUniqueWidget, 
         const Epetra_CombineMode what);

extern int MueLu_ArbitrateAndCommunicate(Epetra_Vector &OrigWt, Epetra_IntVector &ProcWinner,
   Epetra_IntVector *Companion, const Epetra_Map &UniqueMap, const Epetra_Import &Unique2NonUniqueWidget, const bool perturb);

extern int MueLu_RootCandidates(int nvertices, int *Vertex2AggId, MueLu_Graph *Graph,
            int *Candidates, int &NCandidates, int &NCandidatesGlobal);

extern int MueLu_RemoveSmallAggs(MueLu_Aggregate *Aggregates, int min_size,
    Epetra_Vector &Weights, const Epetra_Map &UniqueMap, 
    const Epetra_Import &Unique2NonUniqueWidget);

extern int MueLu_ComputeAggSizes(MueLu_Aggregate *Aggregates, int *AggSizes);




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
/* New aggregates are specified by setting Vertex2AggId. In particular,     */
/* Vertex2AggId[k] = j >= 0 indicates that the kth node resides in the      */
/* jth aggregate. Vertex2AggId[k] == MUELOO_UNAGGREGATED indicates that the */
/* kth node is  unaggregated.                                               */
/*                                                                          */
/* NOTE: This function does not set ProcWinner[]'s. The main trickness      */
/*       with setting ProcWinner[]'s is that one needs to make sure that the*/
/*       ghost ids are properly set on all processors. Since there is no    */
/*       arbitration, this can have easily been done with an import. Instead*/
/*       ProcWinner[] will be set in MueLu_AggregateLeftOvers().           */
/*       MueLu_AggregateLeftOvers() should also function properly when     */
/*       ProcWinner[] is set during Phase 1.                                */
/* ------------------------------------------------------------------------ */

MueLu_Aggregate *MueLu_Aggregate_CoarsenUncoupled(MueLu_AggOptions 
        *AggregateOptions, MueLu_Graph *Graph)
{
   int     i, j, k, m, kk, inode = 0, jnode, length, Nrows;
   int     select_flag, NAggregates, index, mypid, inode2;
   int     *Vertex2AggId = NULL; //, *itmp_array = NULL,
   int     count;
   int     *aggr_stat = NULL, ordering;
   double  printflag;
   int     *randomVector = NULL, *aggr_cnt_array = NULL;
   int     min_nodes_per_aggregate, max_neigh_selected;
   unsigned int nbytes;
   MueLu_Node       *node_head=NULL, *node_tail=NULL, *new_node=NULL;
   MueLu_SuperNode  *aggr_head=NULL, *aggr_curr=NULL, *supernode=NULL;
   MueLu_Aggregate *Aggregates=NULL;

   std::string name = "Uncoupled";

   Aggregates = MueLu_AggregateCreate(Graph, name.c_str());
   Aggregates->Vertex2AggId->ExtractView(&Vertex2AggId);
   
   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid                   = Graph->EGraph->Comm().MyPID();
   min_nodes_per_aggregate = AggregateOptions->min_nodes_per_aggregate;
   max_neigh_selected      = AggregateOptions->max_neigh_already_selected;
   ordering                = AggregateOptions->ordering;
   printflag               = AggregateOptions->print_flag;
   Nrows                   = Graph->NVertices;

   /* ============================================================= */
   /* aggr_stat indicates whether this node has been aggreated, and */
   /* Vertex2AggId stores the aggregate number where this node has    */
   /* been aggregated into.                                         */
   /* ============================================================= */

   nbytes = Nrows * sizeof( int );
   if (nbytes > 0) aggr_stat = (int *) malloc(nbytes);
   for ( i = 0; i < Nrows; i++ ) aggr_stat[i] = MUELOO_AGGR_READY;

   /* ============================================================= */
   /* Set up the data structures for aggregation                    */
   /* ============================================================= */

   NAggregates = 0;
   aggr_head = NULL;
   nbytes = (Nrows+1)*sizeof(int);
   aggr_cnt_array = (int *) malloc(nbytes);
   for ( i = 0; i <= Nrows; i++ ) aggr_cnt_array[i] = 0;

   /* ============================================================= */
   /* Phase 1  :                                                    */
   /*    for all nodes, form a new aggregate with its neighbors     */
   /*    if the number of its neighbors having been aggregated does */
   /*    not exceed a given threshold                               */
   /*    (max_neigh_selected = 0 ===> Vanek's scheme)               */
   /* ============================================================= */

   if ( ordering == 1 )       /* random ordering */
   {
      nbytes = Nrows * sizeof(int);
      randomVector = (int *) malloc(nbytes);
      for (i = 0; i < Nrows; i++) randomVector[i] = i;
      MueLu_RandomReorder(randomVector, Graph->EGraph->DomainMap());
   } 
   else if ( ordering == 2 )  /* graph ordering */
   {
      new_node = (MueLu_Node *) malloc(sizeof(MueLu_Node));      
      new_node->node_id = 0;
      node_head = new_node;
      node_tail = new_node;
      new_node->next = NULL;
   }
   
   inode2 = 0;
   while ( inode2 < Nrows)
   {
      /*------------------------------------------------------ */
      /* pick the next node to aggregate                       */
      /*------------------------------------------------------ */

      if      ( ordering == 0 ) inode = inode2++;
      else if ( ordering == 1 ) inode = randomVector[inode2++];
      else if ( ordering == 2 ) 
      {
         if ( node_head == NULL ) 
         {
            for ( jnode = 0; jnode < Nrows; jnode++ ) 
            {
               if ( aggr_stat[jnode] == MUELOO_AGGR_READY )
               { 
                  new_node = (MueLu_Node *) malloc(sizeof(MueLu_Node));      
                  new_node->node_id = jnode;
                  node_head = new_node;
                  node_tail = new_node;
                  new_node->next = NULL;
                  break;
               }
            }
         }
         if ( node_head == NULL ) break;
         new_node = node_head;
         inode = new_node->node_id;
         node_head = new_node->next;
         free(new_node);
      }

      /*------------------------------------------------------ */
      /* consider further only if the node is in READY mode    */
      /*------------------------------------------------------ */

      if ( aggr_stat[inode] == MUELOO_AGGR_READY ) 
      {
         length = Graph->VertexNeighborsPtr[inode+1] - 
                  Graph->VertexNeighborsPtr[inode] + 1;
         supernode = (MueLu_SuperNode *) malloc(sizeof(MueLu_SuperNode));      
         supernode->list = (int*) malloc(length*sizeof(int));

         if ((supernode->list) == NULL) 
         {
            printf("Error:couldn't allocate memory for supernode! %d\n",
                            length);
            exit(1);
         }

         supernode->maxlength = length;
         supernode->length = 1;
         supernode->list[0] = inode;
         select_flag = 1;

         /*--------------------------------------------------- */
         /* count the no. of neighbors having been aggregated  */
         /*--------------------------------------------------- */

         count = 0;
         for (jnode=Graph->VertexNeighborsPtr[inode];jnode<Graph->VertexNeighborsPtr[inode+1];jnode++) 
         {
            index = Graph->VertexNeighbors[jnode];
            if ( index < Nrows ) 
            {
               if ( aggr_stat[index] == MUELOO_AGGR_READY || 
                    aggr_stat[index] == MUELOO_AGGR_NOTSEL ) 
                  supernode->list[supernode->length++] = index;
               else count++;

            }
         }

         /*--------------------------------------------------- */
         /* if there are too many neighbors aggregated or the  */
         /* number of nodes in the new aggregate is too few,   */
         /* don't do this one                                  */
         /*--------------------------------------------------- */

         if ( count > max_neigh_selected ) select_flag = 0;

         // Note: the supernode length is actually 1 more than the 
         //       number of nodes in the candidate aggregate. The 
         //       root is counted twice. I'm not sure if this is 
         //       a bug or a feature ... so I'll leave it and change
         //       < to <= in the if just below.

         if (select_flag != 1 || 
             supernode->length <= min_nodes_per_aggregate) 
         {
            aggr_stat[inode] = MUELOO_AGGR_NOTSEL;
            free( supernode->list );
            free( supernode );
            if ( ordering == 2 ) /* if graph ordering */
            {
               for (jnode=Graph->VertexNeighborsPtr[inode];jnode<Graph->VertexNeighborsPtr[inode+1];jnode++) 
               {
                  index = Graph->VertexNeighbors[jnode];
                  if ( aggr_stat[index] == MUELOO_AGGR_READY )
                  { 
                     new_node = (MueLu_Node *) malloc(sizeof(MueLu_Node));      
                     new_node->node_id = index;
                     new_node->next = NULL;
                     if ( node_head == NULL )
                     {
                        node_head = new_node;
                        node_tail = new_node;
                     } else {
                        node_tail->next = new_node;
                        node_tail = new_node;
                     }
                  } 
               } 
            } 
         } 
         else 
         {
            Aggregates->IsRoot[inode] = true;
            for ( j = 0; j < supernode->length; j++ ) 
            {
               jnode = supernode->list[j];
               aggr_stat[jnode] = MUELOO_AGGR_SELECTED;
               Vertex2AggId[jnode] = NAggregates;
               if ( ordering == 2 ) /* if graph ordering */
               {
                    for (kk=Graph->VertexNeighborsPtr[jnode];kk<Graph->VertexNeighborsPtr[jnode+1];kk++) 
                  {
                     if ( aggr_stat[(Graph->VertexNeighbors)[kk]] == MUELOO_AGGR_READY )
                     { 
                        new_node = (MueLu_Node *) malloc(sizeof(MueLu_Node));      
                        new_node->node_id = Graph->VertexNeighbors[kk];
                        new_node->next = NULL;
                        if ( node_head == NULL )
                        {
                           node_head = new_node;
                           node_tail = new_node;
                        } else {
                           node_tail->next = new_node;
                           node_tail = new_node;
                        }
                     }
                  } 
               } 
            }
            supernode->next = NULL;
            supernode->index = NAggregates;
            if ( NAggregates == 0 ) 
            {
               aggr_head = supernode;
               aggr_curr = supernode;
            } 
            else 
            {
               aggr_curr->next = supernode;
               aggr_curr = supernode;
            } 
            aggr_cnt_array[NAggregates++] = supernode->length;
         }
      }
   }
   if ( ordering == 1 ) free(randomVector);
   else if ( ordering == 2 ) 
   {
      while ( node_head != NULL )
      {
         new_node = node_head;
         node_head = new_node->next;
         free( new_node );
      }
   }

   m = 0;
   for ( i = 0; i < Nrows; i++ ) 
      if ( aggr_stat[i] == MUELOO_AGGR_READY ) m++;

   Graph->EGraph->Comm().SumAll(&m,&k,1);
   if ( k > 0 && mypid == 0 && printflag  < MueLu_PrintLevel())
      printf("Aggregation(UC) : Phase 1 (WARNING) - %d READY nodes left\n",k);
   m = 0;
   for ( i = 0; i < Nrows; i++ ) 
      if ( aggr_stat[i] == MUELOO_AGGR_SELECTED ) m++;

   Graph->EGraph->Comm().SumAll(&m,&k,1);
   Graph->EGraph->Comm().SumAll(&Nrows,&m,1);
   Graph->EGraph->Comm().SumAll(&NAggregates,&j,1);
   Aggregates->NAggregates = NAggregates;

   if ( mypid == 0 && printflag  < MueLu_PrintLevel()) 
   {
      printf("Aggregation(UC) : Phase 1 - nodes aggregated = %d (%d)\n",k,m);
      printf("Aggregation(UC) : Phase 1 - total aggregates = %d \n",j);
   }

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   if (aggr_stat != NULL) free(aggr_stat);
   free(aggr_cnt_array);
   aggr_curr = aggr_head;
   while ( aggr_curr != NULL ) 
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->maxlength > 0 ) free( supernode->list );
      free( supernode );
   }

   return Aggregates;
}
// Take a partially aggregated graph and complete the aggregation. This is
// typically needed to take care of vertices that are left over after
// creating a bunch of ideal aggregates (root plus immediate neighbors).
//
// On input, the structure Aggregates describes already aggregated vertices.
// The field ProcWinners[] indicates the processor owning the aggregate to
// which a vertex is "definitively" assigned. If on entry 
// ProcWinners[i] == MUELOO_UNASSIGNED, MueLu_ArbitrateAndCommunicate() 
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
//             use a Weight[] = 2 and ghost vertices have Weight[] = 1.
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
//             is less than min_nodes_per_aggregate) and additionally p has
//             local ids which are not shared with any other processors (so
//             that no other processor's aggregate can claim these vertices).
//             
//             Phase 6 looks at the first unassigned vertex and all of its
//             local unassigned neighbors and makes a new aggregate. If this
//             aggregate has at least min_nodes_per_aggregate vertices, 
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
// arbitration occurs (which means the ProcWinner[] is not set as well) for a
// global shared id if and only if all Weights on all processors corresponding
// to this id is zero. Thus, the general idea is that any global id where we
// want arbitration to occur should have at least one associated Weight on 
// one processor which is nonzero. Any global id where we don't want any
// arbitration should have all Weights set to 0.
//
// Note: ProcWinners is also set to MyPID() by MueLu_ArbitrateAndCommunicate()
// for any nonshared gid's with a nonzero weight.
//

int MueLu_AggregateLeftOvers(MueLu_AggOptions *AggregateOptions, 
                                  MueLu_Aggregate *Aggregates,
				  const char *label,
                                  MueLu_Graph *Graph)
{
  int      Nphase1_agg, phase_one_aggregated, i, j, k, kk, nonaggd_neighbors;
  int      AdjacentAgg, total_aggs, *agg_incremented = NULL;
  int      *Mark = NULL, *SumOfMarks = NULL;
  int      best_score, score, best_agg, BestMark, mypid;
  int      NAggregates, *Vertex2AggId, nvertices, exp_Nrows;
  int      *rowi_col = NULL, rowi_N,total_vertices;
  int      Nleftover = 0, Nsingle = 0, Adjacent;
  double   printflag, factor = 1., penalty;

  nvertices    = Graph->NVertices;
  exp_Nrows    = Graph->NVertices + Graph->NGhost;
  mypid        = Graph->EGraph->Comm().MyPID();
  printflag    = AggregateOptions->print_flag;
  NAggregates  = Aggregates->NAggregates;

  int min_nodes_per_aggregate = AggregateOptions->min_nodes_per_aggregate;
  Nphase1_agg = NAggregates;

  const Epetra_Map &NonUniqueMap=(Epetra_Map &) Aggregates->Vertex2AggId->Map();
  const Epetra_Map &UniqueMap=(Epetra_Map &) Graph->EGraph->DomainMap();
  const Epetra_Import Unique2NonUniqueWidget(NonUniqueMap, UniqueMap);

  // Pull stuff out of epetra vectors

  Epetra_IntVector &Vtx2AggId= (Epetra_IntVector &) *(Aggregates->Vertex2AggId);
  Aggregates->Vertex2AggId->ExtractView(&Vertex2AggId);
  Epetra_IntVector &ProcWinner = *(Aggregates->ProcWinner);


  Epetra_Vector Weights(NonUniqueMap);

  // Aggregated vertices not "definitively" assigned to processors are
  // arbitrated by MueLu_ArbitrateAndCommunicate(). There is some
  // additional logic to prevent losing root nodes in arbitration.
  
  Weights.PutScalar(0.);
  for (int i=0;i<NonUniqueMap.NumMyElements();i++) {
     if (ProcWinner[i] == MUELOO_UNASSIGNED) {
        if (Vertex2AggId[i] != MUELOO_UNAGGREGATED) {
           Weights[i] = 1.;
           if ( (Aggregates->IsRoot)[i] == true) Weights[i] = 2.;
         }
     }
  }
  MueLu_ArbitrateAndCommunicate(Weights, ProcWinner, &Vtx2AggId, 
                UniqueMap, Unique2NonUniqueWidget, true);
  Weights.PutScalar(0.);//All tentatively assigned vertices are now definitive

  // Tentatively assign any vertex (ghost or local) which neighbors a root
  // to the aggregate associated with the root.

   for (int i =0; i < nvertices; i++) {
      if ( ( (Aggregates->IsRoot)[i] == true) && (ProcWinner[i] == mypid) ){
         rowi_col = &(Graph->VertexNeighbors[Graph->VertexNeighborsPtr[i]]);
         rowi_N   = Graph->VertexNeighborsPtr[i+1] - Graph->VertexNeighborsPtr[i];
         for (j = 0; j < rowi_N; j++) {
            int colj = rowi_col[j];
            if (Vertex2AggId[colj] == MUELOO_UNAGGREGATED) {
               Weights[colj]= 1.;
               Vertex2AggId[colj] = Vertex2AggId[i];
            }
         }
       }
   }
  MueLu_ArbitrateAndCommunicate(Weights, ProcWinner,&Vtx2AggId, 
                UniqueMap, Unique2NonUniqueWidget, true);
  Weights.PutScalar(0.);//All tentatively assigned vertices are now definitive

  // Record the number of aggregated vertices

  phase_one_aggregated = 0;
  for (i = 0; i < nvertices; i++) {
     if (Vertex2AggId[i] != MUELOO_UNAGGREGATED)
         phase_one_aggregated++;
  }
  Graph->EGraph->Comm().SumAll(&phase_one_aggregated,&phase_one_aggregated,1);
  Graph->EGraph->Comm().SumAll(&nvertices,&total_vertices,1);


   /* Among unaggregated points, see if we can make a reasonable size    */
   /* aggregate out of it. We do this by looking at neighbors and seeing */
   /* how many are unaggregated and on my processor. Loosely,            */
   /* base the number of new aggregates created on the percentage of     */
   /* unaggregated nodes.                                                */

   factor = ((double) phase_one_aggregated)/((double)(total_vertices + 1));
   factor = pow(factor, AggregateOptions->phase3_agg_creation);

   for (i = 0; i < nvertices; i++) {
     if (Vertex2AggId[i] == MUELOO_UNAGGREGATED) 
     {
       rowi_col = &(Graph->VertexNeighbors[Graph->VertexNeighborsPtr[i]]);
       rowi_N   = Graph->VertexNeighborsPtr[i+1] - Graph->VertexNeighborsPtr[i];
       nonaggd_neighbors = 0;
       for (j = 0; j < rowi_N; j++) {
         int colj = rowi_col[j];
         if (Vertex2AggId[colj] == MUELOO_UNAGGREGATED && colj < nvertices)
           nonaggd_neighbors++;
       }
       if (  (nonaggd_neighbors > min_nodes_per_aggregate) &&
          (((double) nonaggd_neighbors)/((double) rowi_N) > factor))
       {
         Vertex2AggId[i] = (NAggregates)++;
         for (j = 0; j < rowi_N; j++) {
           int colj = rowi_col[j];
           if (Vertex2AggId[colj]==MUELOO_UNAGGREGATED) {
             Vertex2AggId[colj] = Vertex2AggId[i];
             if (colj < nvertices) Weights[colj] = 2.;
             else                  Weights[colj] = 1.;
           }
         }
         (Aggregates->IsRoot)[i] = true;
         Weights[i] = 2.;
       }
     }
   } /*for (i = 0; i < nvertices; i++)*/

  MueLu_ArbitrateAndCommunicate(Weights, ProcWinner,&Vtx2AggId, 
                UniqueMap, Unique2NonUniqueWidget, true);
  Weights.PutScalar(0.);//All tentatively assigned vertices are now definitive


   if ( printflag < MueLu_PrintLevel()) {

     Graph->EGraph->Comm().SumAll(&Nphase1_agg,&total_aggs,1);
     if (mypid == 0) {
       printf("Aggregation(%s) : Phase 1 - nodes aggregated = %d \n",label,
             phase_one_aggregated);
       printf("Aggregation(%s) : Phase 1 - total aggregates = %d\n",label, total_aggs);
     }
     i = NAggregates - Nphase1_agg;
     Graph->EGraph->Comm().SumAll(&i,&i,1);
     if ( mypid == 0 ) {
       printf("Aggregation(%s) : Phase 3 - additional aggregates = %d\n",label, i);
     }
   }

   // Determine vertices that are not shared by setting Temp to all ones
   // and doing MueLu_NonUnique2NonUnique(..., Add). This sums values of all
   // local copies associated with each Gid. Thus, sums > 1 are shared.

   Epetra_Vector Temp(NonUniqueMap);
   Epetra_Vector TempOutput(NonUniqueMap);

   Temp.PutScalar(1.);  
   TempOutput.PutScalar(0.); 
   MueLu_NonUnique2NonUnique(Temp, TempOutput, UniqueMap, 
                              Unique2NonUniqueWidget, Add);
   
   vector<bool> GidNotShared(exp_Nrows);
   for (int i = 0; i < exp_Nrows; i++) {
      if (TempOutput[i] > 1.) GidNotShared[i] = false; 
      else  GidNotShared[i] = true; 
   }

   // Phase 4. 

   double NAggregatesTarget;
   int    NAggregatesGlobal, MinNAggs, MaxNAggs;

   NAggregatesTarget = ((double) UniqueMap.NumGlobalElements())*
                       (((double) UniqueMap.NumGlobalElements())/
                       ((double)Graph->EGraph->NumGlobalNonzeros()));

   Graph->EGraph->Comm().SumAll(&NAggregates,&NAggregatesGlobal,1);
   Graph->EGraph->Comm().MinAll(&NAggregates,&MinNAggs,1);
   Graph->EGraph->Comm().MaxAll(&NAggregates,&MaxNAggs,1);

   //
   // Only do this phase if things look really bad. THIS
   // CODE IS PRETTY EXPERIMENTAL 
   //
#define MUELOO_PHASE4BUCKETS 6
   if ((NAggregatesGlobal < Graph->EGraph->Comm().NumProc()) &&
        (2.5*NAggregatesGlobal < NAggregatesTarget) &&
        (MinNAggs ==0) && (MaxNAggs <= 1)) {

      Epetra_Util Util;
      Util.SetSeed( (unsigned int) mypid*2 + (int) (11*rand()));
      k = (int)ceil( (10.*mypid)/Graph->EGraph->Comm().NumProc());
      for (i = 0; i < k+7; i++) Util.SetSeed( (unsigned int) Util.RandomInt() );
      Temp.SetSeed( (unsigned int) Util.RandomInt() );
      Temp.Random(); 

      // build a list of candidate root nodes (vertices not adjacent
      // to aggregated vertices)

      int NCandidates = 0, NCandidatesGlobal;
      int *Candidates = new int[nvertices+1];

      double PriorThreshold = 0.;
     for (int kkk = 0 ; kkk < MUELOO_PHASE4BUCKETS; kkk++) {
       MueLu_RootCandidates(nvertices, Vertex2AggId, Graph,
            Candidates, NCandidates, NCandidatesGlobal);


       double NTargetNewGuys =  NAggregatesTarget - NAggregatesGlobal;
       double Threshold      =  PriorThreshold + (1. - PriorThreshold)*
                                NTargetNewGuys/(NCandidatesGlobal + .001);
   

       Threshold = (Threshold*(kkk+1.))/((double) MUELOO_PHASE4BUCKETS);
       PriorThreshold = Threshold;

       for (k = 0; k < NCandidates ; k++ ) {
         i = Candidates[k];
         if ((Vertex2AggId[i] == MUELOO_UNAGGREGATED) && 
             (fabs(Temp[i])  < Threshold)) {
                 // Note: PriorThreshold <= fabs(Temp[i]) <= 1
            rowi_col = &(Graph->VertexNeighbors[Graph->VertexNeighborsPtr[i]]);
            rowi_N   = Graph->VertexNeighborsPtr[i+1] - Graph->VertexNeighborsPtr[i];

            if (rowi_N >= min_nodes_per_aggregate) {
               int count = 0;
               for (j = 0; j < rowi_N; j++) {
                  Adjacent    = rowi_col[j];
                  // This might not be true if someone close to i
                  // is chosen as a root via fabs(Temp[]) < Threshold
                  if (Vertex2AggId[Adjacent] == MUELOO_UNAGGREGATED){
                     count++;
                     Vertex2AggId[Adjacent] = NAggregates;
                     Weights[Adjacent] = 1.;
                  }
               }
               if (count >= min_nodes_per_aggregate) {
                  Vertex2AggId[i] = NAggregates++;
                  Weights[i] = 2.;
                  (Aggregates->IsRoot)[i] = true;
               }
               else { // undo things
                  for (j = 0; j < rowi_N; j++) {
                     Adjacent    = rowi_col[j];
                     if (Vertex2AggId[Adjacent] == NAggregates){
                        Vertex2AggId[Adjacent] = MUELOO_UNAGGREGATED;
                        Weights[Adjacent] = 0.;
                     }
                  }
               }

            }
         }
       }
       MueLu_ArbitrateAndCommunicate(Weights, ProcWinner,&Vtx2AggId, 
                UniqueMap, Unique2NonUniqueWidget, true);
       Weights.PutScalar(0.);//All tentatively assigned vertices are now definitive
       Graph->EGraph->Comm().SumAll(&NAggregates,&NAggregatesGlobal,1);

     // check that there are no aggregates sizes below min_nodes_per_aggregate
     

      Aggregates->NAggregates = NAggregates;

      MueLu_RemoveSmallAggs(Aggregates, min_nodes_per_aggregate, 
                             Weights, UniqueMap, Unique2NonUniqueWidget);
      NAggregates = Aggregates->NAggregates;
     }   // one possibility
   }

   // Initialize things for Phase 5. This includes building the transpose
   // of the matrix ONLY for transposed rows that correspond to unaggregted
   // ghost vertices. Further, the transpose is only a local transpose. 
   // Nonzero edges which exist on other processors are not represented.

   Mark = (int *) malloc(sizeof(int)* (exp_Nrows+1));
   agg_incremented = (int *) malloc(sizeof(int)* (NAggregates+1)); 
   SumOfMarks = (int *) malloc(sizeof(int)*(NAggregates+1));

   for (i = 0; i < exp_Nrows; i++)   Mark[i] = MUELOO_DISTONE_VERTEX_WEIGHT;
   for (i = 0; i < NAggregates; i++) agg_incremented[i] = 0;
   for (i = 0; i < NAggregates; i++) SumOfMarks[i] = 0;

   // Grab the transpose matrix graph for unaggregated ghost vertices.
   //     a) count the number of nonzeros per row in the transpose

   vector<int> RowPtr(exp_Nrows+1-nvertices);
   for (i = nvertices; i < exp_Nrows ;  i++) RowPtr[i-nvertices] = 0;
   for (i = 0; i < nvertices;  i++) {
      rowi_col = &(Graph->VertexNeighbors[Graph->VertexNeighborsPtr[i]]);
      rowi_N   = Graph->VertexNeighborsPtr[i+1] - Graph->VertexNeighborsPtr[i];
      for (int k = 0; k < rowi_N; k++) {
         j = rowi_col[k];
         if ( (j >= nvertices) && (Vertex2AggId[j] == MUELOO_UNAGGREGATED)){
            RowPtr[j-nvertices]++;
         }
      }
   }

   //     b) Convert RowPtr[i] to point to 1st first nnz spot in row i.

   int isum = 0, itemp;
   for (i = nvertices; i < exp_Nrows ;  i++) {
      itemp = RowPtr[i-nvertices];
      RowPtr[i-nvertices] = isum;
      isum += itemp;
   }
   RowPtr[exp_Nrows-nvertices] = isum;
   vector<int> Cols(isum+1);
   
   //     c) Traverse matrix and insert entries in proper location.
   for (i = 0; i < nvertices;  i++) {
      rowi_col = &(Graph->VertexNeighbors[Graph->VertexNeighborsPtr[i]]);
      rowi_N   = Graph->VertexNeighborsPtr[i+1] - Graph->VertexNeighborsPtr[i];
      for (int k = 0; k < rowi_N; k++) {
         j = rowi_col[k];
         if ( (j >= nvertices) && (Vertex2AggId[j] == MUELOO_UNAGGREGATED)){
            Cols[RowPtr[j-nvertices]++] = i;
         }
      }
   }
   //     d) RowPtr[i] points to beginning of row i+1 so shift by one location.

   for (i = exp_Nrows; i > nvertices;  i--)
      RowPtr[i-nvertices] = RowPtr[i-1-nvertices];
   RowPtr[0] = 0;
  
   int best_score_cutoff;
   int thresholds[10] = {300,200,100,50,25,13,7,4,2,0};

   // Stick unaggregated vertices into existing aggregates as described above. 
   
   bool CannotLoseAllFriends; // Used to address possible loss of vertices in 
                              // arbitration of shared nodes discussed above.
   for (kk = 0; kk < 10; kk += 2) {
     best_score_cutoff = thresholds[kk];
     for (i = 0; i < exp_Nrows; i++) {
       if (Vertex2AggId[i] == MUELOO_UNAGGREGATED) {

         // Grab neighboring vertices which is either in graph for local ids
         // or sits in transposed fragment just constructed above for ghosts.
         if (i < nvertices) {
            rowi_col = &(Graph->VertexNeighbors[Graph->VertexNeighborsPtr[i]]);
            rowi_N   = Graph->VertexNeighborsPtr[i+1] - Graph->VertexNeighborsPtr[i];
         }
         else {
            rowi_col = &(Cols[RowPtr[i-nvertices]]);
            rowi_N   = RowPtr[i+1-nvertices] - RowPtr[i-nvertices];
         }
         for (j = 0; j < rowi_N; j++) {
           Adjacent    = rowi_col[j];
           AdjacentAgg = Vertex2AggId[Adjacent];

          //Adjacent is aggregated and either I own the aggregate
          // or I could own the aggregate after arbitration.
          if ((AdjacentAgg != MUELOO_UNAGGREGATED) && 
              ((ProcWinner[Adjacent] == mypid) ||     
               (ProcWinner[Adjacent] == MUELOO_UNASSIGNED))){
            SumOfMarks[AdjacentAgg] += Mark[Adjacent];
          }
         }
         best_score = MUELOO_NOSCORE;
         for (j = 0; j < rowi_N; j++) {
           Adjacent    = rowi_col[j];
           AdjacentAgg = Vertex2AggId[Adjacent];
           //Adjacent is unaggregated, has some value and no
           //other processor has definitively claimed him
           if ((AdjacentAgg != MUELOO_UNAGGREGATED) && 
               (SumOfMarks[AdjacentAgg] != 0) &&
               ((ProcWinner[Adjacent] == mypid) ||
                (ProcWinner[Adjacent] == MUELOO_UNASSIGNED ))) {

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
               CannotLoseAllFriends = false;
 
               // This address issue mentioned above by checking whether
               // Adjacent could be lost in arbitration. Weight==0 means that
               // Adjacent was not set during this loop of Phase 5 (and so it
               // has already undergone arbitration). GidNotShared == true 
               // obviously implies that Adjacent cannot be lost to arbitration
               if ((Weights[Adjacent]== 0.) || (GidNotShared[Adjacent] == true))
                  CannotLoseAllFriends = true;
             }
             // Another vertex within current best aggregate found.
             // We should have (best_score == score). We need to see
             // if we can improve BestMark and CannotLoseAllFriends.
             else if (best_agg == AdjacentAgg) {
               if ((Weights[Adjacent]== 0.) || (GidNotShared[Adjacent] == true))
                  CannotLoseAllFriends = true;
               if (Mark[Adjacent] > BestMark) BestMark = Mark[Adjacent];
             }
           }
         }
         // Clean up
         for (j = 0; j < rowi_N; j++) {
           AdjacentAgg = Vertex2AggId[rowi_col[j]];
           if (AdjacentAgg >= 0) SumOfMarks[AdjacentAgg] = 0;
         }
         // Tentatively assign vertex to best_agg. 
         if ( (best_score >= best_score_cutoff) && (CannotLoseAllFriends)) { 
           Vertex2AggId[i] = best_agg;
           Weights[i] = best_score;
           agg_incremented[best_agg]++;
           Mark[i] = (int) ceil(   ((double) BestMark)/2.);
         }
       }
     }
     MueLu_ArbitrateAndCommunicate(Weights, ProcWinner,&Vtx2AggId, 
                UniqueMap, Unique2NonUniqueWidget, true);
     Weights.PutScalar(0.);//All tentatively assigned vertices are now definitive
   }

   // Phase 6: Aggregate remain unaggregated vertices and try at all costs
   //          to avoid small aggregates.
   //          One case where we can find ourselves in this situation
   //          is if all vertices vk adjacent to v have already been
   //          put in other processor's aggregates and v does not have
   //          a direct connection to a local vertex in any of these
   //          aggregates.

   int count = 0;
   for (i = 0; i < nvertices; i++) { 
     if ((Vertex2AggId[i] == MUELOO_UNAGGREGATED) ) {
       Nleftover++;
       rowi_col = &(Graph->VertexNeighbors[Graph->VertexNeighborsPtr[i]]);
       rowi_N   = Graph->VertexNeighborsPtr[i+1] - Graph->VertexNeighborsPtr[i];

	 // We don't want too small of an aggregate. So lets see if there is an
         // unaggregated neighbor that we can also put with this vertex

         Vertex2AggId[i] = NAggregates;
         Weights[i] = 1.;
         if (count ==0) Aggregates->IsRoot[i] = true;
         count++;
	 for (j = 0; j < rowi_N; j++) {
	   if ((rowi_col[j] != i)&&(Vertex2AggId[rowi_col[j]] == MUELOO_UNAGGREGATED)&&
               (rowi_col[j] < nvertices)) {
	      Vertex2AggId[rowi_col[j]] = NAggregates;
              Weights[rowi_col[j]] = 1.;
              count++;
	   }
	 }
         if ( count >= min_nodes_per_aggregate) {
            NAggregates++; 
            count = 0;
         }
     }
   }
   // We have something which is under min_nodes_per_aggregate when 
   if (count != 0) {
      // Can stick small aggregate with 0th aggregate?
      if (NAggregates > 0) {
         for (i = 0; i < nvertices; i++) {
            if ((Vertex2AggId[i] == NAggregates) && (ProcWinner[i] == mypid)){
               Vertex2AggId[i] = 0;
               (Aggregates->IsRoot)[i] = false;
            }
         }
      }
      else {
         Nsingle++;
         NAggregates++;
      }
   }
   MueLu_ArbitrateAndCommunicate(Weights, ProcWinner,&Vtx2AggId, 
                 UniqueMap, Unique2NonUniqueWidget, false);


   if (printflag < MueLu_PrintLevel()) {
     Graph->EGraph->Comm().SumAll(&Nsingle,&Nsingle,1);
     Graph->EGraph->Comm().SumAll(&Nleftover,&Nleftover,1);
     Graph->EGraph->Comm().SumAll(&NAggregates,&total_aggs,1);
     if ( mypid == 0 ) {
       printf("Aggregation(%s) : Phase 3 - total aggregates = %d\n",label, total_aggs);
       printf("Aggregation(%s) : Phase 6 - leftovers = %d and singletons = %d\n",label ,Nleftover, Nsingle);
     }
   }
  if (agg_incremented != NULL) free(agg_incremented);
  if (Mark != NULL) free(Mark);
  if (SumOfMarks != NULL) free(SumOfMarks);
  Aggregates->NAggregates = NAggregates;

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
int MueLu_RandomReorder(int *list, const Epetra_BlockMap &Map)
{
   Epetra_Vector     RandVec(Map);
   Epetra_IntVector iRandVec(Map);

   double *ptr;
   int  *iptr;

   RandVec.Random(); RandVec.ExtractView(&ptr); iRandVec.ExtractView(&iptr);
   for (int i=0; i <  Map.NumMyElements(); i++) iptr[i] = (int) (10000.*ptr[i]);
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
//     UniqueMap                A subset of source.Map() where each GlobalId 
//                              has only one unique copy on one processor.
//                              Normally, source.Map() would have both locals
//                              and ghost elements while UniqueMap would just
//                              have the locals. It should be possible to
//                              remove this or make it an optional argument
//                              and use some existing Epetra capability to 
//                              make a UniqueMap.
//
//     Unique2NonUniqueWidget   This corresponds precisely to 
//                                   Epetra_Import Unique2NonUniqueWidget(
//                                           source.Map(), UniqueMap);
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
     Epetra_Vector &dest, const Epetra_Map &UniqueMap, 
     const Epetra_Import &Unique2NonUniqueWidget, const Epetra_CombineMode what)
{
  Epetra_Vector temp(UniqueMap);
  temp.Export(source, Unique2NonUniqueWidget, what);
  dest.Import(temp,   Unique2NonUniqueWidget,Insert);

  return 0;
}

//
// For each GlobalId associated with Weight.Map():
//
//      1) find the maximum absolute value of Weight[] distributed across all
//         processors and assign this to all local elements of Weight[] (across 
//         processors) associated with the GlobalId.
//      2) set ProcWinner[] to the MyPID() that had the largest element.
//         ProcWinner[] is still set if only one processor owns a GlobalId. 
//
//         The ONLY CASE when ProcWinner[i] is NOT set corresponds to when
//         all local weights associated with a GlobalId are zero. This allows
//         one to effectively skip the maximum/winner calculation for a subset
//         of GlobalId's.  This might occur when a processor has already
//         claimed ownership for a GlobalId and so all local copies have
//         the same value. We want to skip the maximum calculation with 
//         tiebreaking to avoid another processor claiming ownership.
//
//      3) optionally, set Companion[] (across all relevant processors) to the
//         local Companion value associated with the ProcWinner[] processor.
//
//  Input:
//     Weight                   Vector of weights. ASSUMED TO BE nonnegative.
//
//     ProcWinner               Allocated but contents ignored.
//
//     Companion                Either NULL or allocated but contents ignored.
//                              If NULL, step 3 above is skipped.
//
//     UniqueMap                A subset of Weight.Map() where each GlobalId 
//                              has only one unique copy on one processor.
//                              Normally, Weight.Map() would have both locals
//                              and ghost elements while UniqueMap would just
//                              have the locals. It should be possible to
//                              remove this or make it an optional argument
//                              and use some existing Epetra capability to 
//                              make a UniqueMap.
//
//     Unique2NonUniqueWidget   This corresponds precisely to 
//                                   Epetra_Import Unique2NonUniqueWidget(
//                                           Weight.Map(), UniqueMap);
//                              This could also be eliminated and created
//                              here, but for efficiency user's must pass in.
//
//     perturb                  Optional arguments that is either true or 
//                              false (default: true). Weight is perturbed
//                              and the perturbed values are used in step 1)
//                              above. Returned values reflect the perturbed
//                              data. This option avoids having lots of
//                              tiebreaks where the large MyPID() always wins.
//
//  Output:
//     Weight                   Weight[k] <-- Max(Weight[k1],...,Weight[kn])
//                              where Weight[kj] live on different processors
//                              but have the same GlobalId as Weight[k] on
//                              this processor.
//
//     ProcWinner               ProcWinner[k] <-- MyPID associated with the
//                              kj yielding the max in 
//                                    Max(Weight[k1],...,Weight[kn]) .
//                              See Weight Output comments.
//                              NOTE: If all input Weight[kj]'s are zero,
//                                    then ProcWinner[k] is left untouched.
//
//     Companion                If not null, 
//                                 Companion[k] <-- Companion[kj] where
//                              Companion[kj] lives on processor ProcWinner[k].
//                              and corresponds to the same GlobalId as k.
//                              NOTE: If for a particlar GlobalId, no processor
//                                    has a value of ProcWinner that matches
//                                    its MyPID, the corresponding Companion
//                                    is not altered.
//
int MueLu_ArbitrateAndCommunicate(Epetra_Vector &Weight, 
               Epetra_IntVector &ProcWinner,
               Epetra_IntVector *Companion, const Epetra_Map &UniqueMap, 
               const Epetra_Import &Unique2NonUniqueWidget, const bool perturb)
{
   int MyPid = Weight.Comm().MyPID();

   if (perturb) {
      Epetra_Vector PerturbWt(Weight.Map());

      double LargestGlobalWeight;
      Weight.MaxValue(&LargestGlobalWeight);

      Epetra_Util Util;
      Util.SetSeed( (unsigned int) MyPid*2 + (int) (11*rand()));
      for (int i = 0; i < 10; i++) Util.SetSeed( (unsigned int) Util.RandomInt() );

      PerturbWt.SetSeed( (unsigned int) Util.RandomInt() );
      PerturbWt.Random(); 

      for (int i=0; i < Weight.Map().NumMyElements(); i++) {
         if (Weight[i] == 0.) PerturbWt[i] = 0.;
         else (PerturbWt)[i] = Weight[i] + 1.0e-7*LargestGlobalWeight*fabs(PerturbWt[i]);
      }
      for (int i=0; i < Weight.Map().NumMyElements(); i++) Weight[i]= PerturbWt[i]; 
   }

   // Communicate Weights and store results in PostComm (which will be copied
   // back into Weights later. When multiple processors have different Weights
   // for the same GID, we take the largest Weight. After this fragment every
   // processor should have the same value for PostComm[] even when multiple
   // copies of the same Gid are involved.

   Epetra_Vector  PostComm(Weight.Map());
   PostComm.PutScalar(0.0);

   MueLu_NonUnique2NonUnique(Weight, PostComm, UniqueMap, Unique2NonUniqueWidget, AbsMax);

   
   // Let every processor know who is the ProcWinner. For nonunique
   // copies of the same Gid, this corresponds to the procesosr with
   // the highest Wt[]. When several processors have the same positive value
   // for Weight[] (which is also the maximum value), the highest proc id
   // is declared the ProcWinner.
   //
   // Note:This is accomplished by filling a vector with MyPid+1 if Weight[k] is
   //      nonzero and PostComm[k]==Weight[k]. NonUnique2NonUnique(...,AbsMax)
   //      is invoked to let everyone know the ProcWinner.
   //      One is then subtracted so that ProcWinner[i] indicates the 
   //      Pid of the winning processor.
   //      When all Weight's for a GID are zero, the associated ProcWinner's
   //      are left untouched.

   Epetra_Vector CandidateWinners(Weight.Map());
   CandidateWinners.PutScalar(0.0);

   for (int i=0; i < Weight.Map().NumMyElements(); i++) {
       if (PostComm[i] == Weight[i]) CandidateWinners[i] = (double) MyPid+1;
   }

   for (int i=0; i < Weight.Map().NumMyElements(); i++) Weight[i]=PostComm[i]; 

   MueLu_NonUnique2NonUnique(CandidateWinners, PostComm, UniqueMap, Unique2NonUniqueWidget, AbsMax);

   // Note: 
   //                      associated CandidateWinners[]
   //    Weight[i]!=0  ==> on some proc is equal to its ==> PostComm[i]!=0
   //                      MyPID+1.
   //          
   for (int i=0; i < Weight.Map().NumMyElements(); i++)  {
      if ( Weight[i] != 0.) ProcWinner[i] = ((int) (PostComm[i])) - 1;
   }

   if (Companion != NULL) {
      // Now build a new Map, WinnerMap which just consists of ProcWinners. 
      // This is done by extracting the Gids for Wt, and shoving
      // the subset that correspond to ProcWinners in MyWinners.
      // WinnerMap is then constructed using MyWinners.
   
      int NumMyWinners = 0;
      for (int i = 0; i < Weight.Map().NumMyElements(); i++) {
         if (ProcWinner[i] == MyPid) NumMyWinners++;
      }
   
      int *MyGids    = new int[Weight.Map().NumMyElements()+1];
      int *MyWinners = new int[NumMyWinners+1];
   
      Weight.Map().MyGlobalElements(MyGids);
   
      NumMyWinners = 0;
      for (int i = 0; i < Weight.Map().NumMyElements(); i++) {
         if (ProcWinner[i] == MyPid)
            MyWinners[NumMyWinners++] = MyGids[i];
      }
      Epetra_Map WinnerMap(-1, NumMyWinners, MyWinners,0,Weight.Comm());
   
      // Pull the Winners out of Companion
      //     JustWinners <-- Companion[Winners];
   
      Epetra_IntVector JustWinners(WinnerMap);
      Epetra_Import WinnerImport(WinnerMap,Weight.Map());
      JustWinners.Import(*Companion, WinnerImport, Insert);
   
      // Put the JustWinner values back into Companion so that
      // all nonunique copies of the same Gid have the ProcWinner's
      // version of the Companion.
   
      Epetra_Import PushWinners(Weight.Map(),WinnerMap);
      Companion->Import(JustWinners, PushWinners, Insert);

      delete [] MyWinners;
      delete [] MyGids;
   }

   return 0; //TODO
}

// build a list of candidate root nodes (vertices not adjacent to already
// aggregated vertices)
int MueLu_RootCandidates(int nvertices, int *Vertex2AggId, MueLu_Graph *Graph,
            int *Candidates, int &NCandidates, int &NCandidatesGlobal) {

   int rowi_N, *rowi_col, Adjacent;
   bool NoAggdNeighbors;

   NCandidates = 0;
 
   for (int i = 0; i < nvertices; i++ ) {
      if (Vertex2AggId[i] == MUELOO_UNAGGREGATED) {
         NoAggdNeighbors = true;
         rowi_col = &(Graph->VertexNeighbors[Graph->VertexNeighborsPtr[i]]);
         rowi_N   = Graph->VertexNeighborsPtr[i+1] - Graph->VertexNeighborsPtr[i];

         for (int j = 0; j < rowi_N; j++) {
            Adjacent    = rowi_col[j];
            if (Vertex2AggId[Adjacent] != MUELOO_UNAGGREGATED) 
               NoAggdNeighbors = false;
         }
         if (NoAggdNeighbors == true) Candidates[NCandidates++] = i;
      }
   }
   Graph->EGraph->Comm().SumAll(&NCandidates,&NCandidatesGlobal,1);

   return 0;
}

// Compute sizes of all the aggregates.
int MueLu_ComputeAggSizes(MueLu_Aggregate *Aggregates, int *AggSizes)
{
  int *Vertex2AggId;
  int NAggregates = Aggregates->NAggregates;
  Epetra_IntVector &ProcWinner = *(Aggregates->ProcWinner);
  int N = ProcWinner.Map().NumMyElements();
  int mypid = ProcWinner.Comm().MyPID();

  Aggregates->Vertex2AggId->ExtractView(&Vertex2AggId);
  for (int i = 0; i < NAggregates; i++) AggSizes[i] = 0;
  for (int k = 0; k < N; k++ ) {
     if (ProcWinner[k] == mypid) AggSizes[Vertex2AggId[k]]++;
  }

  return 0; //TODO
}

int MueLu_RemoveSmallAggs(MueLu_Aggregate *Aggregates, int min_size,
    Epetra_Vector &Weights, const Epetra_Map &UniqueMap, 
    const Epetra_Import &Unique2NonUniqueWidget) 
{
  int NAggregates = Aggregates->NAggregates;
  int *AggInfo = new int[NAggregates+1];
  int *Vertex2AggId;

  Epetra_IntVector &ProcWinner = *(Aggregates->ProcWinner);
  int N = ProcWinner.Map().NumMyElements();
  int mypid = ProcWinner.Comm().MyPID();

  Epetra_IntVector &Vtx2AggId= (Epetra_IntVector &) *(Aggregates->Vertex2AggId);
  Aggregates->Vertex2AggId->ExtractView(&Vertex2AggId);

  MueLu_ComputeAggSizes(Aggregates, AggInfo);

  // Make a list of all aggregates indicating New AggId
  // Use AggInfo array for this.

  int NewNAggs = 0; 
  for (int i = 0; i < NAggregates; i++) {
     if ( AggInfo[i] < min_size) { 
            AggInfo[i] =  MUELOO_UNAGGREGATED;
     }
     else AggInfo[i] = NewNAggs++;
  }

  for (int k = 0; k < N; k++ ) {
     if (ProcWinner[k] == mypid) {
        if (Vertex2AggId[k] !=  MUELOO_UNAGGREGATED) {
           Vertex2AggId[k] = AggInfo[Vertex2AggId[k]];
           Weights[k] = 1.;
        }
        if (Vertex2AggId[k] ==  MUELOO_UNAGGREGATED) 
           (Aggregates->IsRoot)[k] = false;
     }
  }
  NAggregates = NewNAggs;

  MueLu_ArbitrateAndCommunicate(Weights, ProcWinner,&Vtx2AggId, 
                UniqueMap, Unique2NonUniqueWidget, true);
  Weights.PutScalar(0.);//All tentatively assigned vertices are now definitive

  // ProcWinner is not set correctly for aggregates which have 
  // been eliminated
  for (int i = 0; i < N; i++) {
     if (Vertex2AggId[i] == MUELOO_UNAGGREGATED) 
            ProcWinner[i] = MUELOO_UNASSIGNED;
  }
  Aggregates->NAggregates = NAggregates;

  return 0; //TODO
}
