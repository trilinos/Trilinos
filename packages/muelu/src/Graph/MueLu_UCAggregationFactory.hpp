#ifndef MUELU_UCAGGREGATIONFACTORY_HPP
#define MUELU_UCAGGREGATIONFACTORY_HPP

#include "Cthulhu_VectorFactory.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AggregationOptions.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AggAlgorithm2Comm.hpp"

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <iostream>

namespace MueLu {

using Teuchos::ArrayView;
using Teuchos::ArrayRCP;

// MPI helper
#define sumAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out) \
  Teuchos::reduceAll<int>(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));


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

#define MUELU_NOSCORE       -100  /* indicates that a quality score has not  */
                                  /* yet been assigned when determining to   */
                                  /* which existing aggregate a vertex       */
                                  /* should be assigned.                     */

#define MUELU_DISTONE_VERTEX_WEIGHT 100  /* Weights associated with all      */
                                  /* vertices that have a direct connection  */
                                  /* to the aggregate root.                  */

#define INCR_SCALING 3            /* Determines how much of a penalty should */
                                  /* be deduced from a score during Phase 5  */
                                  /* for each Phase 5 vertex already added   */
                                  /* to this aggregate. Specifically the     */
                                  /* penalty associated with aggregate y is  */
                                  /*   max (INCR_SCALING*NNewVtx,            */
                                  /*        UnpenalizedScore*(1-             */
                                  /*              MUELU_PENALTYFACTOR))      */
                                  /* where NNewVtx is the number of phase 5  */
                                  /* vertices already assigned to y.         */

#define MUELU_PENALTYFACTOR .30   /* Determines maximum allowable            */
                                  /* percentage of a score that can be       */
                                  /* deducted from this score for having     */
                                  /* already enlargened an aggregate to      */
                                  /* which we are contemplated adding another*/
                                  /* vertex.  Should be between 0 and 1.     */


/*!
  @class UCAggregationFactory class.
  @brief Factory for coarsening a graph with uncoupled aggregation.

  This method has two phases.  The first is a local clustering algorithm.  The second creates aggregates
  that can include unknowns from more than one process.
*/

template<class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class UCAggregationFactory : public Teuchos::Describable {
#include "MueLu_UseShortNames_Graph.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    UCAggregationFactory() {}

    //! Destructor.
    virtual ~UCAggregationFactory() {}
    //@}

    //! @name Build methods.
    //@{

    //! Build aggregates.
    Teuchos::RCP<Aggregates> Build(Graph const &graph, AggregationOptions const &options) const;
    //@}

  private:
      //! aggregation algorithm type
      std::string Algorithm_;

      //! @name Aggregation methods.
      //@{

      /*! @brief Local aggregation.
      */
      RCP<MueLu::Aggregates<LO,GO> >
      CoarsenUncoupled(MueLu::AggregationOptions const & aggOptions, MueLu::Graph<LO,GO> const & graph) const;


      /*! @brief Take a partially aggregated graph and complete the aggregation.
      
         This is typically needed to take care of vertices that are left over after
         creating a bunch of ideal aggregates (root plus immediate neighbors).
        
         On input, the structure Aggregates describes already aggregated vertices.
         The field procWinners[] indicates the processor owning the aggregate to
         which a vertex is "definitively" assigned. If on entry 
         procWinners[i] == MUELU_UNASSIGNED, MueLu_ArbitrateAndCommunicate() 
         will arbitrate and decide which processor's aggregate really has
         the vertex. If only one processor claims ownership (as in
         the Uncoupled case), no real arbitration is needed. Otherwise,
         random arbitration is done.
        
         This cleanup has many phases:
           
           Phase 1b: Invoke MueLu_ArbitrateAndCommunicate() to ensure that
                     all processors have the same view of aggregated vertices
                     (e.g., to which aggregate they have been assigend and
                     which processor owns that aggregate).
           Phase 2:  Check for vertices (local or nonlocal) which are Adjacent
                     to root nodes. Tentatively assign these to the aggregate
                     associated with the root. Arbitrate any cases where 
                     several processors claim the same vertex for one of 
                     their aggregates via MueLu_ArbitrateAndCommunicate().
           Phase 3:  Try to create new aggregates if it looks like there are
                     root node candidates which have many unaggregated neighbors.
                     This decision to make a new aggregate is based on only local
                     information. However, the new aggregate will be tentatively
                     assigned any unaggregated ghost vertices. Arbitration is again
                     done by MueLu_ArbitrateAndCommunicate() where local vertices
                     use a weight[] = 2 and ghost vertices have weight[] = 1.
                     The basic idea is that after arbitration, each aggregate
                     is guaranteed to keep all local vertices assigned in
                     this phase. Thus, by basing the aggregation creation logic 
                     on local information, we are guarantee to have a sufficiently
                     large aggregation. The only local vertices that might be
                     assigned to another processor's aggregates are unclaimed
                     during this phase of the aggregation.
           Phase 5:  Sweep new points into existing aggregates. Each processor tries
                     to assign any (whether it is a ghost or local) unaggregated
                     vertex that it has to an aggregate that it owns. In other words,
                     processor p attempts to assign vertex v to aggregate y where
                     y is owned by p but v may be a ghost vertex (and so really 
                     assigned to another processor). Deciding which aggregate
                     a vertex is assigned to is done by scoring. Roughly, we want 
        
                          a) larger scores for y if v is is close (graph distance)
                             to y's root.
                          b) larger scores for y if v has direct connections to 
                             several different vertices already assigned to y.
                          c) lower scores for y if several vertices have already
                             been swept into y during this phase.
        
                     Some care must be taken for vertices that are shared (either
                     local vertices that are sent to other processors or ghost
                     vertices) in that each processor sharing the vertex
                     will attempt to assign it to different aggregates. 
                     MueLu_ArbitrateAndCommunicate() is again used for arbitration
                     with the score being given as the weight. 
        
                     The main tricky thing occurs when v is tentatively added to y.
                     When assigning vprime to y, the assumed connection with v should
                     not be the sole basis of this decisioin if there is some chance
                     that v might be lost in arbitration. This could actually lead to
                     vprime being disconnected from the rest of the aggregate.  This
                     is by building a list of shared ids and checking that there is
                     at least one vertex in the scoring that 
                     is either not shared or has already been definitively
                     assigned to this processor's aggregate (i.e. have been assigned
                     to a local aggregate and have been through arbitration).
                     
                     Scoring is done by first giving a mark to vertices that have been
                     already been assigned to aggregates. This mark essentially
                     reflects the distance of this point to the root. Specifically,
        
                       mark(v) <-- MUELU_DISTONE_VERTEX_WEIGHT if v assigned to 
                                                                aggregate prior
                                                                to this phase.
        
                       mark(v) <-- max(mark(vk))/2              otherwise
        
                     where max(mark(vk)) considers all vertices definitively
                     assigned to y that have direct connections to v.
        
                     Finally,
                       score(vtilde,y)<--sum(mark(vkhat)) - AggregateIncrementPenalty
        
                     where vtilde is an unaggregated vertex being considered for
                     assignment in aggregate y and vkhat are all vertices in y
                     with a direct connection to vtilde. AggregateIncrementPenalty
                     is equal to 
                         max (INCR_SCALING*NNewVtx,
                              sum(mark(vkhat))*(1-MUELU_PENALTYFACTOR))
                     where NNewVtx is the number of phase 5 vertices already
                     assigned to y.
        
                     One last wrinkle, is that we have wrapped the whole 
                     scoring/assigning of vertices around a big loop that
                     looks something like
                        for ( Threshold = big; Threshold >= 0; Reduce(Threshold)){
                                 .
                                 .
                                 .
                             MueLu_ArbitrateAndCommunicate() i
                        }
        
                     New vertices are swept into aggregates only if their best
                     score is >= a Threshold.  This encourages only very good
                     vertices to be assigned first followed by somewhat less
                     well connected ones in later iterations of the loop.
                     It also helps minimize the number of exclusions that would
                     occur to address the issue mentioned above where we don't want
                     to make assignment decisions based on connections to vertices
                     that might be later lost in arbitration.
           Phase 6:  Aggregate remaining vertices and avoid small aggregates (e.g.,
                     singletons) at all costs. Typically, most everything should
                     be aggregated by Phase's 1-5.  One way that we could still have 
                     unaggegated vertices is if processor p was never assigned a
                     root node (perhaps because the number of local vertices on p
                     is less than minNodesPerAggregate) and additionally p has
                     local ids which are not shared with any other processors (so
                     that no other processor's aggregate can claim these vertices).
                     
                     Phase 6 looks at the first unassigned vertex and all of its
                     local unassigned neighbors and makes a new aggregate. If this
                     aggregate has at least minNodesPerAggregate vertices, 
                     we continue this process of creating new aggregates by 
                     examining other unassigned vertices. If the new aggregate
                     is too small, we try add the next unassigned vertex
                     and its neighbors to the same newly created aggregate. 
                     Once again, we check the size of this new aggregate to
                     decide whether other unassigned vertices should be added
                     to this aggregate or used to create a new aggregate. 
                     If the last newly created aggregate (of course there may be just
                     one newly created aggregate) is too small, we then see if
                     there is at least one other aggregate that this processor owns.
                     If so, we merge the small newly created aggregate with aggregate
                     0. If not, we accept the fact that a small aggregate has been
                     created.

          
         One final note about the use of MueLu_ArbitrateAndCommunicate(). No
         arbitration occurs (which means the procWinner[] is not set as well) for a
         global shared id if and only if all weights on all processors corresponding
         to this id is zero. Thus, the general idea is that any global id where we
         want arbitration to occur should have at least one associated weight on 
         one processor which is nonzero. Any global id where we don't want any
         arbitration should have all weights set to 0.
        
         Note: procWinners is also set to MyPid() by MueLu_ArbitrateAndCommunicate()
         for any nonshared gid's with a nonzero weight.
      */
      int AggregateLeftOvers(AggregationOptions const &aggOptions, Aggregates &aggregates,
                             std::string const & label, Graph const &graph) const;

      //@}

      //! @name Utilities
      //@{

      /*! @brief Build a list of candidate root nodes.

         Candidates are vertices not adjacent to already aggregated vertices.
      */
      int MueLu_RootCandidates(int nVertices, ArrayView<const int> & vertex2AggId, const Graph graph,
                               ArrayRCP<int> &candidates, int &nCandidates, int &nCandidatesGlobal) const;

      //! @brief Compute sizes of all the aggregates.
      int MueLu_ComputeAggSizes(Aggregates & aggregates, ArrayRCP<int> & aggSizes) const;

      //! @brief Attempt to clean up aggregates that are too small.
      int MueLu_RemoveSmallAggs(Aggregates& aggregates, int min_size,
                                RCP<Cthulhu::Vector<double> > & distWeights, const AggAlgorithm2Comm & myWidget) const;

      /*! @brief Utility to take a list of integers (which should be the same 
         length as the number of local ids in Map) and reorder them randomly.

         @param list      On input, a bunch of integers. On output, the same integers in a different order
                          that is determined randomly.
         @param map       ?????????????????????????
      */
      int MueLu_RandomReorder(Teuchos::ArrayRCP<int> list, const Map &map) const;

      //@}


      // JG TODO: rename variables:
      //  Adjacent-> adjacent
      //  homogenization of variables names :
      //  - colj and j
      //  - i and iNode
      //  - k->kNode
      //  - ...

}; //class UCAggregationFactory

} //namespace MueLu

#define MUELU_UCAGGREGATIONFACTORY_SHORT

#endif //ifndef MUELU_UCAGGREGATIONFACTORY_HPP
