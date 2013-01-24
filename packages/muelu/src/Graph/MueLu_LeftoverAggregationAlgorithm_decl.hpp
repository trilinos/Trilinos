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
#ifndef MUELU_LEFTOVERAGGREGATIONALGORITHM_DECL_HPP
#define MUELU_LEFTOVERAGGREGATIONALGORITHM_DECL_HPP

#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_LeftoverAggregationAlgorithm_fwd.hpp"

#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_GraphBase.hpp"

#include "MueLu_CoupledAggregationCommHelper_fwd.hpp"

namespace MueLu {

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class LeftoverAggregationAlgorithm : public BaseClass {
#undef MUELU_LEFTOVERAGGREGATIONALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

    typedef GO global_size_t; //TODO
    typedef LO my_size_t;     //TODO

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    LeftoverAggregationAlgorithm();

    //! Destructor.
    virtual ~LeftoverAggregationAlgorithm() { }

    //@}

    //! @name Set/get methods.
    //@{

    void SetMinNodesPerAggregate(int minNodesPerAggregate) { minNodesPerAggregate_ = minNodesPerAggregate; }
    void SetPhase3AggCreation(double phase3AggCreation) { phase3AggCreation_ = phase3AggCreation; }

    double GetPhase3AggCreation() const { return phase3AggCreation_; }
    int GetMinNodesPerAggregate() const { return minNodesPerAggregate_; }

    // TODO: Set/GetGraphName
    //@}

    //! @name Aggregation methods.
    //@{

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


    /*! @brief Take a partially aggregated graph and complete the aggregation.

    This is typically needed to take care of vertices that are left over after
    creating a bunch of ideal aggregates (root plus immediate neighbors).

    On input, the structure Aggregates describes already aggregated vertices.
    The field procWinners[] indicates the processor owning the aggregate to
    which a vertex is "definitively" assigned. If on entry
    procWinners[i] == MUELU_UNASSIGNED, ArbitrateAndCommunicate()
    will arbitrate and decide which processor's aggregate really has
    the vertex. If only one processor claims ownership (as in
    the Uncoupled case), no real arbitration is needed. Otherwise,
    random arbitration is done.

    This cleanup has many phases:

     - Phase 1b: Invoke ArbitrateAndCommunicate() to ensure that
    all processors have the same view of aggregated vertices
    (e.g., to which aggregate they have been assigned and
    which processor owns that aggregate).

     - Phase 2:  Check for vertices (local or nonlocal) which are Adjacent
    to root nodes. Tentatively assign these to the aggregate
    associated with the root. Arbitrate any cases where
    several processors claim the same vertex for one of
    their aggregates via ArbitrateAndCommunicate().

    - Phase 3:  Try to create new aggregates if it looks like there are
    root node candidates which have many unaggregated neighbors.
    This decision to make a new aggregate is based on only local
    information. However, the new aggregate will be tentatively
    assigned any unaggregated ghost vertices. Arbitration is again
    done by ArbitrateAndCommunicate() where local vertices
    use a weight[] = 2 and ghost vertices have weight[] = 1.
    The basic idea is that after arbitration, each aggregate
    is guaranteed to keep all local vertices assigned in
    this phase. Thus, by basing the aggregation creation logic
    on local information, we are guaranteed to have a sufficiently
    large aggregation. The only local vertices that might be
    assigned to another processor's aggregates are unclaimed
    during this phase of the aggregation.

    - Phase 4: EXPERIMENTAL

    - Phase 5:  Sweep new points into existing aggregates. Each processor tries
    to assign any (whether it is a ghost or local) unaggregated
    vertex that it has to an aggregate that it owns. In other words,
    processor p attempts to assign vertex \f$v\f$ to aggregate \f$y\f$ where
    \f$y\f$ is owned by p but \f$v\f$ may be a ghost vertex (and so really
    assigned to another processor). Deciding which aggregate
    a vertex is assigned to is done by scoring. Roughly, we want

     -a) larger scores for \f$y\f$ if \f$v\f$ is is close (graph distance)
    to \f$y\f$'s root.
     -b) larger scores for \f$y\f$ if \f$v\f$ has direct connections to
    several different vertices already assigned to \f$y\f$.
     -c) lower scores for \f$y\f$ if several vertices have already
    been swept into \f$y\f$ during this phase.

    Some care must be taken for vertices that are shared (either
    local vertices that are sent to other processors or ghost
    vertices) in that each processor sharing the vertex
    will attempt to assign it to different aggregates.
    ArbitrateAndCommunicate() is again used for arbitration
    with the score being given as the weight.

    The main tricky thing occurs when \f$v\f$ is tentatively added to \f$y\f$.
    When assigning \f$v'\f$ to \f$y\f$, the assumed connection with \f$v\f$ should
    not be the sole basis of this decision if there is some chance
    that \f$v\f$ might be lost in arbitration. This could actually lead to
    \f$v'\f$ being disconnected from the rest of the aggregate.  This
    is by building a list of shared ids and checking that there is
    at least one vertex in the scoring that
    is either not shared or has already been definitively
    assigned to this processor's aggregate (i.e. have been assigned
    to a local aggregate and have been through arbitration).

    Scoring is done by first giving a mark to vertices that have been
    already been assigned to aggregates. This mark essentially
    reflects the distance of this point to the root. Specifically,

    \f[
    mark(v) \leftarrow MUELU\_DISTONE\_VERTEX\_WEIGHT
    \f]
    
    if \f$v\f$ was assigned to an aggregate prior to this phase,

    \f[
    mark(v) \leftarrow max(mark(v_k))/2
    \f]

    otherwise, where \f$max(mark(v_k))\f$ considers all vertices definitively
    assigned to \f$y\f$ that have direct connections to \f$v\f$.

    Finally,

    \f[
    score(\tilde{v},y) \leftarrow \Sigma(mark(\hat{v}_k)) - AggregateIncrementPenalty
    \f]

    where \f$\tilde{v}\f$ is an unaggregated vertex being considered for
    assignment in aggregate \f$y\f$ and \f$hat{v}_k\f$ are all vertices in \f$y\f$
    with a direct connection to \f$\tilde{v}\f$. AggregateIncrementPenalty
    is equal to

    \f[
    \max (\mbox{INCR_SCALING}*NNewVtx, \Sigma(mark(\hat{v}_k))*(1-\mbox{MUELU_PENALTYFACTOR}))
    \f]

    where \f$ NNewVtx \f$ is the number of phase 5 vertices already
    assigned to \f$y\f$.

    One last wrinkle, is that we have wrapped the whole
    scoring/assigning of vertices around a big loop that
    looks something like

    for ( Threshold = big; Threshold >= 0; Reduce(Threshold));

    New vertices are swept into aggregates only if their best
    score is >= a Threshold.  This encourages only very good
    vertices to be assigned first followed by somewhat less
    well connected ones in later iterations of the loop.
    It also helps minimize the number of exclusions that would
    occur to address the issue mentioned above where we don't want
    to make assignment decisions based on connections to vertices
    that might be later lost in arbitration.

    - Phase 6:  Aggregate remaining vertices and avoid small aggregates (e.g.,
    singletons) at all costs. Typically, most everything should
    be aggregated by Phase's 1-5.  One way that we could still have
    unaggregated vertices is if processor p was never assigned a
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


    One final note about the use of ArbitrateAndCommunicate(). No
    arbitration occurs (which means the procWinner[] is not set as well) for a
    global shared id if and only if all weights on all processors corresponding
    to this id is zero. Thus, the general idea is that any global id where we
    want arbitration to occur should have at least one associated weight on
    one processor which is nonzero. Any global id where we don't want any
    arbitration should have all weights set to 0.

    Note: procWinners is also set to MyPid() by ArbitrateAndCommunicate()
    for any nonshared gid's with a nonzero weight.
    */
    void AggregateLeftovers(GraphBase const &graph, Aggregates &aggregates) const; //AggregateLeftovers

      //@}

      //! @name Utilities
      //@{

      /*! @brief Build a list of candidate root nodes.

      Candidates are vertices not adjacent to already aggregated vertices.
      */
    void RootCandidates(my_size_t nVertices, ArrayView<const LO> & vertex2AggId, GraphBase const &graph,
                        ArrayRCP<LO> &candidates, my_size_t &nCandidates, global_size_t &nCandidatesGlobal) const; //RootCandidates

      //! @brief Attempt to clean up aggregates that are too small.
    int RemoveSmallAggs(Aggregates& aggregates, int min_size,
                        RCP<Xpetra::Vector<double,LO,GO,NO> > & distWeights, const MueLu::CoupledAggregationCommHelper<LO,GO,NO,LMO> & myWidget) const; //RemoveSmallAggs

    //@}

  private:
    double phase3AggCreation_;
    int    minNodesPerAggregate_;

    // JG TODO: rename variables:
    //  Adjacent-> adjacent
    //  homogenization of variables names :
    //  - colj and j
    //  - i and iNode
    //  - k->kNode
    //  - ...

  }; //class LeftoverAggregationAlgorithm

} //namespace MueLu

//      graphName_("UC_CleanUp")

#define MUELU_LEFTOVERAGGREGATIONALGORITHM_SHORT
#endif // MUELU_LEFTOVERAGGREGATIONALGORITHM_DECL_HPP
