#ifndef MUELU_UCAGGREGATIONFACTORY_HPP
#define MUELU_UCAGGREGATIONFACTORY_HPP

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_LocalAggregationFactory.hpp"
#include "MueLu_UCAggregationCommHelper.hpp"

namespace MueLu {

  /*!
    @class UCAggregationFactory class.
    @brief Factory for coarsening a graph with uncoupled aggregation.

    This method has two phases.  The first is a local clustering algorithm.  The second creates aggregates
    that can include unknowns from more than one process.

    - TODO remove template dependence on Scalar
  */

  template <class Scalar        = double,
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
  class UCAggregationFactory : public Teuchos::Describable {
    //#include "MueLu_UseShortNamesOrdinal.hpp"
#include "MueLu_UseShortNames.hpp"

    typedef GO global_size_t; //TODO
    typedef LO my_size_t;     //TODO

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    UCAggregationFactory() :
      localAggregationFactory_(),
      phase3AggCreation_(.5), 
      graphName_("unnamed") //, Algorithm_("notSpecified")
    {

    }

    //! Constructor.
    UCAggregationFactory(RCP<CoalesceDropFactory> const & cdFact) :
      localAggregationFactory_(cdFact),
      phase3AggCreation_(.5), 
      graphName_("unnamed") //, Algorithm_("notSpecified")
    {

    }

    //! Destructor.
    virtual ~UCAggregationFactory() {}
    //@}

    //! @name Set/get methods.
    //@{

    // Options common to UCAggregationFactory and LocalAggregationFactory are hold by LocalAggregationFactory
    // It's not perfect but it's the way it is. MinNodesPerAggregate is the only option used on both factory.
    void SetPrintFlag(int printFlag)                             { localAggregationFactory_.SetPrintFlag(printFlag);                             } //TODO: to be removed
    void SetOrdering(Ordering ordering)                          { localAggregationFactory_.SetOrdering(ordering);                               }
    void SetMinNodesPerAggregate(int minNodesPerAggregate)       { localAggregationFactory_.SetMinNodesPerAggregate(minNodesPerAggregate);       }
    void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) { localAggregationFactory_.SetMaxNeighAlreadySelected(maxNeighAlreadySelected); }
    void SetGraphName(std::string const &graphName)              { localAggregationFactory_.SetGraphName(graphName);                             }
    
    double      GetPrintFlag()               const { return localAggregationFactory_.GetPrintFlag();               } //TODO: to be removed
    Ordering    GetOrdering()                const { return localAggregationFactory_.GetOrdering();                }
    int         GetMinNodesPerAggregate()    const { return localAggregationFactory_.GetMinNodesPerAggregate();    }
    int         GetMaxNeighAlreadySelected() const { return localAggregationFactory_.GetmaxNeighAlreadySelected(); }
    std::string GetGraphName()               const { return localAggregationFactory_.GetGraphName();               }

    // Option specific to UCAggregationFactory
    void   SetPhase3AggCreation(double phase3AggCreation) { phase3AggCreation_ = phase3AggCreation; }
    double GetPhase3AggCreation() const { return phase3AggCreation_; }

    //@}

    //! @name Build methods.
    //@{

    /*! @brief Build aggregates.

    - TODO reuse of aggregates
    - TODO check if called twice (bug TEUCHOS_TEST_EQUALITY)
    */
    void Build(Level &currentLevel) const
    {
      //TODO check for reuse of aggregates here
      //FIXME should there be some way to specify the name of the graph in the needs table, i.e., could
      //FIXME there ever be more than one graph?
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("UCAggregationFactory::Build"));
      timer->start(true);

      currentLevel.Request("Graph");
      currentLevel.Request("Aggregates");

      // Graph creation is delegated to localAggregationFactory_
      // Maybe it is not the best way to do it.
      localAggregationFactory_.Build(currentLevel);

      RCP<Graph> graph;
      currentLevel.Get("Graph",graph);
      currentLevel.Release("Graph");
      
      RCP<Aggregates> aggregates;
      currentLevel.Get("Aggregates",aggregates);
      currentLevel.Release("Aggregates");
     
      std::string name =  "UC_CleanUp";
      AggregateLeftOvers(*aggregates, name, *graph);

      currentLevel.Set("Aggregates",aggregates);

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(graph->GetComm()));
    }

    RCP<Aggregates> Build(const Graph& graph) const
    {
      RCP<Aggregates> aggregates = localAggregationFactory_.Build(graph);
      std::string name =  "UC_CleanUp";
      AggregateLeftOvers(*aggregates, name, graph);
      return aggregates;
    }

  private:

    //! 
    LocalAggregationFactory localAggregationFactory_;

    //! Aggregation options
    double phase3AggCreation_;    /**<  Steers how the MIS  and Uncoupled 
                                     handle phase 3 of aggregation.     
                                     Values near 0 create few additional
                                     aggregates.Large values create many
                                     additional aggregates. Convergence 
                                     can be improve by new  
                                     aggregates but nonzero fill-in     
                                     increases on coarse meshes.        
                                     Default: .5                         */

    // unused 
    // aggregation algorithm type
    // std::string Algorithm_;

    //! user-defined graph label
    std::string graphName_;

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
           
    Phase 1b: Invoke ArbitrateAndCommunicate() to ensure that
    all processors have the same view of aggregated vertices
    (e.g., to which aggregate they have been assigned and
    which processor owns that aggregate).
    Phase 2:  Check for vertices (local or nonlocal) which are Adjacent
    to root nodes. Tentatively assign these to the aggregate
    associated with the root. Arbitrate any cases where 
    several processors claim the same vertex for one of 
    their aggregates via ArbitrateAndCommunicate().
    Phase 3:  Try to create new aggregates if it looks like there are
    root node candidates which have many unaggregated neighbors.
    This decision to make a new aggregate is based on only local
    information. However, the new aggregate will be tentatively
    assigned any unaggregated ghost vertices. Arbitration is again
    done by ArbitrateAndCommunicate() where local vertices
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
    ArbitrateAndCommunicate() is again used for arbitration
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
    ArbitrateAndCommunicate() i
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
    int AggregateLeftOvers(Aggregates &aggregates,  std::string const & label, Graph const &graph) const
    {
      my_size_t nVertices = graph.GetNodeNumVertices();
      int exp_nRows    = aggregates.GetMap()->getNodeNumElements(); // Tentative fix... was previously exp_nRows = nVertices + graph.GetNodeNumGhost();
      int myPid        = graph.GetComm()->getRank();
      my_size_t nAggregates  = aggregates.GetNumAggregates();

      int minNodesPerAggregate = GetMinNodesPerAggregate();

      const RCP<const Map> nonUniqueMap = aggregates.GetMap();
      const RCP<const Map> uniqueMap = graph.GetDomainMap(); // Q: DomainMap or RowMap??

      MueLu::UCAggregationCommHelper<double,LO,GO,NO,LMO> myWidget(uniqueMap, nonUniqueMap);

      RCP<Xpetra::Vector<double,LO,GO,NO> > distWeights = Xpetra::VectorFactory<double,LO,GO,NO>::Build(nonUniqueMap);

      // Aggregated vertices not "definitively" assigned to processors are
      // arbitrated by ArbitrateAndCommunicate(). There is some
      // additional logic to prevent losing root nodes in arbitration.
      {
        ArrayRCP<const LO> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
        ArrayRCP<const LO> procWinner   = aggregates.GetProcWinner()->getData(0);
        ArrayRCP<double>    weights     = distWeights->getDataNonConst(0);
          
        distWeights->putScalar(0.);
        for (size_t i=0;i<nonUniqueMap->getNodeNumElements();i++) {
          if (procWinner[i] == MUELU_UNASSIGNED) {
            if (vertex2AggId[i] != MUELU_UNAGGREGATED) {
              weights[i] = 1.;
              if (aggregates.IsRoot(i)) weights[i] = 2.;
            }
          }
        }
          
        // views on distributed vectors are freed here.
      }

      myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
      distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive

      // Tentatively assign any vertex (ghost or local) which neighbors a root
      // to the aggregate associated with the root.
      {
        ArrayRCP<LO>       vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
        ArrayRCP<const LO> procWinner   = aggregates.GetProcWinner()->getData(0);
        ArrayRCP<double>   weights      = distWeights->getDataNonConst(0);

        for (my_size_t i = 0; i < nVertices; i++) { 
          if ( aggregates.IsRoot(i) && (procWinner[i] == myPid) ) {
              
            // neighOfINode is the neighbor node list of node 'i'.
            ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);
              
            for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int colj = *it;
              if (vertex2AggId[colj] == MUELU_UNAGGREGATED) {
                weights[colj]= 1.;
                vertex2AggId[colj] = vertex2AggId[i];
              }
            }
          }
        }

        // views on distributed vectors are freed here.
      }

      myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
      distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive

      // Record the number of aggregated vertices
      GO total_phase_one_aggregated = 0;
      {
        ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
            
        GO phase_one_aggregated = 0;
        for (my_size_t i = 0; i < nVertices; i++) {
          if (vertex2AggId[i] != MUELU_UNAGGREGATED)
            phase_one_aggregated++;
        }
          
        sumAll(graph.GetComm(), phase_one_aggregated, total_phase_one_aggregated);

        GO local_nVertices = nVertices, total_nVertices = 0;
        sumAll(graph.GetComm(), local_nVertices, total_nVertices);
          
        /* Among unaggregated points, see if we can make a reasonable size    */
        /* aggregate out of it. We do this by looking at neighbors and seeing */
        /* how many are unaggregated and on my processor. Loosely,            */
        /* base the number of new aggregates created on the percentage of     */
        /* unaggregated nodes.                                                */

        ArrayRCP<double>    weights      = distWeights->getDataNonConst(0);

        double factor = 1.;
        factor = ((double) total_phase_one_aggregated)/((double)(total_nVertices + 1));
        factor = pow(factor, GetPhase3AggCreation());
          
        for (my_size_t i = 0; i < nVertices; i++) {
          if (vertex2AggId[i] == MUELU_UNAGGREGATED) 
            {
                
              // neighOfINode is the neighbor node list of node 'iNode'.
              ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);
              int rowi_N = neighOfINode.size();
                
              int nonaggd_neighbors = 0;
              for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
                int colj = *it;
                if (vertex2AggId[colj] == MUELU_UNAGGREGATED && colj < nVertices)
                  nonaggd_neighbors++;
              }
              if (  (nonaggd_neighbors > minNodesPerAggregate) &&
                    (((double) nonaggd_neighbors)/((double) rowi_N) > factor))
                {
                  vertex2AggId[i] = (nAggregates)++;
                  for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
                    int colj = *it;
                    if (vertex2AggId[colj]==MUELU_UNAGGREGATED) {
                      vertex2AggId[colj] = vertex2AggId[i];
                      if (colj < nVertices) weights[colj] = 2.;
                      else                  weights[colj] = 1.;
                    }
                  }
                  aggregates.SetIsRoot(i);
                  weights[i] = 2.;
                }
            }
        } // for (i = 0; i < nVertices; i++)

          // views on distributed vectors are freed here.
      }

      myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
      distWeights->putScalar(0.);//All tentatively assigned vertices are now definitive

      if ( GetPrintFlag() < 7) { //FIXME
        GO Nphase1_agg = nAggregates;
        GO total_aggs;

        sumAll(graph.GetComm(), Nphase1_agg, total_aggs);

        if (myPid == 0) {
          std::cout << "Aggregation(" << label << ") : Phase 1 - nodes aggregated = " << total_phase_one_aggregated << std::endl;
          std::cout << "Aggregation(" << label << ") : Phase 1 - total aggregates = " << total_aggs << std::endl;
        }
        GO i = nAggregates - Nphase1_agg;

        { GO ii; sumAll(graph.GetComm(),i,ii); i = ii; }
        if ( myPid == 0 ) {
          std::cout << "Aggregation(" << label << ") : Phase 3 - additional aggregates = " << i << std::endl;
        }
      }

      // Determine vertices that are not shared by setting Temp to all ones
      // and doing NonUnique2NonUnique(..., ADD). This sums values of all
      // local copies associated with each Gid. Thus, sums > 1 are shared.

      //         std::cout << "exp_nrows=" << exp_nRows << " (nVertices= " << nVertices << ", numGhost=" << graph.GetNodeNumGhost() << ")" << std::endl;
      //         std::cout << "nonUniqueMap=" << nonUniqueMap->getNodeNumElements() << std::endl;

      RCP<Xpetra::Vector<double,LO,GO,NO> > temp_ = Xpetra::VectorFactory<double,LO,GO,NO> ::Build(nonUniqueMap);
      temp_->putScalar(1.);  

      RCP<Xpetra::Vector<double,LO,GO,NO> > tempOutput_ = Xpetra::VectorFactory<double,LO,GO,NO> ::Build(nonUniqueMap);
      tempOutput_->putScalar(0.); 

      myWidget.NonUnique2NonUnique(*temp_, *tempOutput_, Xpetra::ADD);

      std::vector<bool> gidNotShared(exp_nRows);
      {
        ArrayRCP<const double> tempOutput = tempOutput_->getData(0);
        for (int i = 0; i < exp_nRows; i++) {
          if (tempOutput[i] > 1.) 
            gidNotShared[i] = false; 
          else 
            gidNotShared[i] = true; 
        }
      }

      // Phase 4. 
      double nAggregatesTarget; 
      nAggregatesTarget = ((double)  uniqueMap->getGlobalNumElements())* (((double) uniqueMap->getGlobalNumElements())/ ((double) graph.GetGlobalNumEdges())); 
        
      GO nAggregatesLocal=nAggregates, nAggregatesGlobal; sumAll(graph.GetComm(), nAggregatesLocal, nAggregatesGlobal); 

      LO minNAggs; minAll(graph.GetComm(), nAggregates, minNAggs); 
      LO maxNAggs; maxAll(graph.GetComm(), nAggregates, maxNAggs);

      //
      // Only do this phase if things look really bad. THIS
      // CODE IS PRETTY EXPERIMENTAL 
      //
#define MUELU_PHASE4BUCKETS 6
      if ((nAggregatesGlobal < graph.GetComm()->getSize()) &&
          (2.5*nAggregatesGlobal < nAggregatesTarget) &&
          (minNAggs ==0) && (maxNAggs <= 1)) {

        // Modify seed of the random algorithm used by temp_->randomize()
        {
          typedef Teuchos::ScalarTraits<double> scalarTrait; // temp_ is of type double.
          scalarTrait::seedrandom(static_cast<unsigned int>(myPid*2 + (int) (11*scalarTrait::random())));
          int k = (int)ceil( (10.*myPid)/graph.GetComm()->getSize());
          for (int i = 0; i < k+7; i++) scalarTrait::random();
          temp_->setSeed(static_cast<unsigned int>(scalarTrait::random()));
        }

        temp_->randomize(); 

        ArrayRCP<double> temp = temp_->getDataNonConst(0);

        // build a list of candidate root nodes (vertices not adjacent
        // to aggregated vertices)

        my_size_t nCandidates = 0;
        global_size_t nCandidatesGlobal;

        ArrayRCP<LO> candidates = Teuchos::arcp<LO>(nVertices+1);

        double priorThreshold = 0.;
        for (int kkk = 0 ; kkk < MUELU_PHASE4BUCKETS; kkk++) {

          {
            ArrayRCP<const LO> vertex2AggId = aggregates.GetVertex2AggId()->getData(0);
            ArrayView<const LO> vertex2AggIdView = vertex2AggId();
            RootCandidates(nVertices, vertex2AggIdView, graph, candidates, nCandidates, nCandidatesGlobal);
            // views on distributed vectors are freed here.
          }

          double nTargetNewGuys =  nAggregatesTarget - nAggregatesGlobal;
          double threshold      =  priorThreshold + (1. - priorThreshold)*nTargetNewGuys/(nCandidatesGlobal + .001);
         
          threshold = (threshold*(kkk+1.))/((double) MUELU_PHASE4BUCKETS);
          priorThreshold = threshold;

          {
            ArrayRCP<LO>     vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
            ArrayRCP<double> weights      = distWeights->getDataNonConst(0);

            for (int k = 0; k < nCandidates ; k++ ) { 
              int i = candidates[k];                  
              if ((vertex2AggId[i] == MUELU_UNAGGREGATED) && (fabs(temp[i])  < threshold)) {
                // Note: priorThreshold <= fabs(temp[i]) <= 1

                // neighOfINode is the neighbor node list of node 'iNode'.
                ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

                if (neighOfINode.size() > minNodesPerAggregate) { //TODO: check if this test is exactly was we want to do
                  int count = 0;
                  for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
                    int Adjacent    = *it;
                    // This might not be true if someone close to i
                    // is chosen as a root via fabs(temp[]) < Threshold
                    if (vertex2AggId[Adjacent] == MUELU_UNAGGREGATED){
                      count++;
                      vertex2AggId[Adjacent] = nAggregates;
                      weights[Adjacent] = 1.;
                    }
                  }
                  if (count >= minNodesPerAggregate) {
                    vertex2AggId[i] = nAggregates++;
                    weights[i] = 2.;
                    aggregates.SetIsRoot(i);
                  }
                  else { // undo things
                    for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
                      int Adjacent    = *it;
                      if (vertex2AggId[Adjacent] == nAggregates){
                        vertex2AggId[Adjacent] = MUELU_UNAGGREGATED;
                        weights[Adjacent] = 0.;
                      }
                    }
                  }
                }
              }
            }
            // views on distributed vectors are freed here.
          }
          myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
          distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive
          nAggregatesLocal=nAggregates;
          sumAll(graph.GetComm(), nAggregatesLocal, nAggregatesGlobal);

          // check that there are no aggregates sizes below minNodesPerAggregate
           
          aggregates.SetNumAggregates(nAggregates);

          RemoveSmallAggs(aggregates, minNodesPerAggregate, distWeights, myWidget);
            
          nAggregates = aggregates.GetNumAggregates();
        }   // one possibility
      }

      // Initialize things for Phase 5. This includes building the transpose
      // of the matrix ONLY for transposed rows that correspond to unaggregted
      // ghost vertices. Further, the transpose is only a local transpose. 
      // Nonzero edges which exist on other processors are not represented.


      int observedNAgg=-1; //number of aggregates that contain vertices on this process

      {
        ArrayRCP<LO>       vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
        ArrayRCP<const LO> procWinner   = aggregates.GetProcWinner()->getData(0);
        for(LO k = 0; k < vertex2AggId.size(); ++k )
          if(vertex2AggId[k]>observedNAgg)
            observedNAgg=vertex2AggId[k];
        observedNAgg++; 
      }

      ArrayRCP<int> Mark = Teuchos::arcp<int>(exp_nRows+1);
      ArrayRCP<int> agg_incremented = Teuchos::arcp<int>(observedNAgg);
      ArrayRCP<int> SumOfMarks = Teuchos::arcp<int>(observedNAgg);

      for (int i = 0; i < exp_nRows; i++)   Mark[i] = MUELU_DISTONE_VERTEX_WEIGHT;
      for (int i = 0; i < agg_incremented.size(); i++) agg_incremented[i] = 0;
      for (int i = 0; i < SumOfMarks.size(); i++) SumOfMarks[i] = 0;

      // Grab the transpose matrix graph for unaggregated ghost vertices.
      //     a) count the number of nonzeros per row in the transpose
      std::vector<int> RowPtr(exp_nRows+1-nVertices);
      //{
      ArrayRCP<const LO> vertex2AggIdCst = aggregates.GetVertex2AggId()->getData(0);

      for (int i = nVertices; i < exp_nRows ;  i++) RowPtr[i-nVertices] = 0;
      for (int i = 0; i < nVertices;  i++) {

        // neighOfINode is the neighbor node list of node 'iNode'.
        ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

        for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
          int j = *it;
          if ( (j >= nVertices) && (vertex2AggIdCst[j] == MUELU_UNAGGREGATED)){
            RowPtr[j-nVertices]++;
          }
        }
      }

      //     b) Convert RowPtr[i] to point to 1st first nnz spot in row i.

      int iSum = 0, iTemp;
      for (int i = nVertices; i < exp_nRows ;  i++) {
        iTemp = RowPtr[i-nVertices];
        RowPtr[i-nVertices] = iSum;
        iSum += iTemp;
      }
      RowPtr[exp_nRows-nVertices] = iSum;
      std::vector<LO> cols(iSum+1);
         
      //     c) Traverse matrix and insert entries in proper location.
      for (int i = 0; i < nVertices;  i++) {

        // neighOfINode is the neighbor node list of node 'iNode'.
        ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);
          
        for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
          int j = *it;
          if ( (j >= nVertices) && (vertex2AggIdCst[j] == MUELU_UNAGGREGATED)){
            cols[RowPtr[j-nVertices]++] = i;
          }
        }
      }

      //     d) RowPtr[i] points to beginning of row i+1 so shift by one location.
      for (int i = exp_nRows; i > nVertices;  i--)
        RowPtr[i-nVertices] = RowPtr[i-1-nVertices];
      RowPtr[0] = 0;
          
      // views on distributed vectors are freed here.
      vertex2AggIdCst = Teuchos::null;
      //}

      int bestScoreCutoff;
      int thresholds[10] = {300,200,100,50,25,13,7,4,2,0};

      // Stick unaggregated vertices into existing aggregates as described above. 
      for (int kk = 0; kk < 10; kk += 2) {
        bestScoreCutoff = thresholds[kk];
        for (int i = 0; i < exp_nRows; i++) {

          ArrayRCP<LO> vertex2AggId     = aggregates.GetVertex2AggId()->getDataNonConst(0);
          ArrayRCP<const LO> procWinner = aggregates.GetProcWinner()->getData(0);
          ArrayRCP<double> weights       = distWeights->getDataNonConst(0);

          if (vertex2AggId[i] == MUELU_UNAGGREGATED) {

            // neighOfINode is the neighbor node list of node 'iNode'.
            ArrayView<const LO> neighOfINode;

            // Grab neighboring vertices which is either in graph for local ids
            // or sits in transposed fragment just constructed above for ghosts.
            if (i < nVertices) {
              neighOfINode = graph.getNeighborVertices(i);
            }
            else {
              LO *rowi_col = NULL, rowi_N;
              rowi_col = &(cols[RowPtr[i-nVertices]]);
              rowi_N   = RowPtr[i+1-nVertices] - RowPtr[i-nVertices];
                
              neighOfINode = ArrayView<const LO>(rowi_col, rowi_N);
            }
            for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int Adjacent    = *it;
              int AdjacentAgg = vertex2AggId[Adjacent];

              //Adjacent is aggregated and either I own the aggregate
              // or I could own the aggregate after arbitration.
              if ((AdjacentAgg != MUELU_UNAGGREGATED) && 
                  ((procWinner[Adjacent] == myPid) ||     
                   (procWinner[Adjacent] == MUELU_UNASSIGNED))){
                SumOfMarks[AdjacentAgg] += Mark[Adjacent];
              }
            }
            int best_score = MUELU_NOSCORE;
            int best_agg = -1;
            int BestMark = -1;
            bool cannotLoseAllFriends=false; // Used to address possible loss of vertices in arbitration of shared nodes discussed above. (Initialized to false only to avoid a compiler warning).

            for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int Adjacent    = *it;
              int AdjacentAgg = vertex2AggId[Adjacent];
              //Adjacent is unaggregated, has some value and no
              //other processor has definitively claimed him
              if ((AdjacentAgg != MUELU_UNAGGREGATED) && 
                  (SumOfMarks[AdjacentAgg] != 0) &&
                  ((procWinner[Adjacent] == myPid) ||
                   (procWinner[Adjacent] == MUELU_UNASSIGNED ))) {

                // first figure out the penalty associated with
                // AdjacentAgg having already been incremented 
                // during this phase, then compute score.

                double penalty = (double) (INCR_SCALING*agg_incremented[AdjacentAgg]);
                if (penalty > MUELU_PENALTYFACTOR*((double)SumOfMarks[AdjacentAgg]))
                  penalty = MUELU_PENALTYFACTOR*((double)SumOfMarks[AdjacentAgg]);
                int score = SumOfMarks[AdjacentAgg]- ((int) floor(penalty));

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
            for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int Adjacent    = *it;
              int AdjacentAgg = vertex2AggId[Adjacent];
              if (AdjacentAgg >= 0) SumOfMarks[AdjacentAgg] = 0;
            }
            // Tentatively assign vertex to best_agg. 
            if ( (best_score >= bestScoreCutoff) && (cannotLoseAllFriends)) { 

              TEST_FOR_EXCEPTION(best_agg == -1 || BestMark == -1, MueLu::Exceptions::RuntimeError, "MueLu::UCAggregationFactory internal error"); // should never happen

              vertex2AggId[i] = best_agg;
              weights[i] = best_score;
              agg_incremented[best_agg]++;
              Mark[i] = (int) ceil(   ((double) BestMark)/2.);
            }
          }

          // views on distributed vectors are freed here.
        }

        myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
        distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive
      }

      // Phase 6: Aggregate remain unaggregated vertices and try at all costs
      //          to avoid small aggregates.
      //          One case where we can find ourselves in this situation
      //          is if all vertices vk adjacent to v have already been
      //          put in other processor's aggregates and v does not have
      //          a direct connection to a local vertex in any of these
      //          aggregates.

      int Nleftover = 0, Nsingle = 0;
      {

        ArrayRCP<LO> vertex2AggId     = aggregates.GetVertex2AggId()->getDataNonConst(0);
        ArrayRCP<double> weights       = distWeights->getDataNonConst(0);
        ArrayRCP<const LO> procWinner = aggregates.GetProcWinner()->getData(0);

        int count = 0;
        for (my_size_t i = 0; i < nVertices; i++) { 
          if ((vertex2AggId[i] == MUELU_UNAGGREGATED) ) {
            Nleftover++;

            // neighOfINode is the neighbor node list of node 'iNode'.
            ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

            // We don't want too small of an aggregate. So lets see if there is an
            // unaggregated neighbor that we can also put with this vertex

            vertex2AggId[i] = nAggregates;
            weights[i] = 1.;
            if (count == 0) aggregates.SetIsRoot(i);
            count++;
            for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
              int j = *it;
              if ((j != i)&&(vertex2AggId[j] == MUELU_UNAGGREGATED)&&
                  (j < nVertices)) {
                vertex2AggId[j] = nAggregates;
                weights[j] = 1.;
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
            for (my_size_t i = 0; i < nVertices; i++) {
              if ((vertex2AggId[i] == nAggregates) && (procWinner[i] == myPid)){
                vertex2AggId[i] = 0;
                aggregates.SetIsRoot(i,false);
              }
            }
          }
          else {
            Nsingle++;
            nAggregates++;
          }
        }

        // views on distributed vectors are freed here.
      }

      myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, false);

      if (GetPrintFlag() < 7) { //FIXME
        GO total_Nsingle=0;   sumAll(graph.GetComm(), (GO)Nsingle,     total_Nsingle);    
        GO total_Nleftover=0; sumAll(graph.GetComm(), (GO)Nleftover,   total_Nleftover);
        GO total_aggs;        sumAll(graph.GetComm(), (GO)nAggregates, total_aggs);
        if ( myPid == 0 ) {
          std::cout << "Aggregation(" << label << ") : Phase 3 - total aggregates = " << total_aggs << std::endl;
          std::cout << "Aggregation(" << label << ") : Phase 6 - leftovers = " <<  total_Nleftover << " and singletons = " << total_Nsingle << std::endl;
        }
      }

      aggregates.SetNumAggregates(nAggregates);

      return 0;
    } //AggregateLeftOvers

      //@}

      //! @name Utilities
      //@{

      /*! @brief Build a list of candidate root nodes.

      Candidates are vertices not adjacent to already aggregated vertices.
      */
    void RootCandidates(my_size_t nVertices, ArrayView<const LO> & vertex2AggId, Graph const &graph,
                        ArrayRCP<LO> &candidates, my_size_t &nCandidates, global_size_t &nCandidatesGlobal) const
    {
      nCandidates = 0;
       
      for (my_size_t i = 0; i < nVertices; i++ ) {
        if (vertex2AggId[i] == MUELU_UNAGGREGATED) {
          bool noAggdNeighbors = true;

          // neighOfINode is the neighbor node list of node 'iNode'.
          ArrayView<const LO> neighOfINode = graph.getNeighborVertices(i);

          for (typename ArrayView<const LO>::const_iterator it = neighOfINode.begin(); it != neighOfINode.end(); ++it) {
            int adjacent    = *it;
            if (vertex2AggId[adjacent] != MUELU_UNAGGREGATED) 
              noAggdNeighbors = false;
          }
          if (noAggdNeighbors == true) candidates[nCandidates++] = i;
        }
      }

      sumAll(graph.GetComm(), (GO)nCandidates, nCandidatesGlobal);

    } //RootCandidates

      //! @brief Attempt to clean up aggregates that are too small.
    int RemoveSmallAggs(Aggregates& aggregates, int min_size,
                        RCP<Xpetra::Vector<double,LO,GO,NO> > & distWeights, const MueLu::UCAggregationCommHelper<double,LO,GO,NO,LMO> & myWidget) const
    {
      int myPid = aggregates.GetMap()->getComm()->getRank();
        
      LO nAggregates = aggregates.GetNumAggregates();

      ArrayRCP<LO> procWinner   = aggregates.GetProcWinner()->getDataNonConst(0);
      ArrayRCP<LO> vertex2AggId = aggregates.GetVertex2AggId()->getDataNonConst(0);
      LO size = procWinner.size();

      //ArrayRCP<int> AggInfo = Teuchos::arcp<int>(nAggregates+1);
      //aggregates.ComputeAggSizes(AggInfo);
      ArrayRCP<LO> AggInfo = aggregates.ComputeAggregateSizes();

      ArrayRCP<double> weights = distWeights->getDataNonConst(0);

      // Make a list of all aggregates indicating New AggId
      // Use AggInfo array for this.

      LO NewNAggs = 0; 
      for (LO i = 0; i < nAggregates; i++) {
        if ( AggInfo[i] < min_size) { 
          AggInfo[i] =  MUELU_UNAGGREGATED;
        }
        else AggInfo[i] = NewNAggs++;
      }

      for (LO k = 0; k < size; k++ ) {
        if (procWinner[k] == myPid) {
          if (vertex2AggId[k] !=  MUELU_UNAGGREGATED) {
            vertex2AggId[k] = AggInfo[vertex2AggId[k]];
            weights[k] = 1.;
          }
          if (vertex2AggId[k] ==  MUELU_UNAGGREGATED) 
            aggregates.SetIsRoot(k,false);
        }
      }
      nAggregates = NewNAggs;

      myWidget.ArbitrateAndCommunicate(*distWeights, aggregates, true);
      distWeights->putScalar(0.); // All tentatively assigned vertices are now definitive

      // procWinner is not set correctly for aggregates which have 
      // been eliminated
      for (LO i = 0; i < size; i++) {
        if (vertex2AggId[i] == MUELU_UNAGGREGATED) 
          procWinner[i] = MUELU_UNASSIGNED;
      }
      aggregates.SetNumAggregates(nAggregates);

      return 0; //TODO
    } //RemoveSmallAggs

  
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
