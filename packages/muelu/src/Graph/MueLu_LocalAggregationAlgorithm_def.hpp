#ifndef MUELU_LOCALAGGREGATIONFACTORY_HPP_DEF
#define MUELU_LOCALAGGREGATIONFACTORY_HPP_DEF

#include "MueLu_LocalAggregationAlgorithm_decl.hpp"

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

namespace MueLu {

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


  /*!
    @class LocalAggregationAlgorithm class.
    @brief Algorithm for coarsening a graph with uncoupled aggregation.

    This method has two phases.  The first is a local clustering algorithm.  The second creates aggregates
    that can include unknowns from more than one process.

  */


  //! Constructor.
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::LocalAggregationAlgorithm(RCP<FactoryBase> const &graphFact)
    : ordering_(NATURAL), minNodesPerAggregate_(1), maxNeighAlreadySelected_(0)
  { }

  //! Destructor.
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~LocalAggregationAlgorithm() { }

  //@}

  //! @name Set/get methods.
  //@{

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetOrdering(Ordering ordering)                          { ordering_                = ordering;                }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetMinNodesPerAggregate(int minNodesPerAggregate)       { minNodesPerAggregate_    = minNodesPerAggregate;    }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) { maxNeighAlreadySelected_ = maxNeighAlreadySelected; }
    
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Ordering LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetOrdering()                const { return ordering_;                }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  int      LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMinNodesPerAggregate()    const { return minNodesPerAggregate_;    }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  int      LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMaxNeighAlreadySelected() const { return maxNeighAlreadySelected_; }

  //@}

  //! @name Aggregation methods.
  //@{

  /*! @brief Local aggregation. */
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoarsenUncoupled(Graph const & graph, Aggregates & aggregates) const
  {
    Monitor m(*this, "Coarsen Uncoupled");

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
                TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoarsenUncoupled(): Error: couldn't allocate memory for supernode! length=" + Teuchos::toString(length));
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

  //   private:
  //     //! Aggregation options (TODO: Teuchos::ParameterList?)
  //     Ordering ordering_;                /**<  natural, random, graph           */
  //     int      minNodesPerAggregate_;    /**<  aggregate size control           */
  //     int      maxNeighAlreadySelected_; /**<  complexity control               */

  //! @name Utilities
  //@{

  /*! @brief Utility to take a list of integers and reorder them randomly (by using a local permutation).
    @param list On input, a bunch of integers. On output, the same integers in a different order
    that is determined randomly.
  */
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomReorder(Teuchos::ArrayRCP<LO> list) const
  { //TODO: replace int
    int n = list.size();
    for(int i=0; i<n-1; i++) {
      std::swap(list[i], list[RandomOrdinal(i,n-1)]);
    }
  } 

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  /*! @brief Generate a random number in the range [min, max] */
  int LocalAggregationAlgorithm<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::RandomOrdinal(int min, int max) const
  {
    return min + static_cast<int>((max-min+1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
  }

  //@}
  
} //namespace MueLu

#endif //ifndef MUELU_LOCALAGGREGATIONFACTORY_HPP_DEF
