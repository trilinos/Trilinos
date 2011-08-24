#ifndef MUELU_LOCALAGGREGATIONFACTORY_HPP
#define MUELU_LOCALAGGREGATIONFACTORY_HPP

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_LinkedList.hpp"

#include "MueLu_Memory.hpp"

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

namespace MueLu {

  namespace AggOptions {
    /* Options defining how to pick-up the next root node in the local aggregation procedure */
    enum Ordering {
      NATURAL = 0, /* node ordering   */
      RANDOM  = 1, /* random ordering */
      GRAPH   = 2  /* graph ordering  */
    };
  } // namespace AggOptions

  using namespace AggOptions;

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

  /* In the algorithm, aggStat[]=READY/NOTSEL/SELECTED indicates whether a node has been aggregated. */
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

  /*!
    @class LocalAggregationFactory class.
    @brief Factory for coarsening a graph with uncoupled aggregation.

    This method has two phases.  The first is a local clustering algorithm.  The second creates aggregates
    that can include unknowns from more than one process.

  */

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class LocalAggregationFactory : public SingleLevelFactoryBase {
#include "MueLu_UseShortNamesOrdinal.hpp"

    typedef GO global_size_t; //TODO
    typedef LO my_size_t; //TODO

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    LocalAggregationFactory(RCP<SingleLevelFactoryBase> const &coalesceDropFact=Teuchos::null) :
      printFlag_(0), //TODO: to be removed
      ordering_(NATURAL), minNodesPerAggregate_(1), maxNeighAlreadySelected_(0),
      coalesceDropFact_(coalesceDropFact), 
      graphName_("unnamed") //, Algorithm_("notSpecified"), 
    {

    }

  //! Destructor.
  virtual ~LocalAggregationFactory() {}
  //@}

  //! @name Set/get methods.
  //@{
  void SetPrintFlag(int printFlag)                             { printFlag_               = printFlag;               } //TODO: to be removed
  void SetOrdering(Ordering ordering)                          { ordering_                = ordering;                }
  void SetMinNodesPerAggregate(int minNodesPerAggregate)       { minNodesPerAggregate_    = minNodesPerAggregate;    }
  void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) { maxNeighAlreadySelected_ = maxNeighAlreadySelected; }
  void SetGraphName(std::string const &graphName)              { graphName_               = graphName;               }
    
  double      GetPrintFlag()               const { return printFlag_;               } //TODO: to be removed
  Ordering    GetOrdering()                const { return ordering_;                }
  int         GetMinNodesPerAggregate()    const { return minNodesPerAggregate_;    }
  int         GetMaxNeighAlreadySelected() const { return maxNeighAlreadySelected_; }
  std::string GetGraphName()               const { return graphName_;               }
  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates.

  - TODO reuse of aggregates
  - TODO check if called twice (bug TEUCHOS_TEST_EQUALITY)
  */
  bool Build(Level &currentLevel) const
  {
    //TODO check for reuse of aggregates here
    //FIXME should there be some way to specify the name of the graph in the needs table, i.e., could
    //FIXME there ever be more than one graph?
    RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("LocalAggregationFactory::Build"));
    timer->start(true);

    currentLevel.Request("Graph");
    if (coalesceDropFact_ != Teuchos::null)
      coalesceDropFact_->SingleLevelBuild(currentLevel);
    else
      currentLevel.GetDefaultFactory("Graph")->SingleLevelBuild(currentLevel);

    RCP<Graph> graph;
    currentLevel.Get("Graph",graph);
    currentLevel.Release("Graph");
    RCP<Aggregates> aggregates = Build(*graph);
    currentLevel.Set("Aggregates",aggregates);

    timer->stop();
    MemUtils::ReportTimeAndMemory(*timer, *(graph->GetComm()));

    return true; //??
  }

  /*! @brief Build aggregates. */
  RCP<Aggregates> Build(const Graph& graph) const
  {
    RCP<Aggregates> aggregates = CoarsenUncoupled(graph); //TODO: remove options_ arg.
    return aggregates;
  }
  //@}
  
private:
  //! Aggregation options
  int printFlag_;
  Ordering ordering_;                /**<  natural, random, graph           */
  int      minNodesPerAggregate_;    /**<  aggregate size control           */
  int      maxNeighAlreadySelected_; /**<  complexity control               */
    
  // unused 
  // aggregation algorithm type
  // std::string Algorithm_;

  //! coalesce and drop factory
  RCP<FactoryBase> coalesceDropFact_;

  //! user-defined graph label
    std::string graphName_;//TODO unused?

  //! @name Aggregation methods.
  //@{

  /*! @brief Local aggregation.
   */
  RCP<Aggregates>
  CoarsenUncoupled(Graph const & graph) const
  {
    /* Create Aggregation object */
    const std::string name = "Uncoupled";
    my_size_t nAggregates = 0;
    RCP<Aggregates> aggregates = rcp(new Aggregates(graph, name));

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
          
      Teuchos::ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0); // output only: contents ignored
          
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
                std::cout << "Error: couldn't allocate memory for supernode! " << length << std::endl;
                exit(1); //TODO: exception instead
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
                  aggregates->SetIsRoot(iNode);
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
    aggregates->SetNumAggregates(nAggregates);

    /* Verbose */
    // TODO: replace AllReduce by Reduce to proc 0
    if (GetPrintFlag() < 7) { //FIXME
      const RCP<const Teuchos::Comm<int> > & comm = graph.GetComm();
      int myPid = comm->getRank();
          
      {
        GO localReady=0, globalReady;
            
        // Compute 'localReady'
        for ( my_size_t i = 0; i < nRows; ++i ) 
          if ( aggStat[i] == READY ) localReady++;
            
        // Compute 'globalReady'
        sumAll(comm, localReady, globalReady);
            
        if (myPid == 0 && globalReady > 0)
          std::cout << "Aggregation(UC) : Phase 1 (WARNING) - " << globalReady << " READY nodes left\n" << std::endl;
      }
          
      {
        // Compute 'localSelected'
        LO localSelected=0;
        for ( my_size_t i = 0; i < nRows; ++i ) 
          if ( aggStat[i] == SELECTED ) localSelected++;
            
        // Compute 'globalSelected'
        GO globalSelected; sumAll(comm, (GO)localSelected, globalSelected);

        // Compute 'globalNRows'
        GO globalNRows; sumAll(comm, (GO)nRows, globalNRows);
            
        if (myPid == 0)
          std::cout << "Aggregation(UC) : Phase 1 - nodes aggregated = " << globalSelected << " (" << globalNRows << ")" << std::endl;
      }
          
      {
        GO nAggregatesGlobal; sumAll(comm, (GO)nAggregates, nAggregatesGlobal);
        if (myPid == 0)
          std::cout << "Aggregation(UC) : Phase 1 - total aggregates = " << nAggregatesGlobal << std::endl;
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
  } //CoarsenUncoupled

    //! @name Utilities
    //@{

    /*! @brief Utility to take a list of integers and reorder them randomly (by using a local permutation).
      @param list On input, a bunch of integers. On output, the same integers in a different order
      that is determined randomly.
    */
  void RandomReorder(Teuchos::ArrayRCP<LO> list) const
  { //TODO: replace int
    int n = list.size();
    for(int i=0; i<n-1; i++) {
      std::swap(list[i], list[RandomOrdinal(i,n-1)]);
    }
  } 

  /*! @brief Generate a random number in the range [min, max] */
  int RandomOrdinal(int min, int max) const
  {
    return min + static_cast<int>((max-min+1) * (static_cast<double>(std::rand()) / (RAND_MAX + 1.0)));
  }

  //@}
  
}; //class LocalAggregationFactory

} //namespace MueLu

#define MUELU_LOCALAGGREGATIONFACTORY_SHORT
#endif //ifndef MUELU_LOCALAGGREGATIONFACTORY_HPP
