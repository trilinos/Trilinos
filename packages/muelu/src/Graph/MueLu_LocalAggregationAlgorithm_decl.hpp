#ifndef MUELU_LOCALAGGREGATIONALGORITHM_DECL_HPP
#define MUELU_LOCALAGGREGATIONALGORITHM_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_LinkedList.hpp"

#include "MueLu_Monitor.hpp"

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
    @class LocalAggregationAlgorithm class.
    @brief Algorithm for coarsening a graph with uncoupled aggregation.

    This method has two phases.  The first is a local clustering algorithm.  The second creates aggregates
    that can include unknowns from more than one process.

  */

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class LocalAggregationAlgorithm : public BaseClass {
#include "MueLu_UseShortNamesOrdinal.hpp"

    typedef GO global_size_t; //TODO
    typedef LO my_size_t; //TODO


  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    LocalAggregationAlgorithm(RCP<FactoryBase> const &graphFact=Teuchos::null)
    ;

    //! Destructor.
    virtual ~LocalAggregationAlgorithm() ;

    //@}

    //! @name Set/get methods.
    //@{

    void SetOrdering(Ordering ordering)                          ;
    void SetMinNodesPerAggregate(int minNodesPerAggregate)       ;
    void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) ;
    
    Ordering GetOrdering()                const ;
    int      GetMinNodesPerAggregate()    const ;
    int      GetMaxNeighAlreadySelected() const ;

    //@}

    //! @name Aggregation methods.
    //@{

    /*! @brief Local aggregation. */
    void CoarsenUncoupled(Graph const & graph, Aggregates & aggregates) const
    ; // CoarsenUncoupled

  private:
    //! Aggregation options (TODO: Teuchos::ParameterList?)
    Ordering ordering_;                /**<  natural, random, graph           */
    int      minNodesPerAggregate_;    /**<  aggregate size control           */
    int      maxNeighAlreadySelected_; /**<  complexity control               */

    //! @name Utilities
    //@{

    /*! @brief Utility to take a list of integers and reorder them randomly (by using a local permutation).
      @param list On input, a bunch of integers. On output, the same integers in a different order
      that is determined randomly.
    */
    void RandomReorder(Teuchos::ArrayRCP<LO> list) const
    ; 

    /*! @brief Generate a random number in the range [min, max] */
    int RandomOrdinal(int min, int max) const
    ;

    //@}
  
  }; //class LocalAggregationFactory

} //namespace MueLu

#define MUELU_LOCALAGGREGATIONALGORITHM_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_LOCALAGGREGATIONALGORITHM_DECL_HPP
