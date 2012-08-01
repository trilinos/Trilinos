/*
 * MueLu_CheapAggregationAlgorithm_decl.hpp
 *
 *  Created on: Jul 25, 2012
 *      Author: wiesner
 */

#ifndef MUELU_CHEAPAGGREGATIONALGORITHM_DECL_HPP_
#define MUELU_CHEAPAGGREGATIONALGORITHM_DECL_HPP_

#include <vector>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_CheapAggregationAlgorithm_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Graph_fwd.hpp"

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

  class Aggregate {
  public:
    int length;                   // current size of aggregate
    int maxLength;                // max size of aggregate
    int index;                    // local aggregate id
    std::vector<int> list;  // list of node ids in aggregate
  };

  /* In the algorithm, aggStat[]=READY/NOTSEL/SELECTED indicates whether a node has been aggregated. */
  enum NodeState {
    READY   = -11,   /* indicates that a node is available to be */
    /* selected as a root node of an aggregate  */

    NOTSEL  = -12,   /* indicates that a node has been rejected  */
    /* as a root node. This could perhaps be    */
    /* because if this node had been selected a */
    /* small aggregate would have resulted.     */

    SELECTED = -13,  /* indicates that a node has been assigned  */
    /* to an aggregate.                         */

    BDRY = -15, /* indicates that a node is a Dirichlet bdry node */

    READY_1PT = -16, /* indicates that a node is ready to be aggregates (as a single point aggregate only) */
    SELECTED_1PT = -17 /* indicates that a 1pt aggregate node is already aggregated */
  };

  /*!
    @class CheapAggregationAlgorithm class.
    @brief Algorithm for coarsening a graph with uncoupled aggregation.
  */

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class CheapAggregationAlgorithm : public BaseClass {
#undef MUELU_CHEAPAGGREGATIONALGORITHM_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

    typedef GO global_size_t; //TODO
    typedef LO my_size_t; //TODO

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    CheapAggregationAlgorithm(RCP<const FactoryBase> const &graphFact = Teuchos::null);

    //! Destructor.
    virtual ~CheapAggregationAlgorithm() { }

    //@}

    //! @name Set/get methods.
    //@{

    void SetOrdering(Ordering ordering)                          { ordering_                = ordering;                }
    void SetMinNodesPerAggregate(int minNodesPerAggregate)       { minNodesPerAggregate_    = minNodesPerAggregate;    }
    void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) { maxNeighAlreadySelected_ = maxNeighAlreadySelected; }

    Ordering GetOrdering()                const { return ordering_;                }
    int      GetMinNodesPerAggregate()    const { return minNodesPerAggregate_;    }
    int      GetMaxNeighAlreadySelected() const { return maxNeighAlreadySelected_; }

    //@}

    //! @name Aggregation methods.
    //@{

    /*! @brief Local aggregation. */


    LocalOrdinal Phase1(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat) const; // local uncoupled coarsening (Phase 1)
    LocalOrdinal Phase1a(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat, Teuchos::ArrayRCP<NodeState> & coarse_aggStat) const; // local uncoupled coarsening (Phase 1a)
    LocalOrdinal Phase2_maxlink(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat) const; // local uncoupled coarsening (Phase 2 [max_link])

    LocalOrdinal Phase3(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat) const; // local uncoupled coarsening (Phase 3)
    LocalOrdinal Phase4(Graph const & graph, Aggregates & aggregates, Teuchos::ArrayRCP<NodeState> & aggStat) const; // local uncoupled coarsening (Phase 4)

  private:
    //! Aggregation options (TODO: Teuchos::ParameterList?)
    Ordering ordering_;                /**<  natural, random, graph           */
    int      minNodesPerAggregate_;    /**<  aggregate size control           */
    int      maxNeighAlreadySelected_; /**<  complexity control               */

    /*! @brief Utility to take a list of integers and reorder them randomly (by using a local permutation).
      @param list On input, a bunch of integers. On output, the same integers in a different order
      that is determined randomly.
    */
    void RandomReorder(Teuchos::ArrayRCP<LO> list) const; 

    /*! @brief Generate a random number in the range [min, max] */
    int RandomOrdinal(int min, int max) const;

    //@}

  }; //class CheapAggregationAlgorithm

} //namespace MueLu

#endif /* MUELU_CHEAPAGGREGATIONALGORITHM_DECL_HPP_ */
