#ifndef MUELU_LOCALAGGREGATIONFACTORY_HPP_DECL
#define MUELU_LOCALAGGREGATIONFACTORY_HPP_DECL

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

    enum Ordering {
      NATURAL = 0, 
      RANDOM  = 1, 
      GRAPH   = 2  
    };
  }

  using namespace AggOptions;

  enum NodeState {
    READY   = -11,
    NOTSEL  = -12,
    SELECTED = -13
  };

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class LocalAggregationAlgorithm : public BaseClass {
#include "MueLu_UseShortNamesOrdinal.hpp"

    typedef GO global_size_t; //TODO
    typedef LO my_size_t; //TODO


  public:
  
    LocalAggregationAlgorithm(RCP<FactoryBase> const &graphFact=Teuchos::null);

    virtual ~LocalAggregationAlgorithm();

    void SetOrdering(Ordering ordering);                          
    void SetMinNodesPerAggregate(int minNodesPerAggregate)  ;     
    void SetMaxNeighAlreadySelected(int maxNeighAlreadySelected) ;
    
    Ordering GetOrdering()                const;
    int      GetMinNodesPerAggregate()    const;
    int      GetMaxNeighAlreadySelected() const;

    void CoarsenUncoupled(Graph const & graph, Aggregates & aggregates) const;

  private:
    Ordering ordering_;                /**<  natural, random, graph           */
    int      minNodesPerAggregate_;    /**<  aggregate size control           */
    int      maxNeighAlreadySelected_; /**<  complexity control               */

    void RandomReorder(Teuchos::ArrayRCP<LO> list) const;

    int RandomOrdinal(int min, int max) const;
  
  };

}

#define MUELU_LOCALAGGREGATIONALGORITHM_SHORT
#endif //ifndef MUELU_LOCALAGGREGATIONFACTORY_HPP_DECL
