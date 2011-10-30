#ifndef MUELU_LEFTOVERAGGREGATIONALGORITHM_HPP_DEF
#define MUELU_LEFTOVERAGGREGATIONALGORITHM_HPP_DEF

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include <assert.h>
#include <math.h>
#include <vector>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_UCAggregationCommHelper.hpp"

namespace MueLu {

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class LeftoverAggregationAlgorithm : public BaseClass {
#include "MueLu_UseShortNamesOrdinal.hpp"

    typedef GO global_size_t; 
    typedef LO my_size_t;     

  public:

    LeftoverAggregationAlgorithm();

    virtual ~LeftoverAggregationAlgorithm();

    void SetMinNodesPerAggregate(int minNodesPerAggregate);
    void SetPhase3AggCreation(double phase3AggCreation);

    double GetPhase3AggCreation() const;
    int GetMinNodesPerAggregate() const;

    void AggregateLeftovers(Graph const &graph, Aggregates &aggregates) const;

    void RootCandidates(my_size_t nVertices, ArrayView<const LO> & vertex2AggId, Graph const &graph, ArrayRCP<LO> &candidates, my_size_t &nCandidates, global_size_t &nCandidatesGlobal) const;

    int RemoveSmallAggs(Aggregates& aggregates, int min_size, RCP<Xpetra::Vector<double,LO,GO,NO> > & distWeights, const MueLu::UCAggregationCommHelper<LO,GO,NO,LMO> & myWidget) const;
  
  private:
    double phase3AggCreation_;
    int    minNodesPerAggregate_;

  };

}

#define MUELU_LEFTOVERAGGREGATIONALGORITHM_SHORT

#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION

#endif //ifndef MUELU_LEFTOVERAGGREGATIONALGORITHM_HPP
