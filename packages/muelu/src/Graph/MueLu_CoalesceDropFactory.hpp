#ifndef MUELU_COALESCEDROPFACTORY_HPP
#define MUELU_COALESCEDROPFACTORY_HPP

#include "Xpetra_Operator.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"

namespace MueLu {

  /*!
    @class CoalesceDropFactory
    @brief Factory for creating a graph base on a given matrix.

    Factory for creating graphs from matrices with entries selectively dropped.
  
    - TODO This factory is very incomplete.
    - TODO The Build method simply builds the matrix graph with no dropping.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class CoalesceDropFactory : public SingleLevelFactoryBase {

#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    CoalesceDropFactory(RCP<FactoryBase> AFact = Teuchos::null)
      : AFact_(AFact)
    { }

    //! Destructor
    virtual ~CoalesceDropFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const {
      currentLevel.Input("A", AFact_());
    }

    //@}

    bool Build(Level &currentLevel) const {
      RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_);

      RCP<Graph> graph = rcp(new Graph(A->getCrsGraph(), "Graph of A"));

      currentLevel.Set("Graph", graph, this);

      return true; //??

    } // Build

  private:
    //! A Factory
    RCP<FactoryBase> AFact_;

  }; //class CoalesceDropFactory

} //namespace MueLu

#define MUELU_COALESCEDROPFACTORY_SHORT

#endif //ifndef MUELU_COALESCEDROPFACTORY_HPP
