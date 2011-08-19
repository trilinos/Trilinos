#ifndef MUELU_COALESCEDROPFACTORY_HPP
#define MUELU_COALESCEDROPFACTORY_HPP

#include "Xpetra_Operator.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"

namespace MueLu {

/*!
  @class CoalesceDropFactory
  @brief Factory for creating a graph base on a given matrix.

  Factory for creating graphs from matrices with entries selectively dropped.
  
  - TODO This factory is very incomplete.
  - TODO The Build method simply builds the matrix graph with no dropping.
  - TODO Derive from TwoLevelFactoryBase so that UCAggregationFactory does not rely on Scalar template.
         This entails removing Scalar dependence in Level class.
  
*/

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class CoalesceDropFactory : public Teuchos::Describable {

#include "MueLu_UseShortNames.hpp"

private:

public:

  //! @name Constructors/Destructors.
  //@{
  CoalesceDropFactory() {};

  //!Destructor
  virtual ~CoalesceDropFactory() {}
  //@}

  void Build(Level &currentLevel) {
    Teuchos::RCP<Operator> A = currentLevel.template Get< Teuchos::RCP<Operator> >("A");
    Teuchos::RCP<Graph> graph = Teuchos::rcp(new Graph(A->getCrsGraph(), "graph of A"));
    currentLevel.Set("Graph",graph);
  } //Build

}; //class CoalesceDropFactory

} //namespace MueLu

#define MUELU_COALESCEDROPFACTORY_SHORT

#endif //ifndef MUELU_COALESCEDROPFACTORY_HPP
