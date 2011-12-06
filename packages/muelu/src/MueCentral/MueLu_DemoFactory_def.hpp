#ifndef MUELU_DEMOFACTORY_DEF_HPP
#define MUELU_DEMOFACTORY_DEF_HPP

#include "MueLu_DemoFactory_decl.hpp"

// #include <Xpetra_Operator.hpp>

#include "MueLu_Level.hpp"
// #include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  DemoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DemoFactory()
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  DemoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~DemoFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DemoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    // TODO: declare input for factory
    //currentLevel.DeclareInput(varName_,factory_,this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DemoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    // TODO: implement factory
  }

} // namespace MueLu

#define MUELU_DEMOFACTORY_SHORT
#endif // MUELU_DEMOFACTORY_DEF_HPP
