/*
 * MueLu_DemoFactory_def.hpp
 *
 *  Created on: 06.12.2011
 *      Author: tobias
 */

#ifndef MUELU_DEMOFACTORY_DEF_HPP_
#define MUELU_DEMOFACTORY_DEF_HPP_


#include "MueLu_ThresholdAFilterFactory_decl.hpp"

#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

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

#define MUELU_DEMOFACTORY_SHORT
} // namespace MueLu

#endif /* MUELU_DEMOFACTORY_DEF_HPP_ */
