#ifndef MUELU_HIERARCHYFACTORY_DEF_HPP
#define MUELU_HIERARCHYFACTORY_DEF_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "MueLu_HierarchyFactory_decl.hpp"

#include "MueLu_Hierarchy.hpp"

#include "MueLu_RAPFactory.hpp" //TMP

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  HierarchyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::HierarchyFactory(Teuchos::ParameterList & paramList) {
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  HierarchyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::HierarchyFactory(const std::string & xmlFileName) {
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);
    SetParameterList(*paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void HierarchyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetParameterList(Teuchos::ParameterList & paramList) {
    std::cout << "Parameter List:" << std::endl
              << paramList
              << std::endl;
  }
    
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > HierarchyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CreateHierarchy() const {
    return rcp(new Hierarchy());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void HierarchyFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupHierarchy(MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> & H) const {
    FactoryManager M;                         // -
    M.SetFactory("A", rcp(new RAPFactory())); // TODO: to be remove, but will require some work
    H.Setup(M);                               // -
  }
  
} // namespace MueLu

#endif // MUELU_HIERARCHYFACTORY_DEF_HPP
