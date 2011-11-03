#ifndef MUELU_SMOOTHERPROTOTYPE_DEF_HPP
#define MUELU_SMOOTHERPROTOTYPE_DEF_HPP

#include "MueLu_SmootherPrototype_decl.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SmootherPrototype() : isSetup_(false) {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~SmootherPrototype() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::IsSetup() const { return isSetup_; }
    
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::IsSetup(bool const &ToF) { isSetup_ = ToF; }

} // namespace MueLu

//TODO: private copy constructor
//TODO: update comments
#endif // MUELU_SMOOTHERPROTOTYPE_DEF_HPP
