#include "Kokkos_DefaultNode.hpp"

Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> Kokkos::DefaultNode::node_ = Teuchos::null;

namespace Kokkos {

  DefaultNode::DefaultNodeType &DefaultNode::getDefaultNode() {
    if (!node_.get()) {
#ifdef HAVE_KOKKOS_TBB
      node_ = Teuchos::rcp<TBBNode>(new TBBNode(0));
#else
      node_ = Teuchos::rcp<SerialNode>(new SerialNode());
#endif
    }
    return *node_;
  }

}
