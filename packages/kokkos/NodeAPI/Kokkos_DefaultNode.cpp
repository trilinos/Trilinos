#include "Kokkos_DefaultNode.hpp"
#include <Teuchos_ParameterList.hpp>

Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> Kokkos::DefaultNode::node_ = Teuchos::null;

namespace Kokkos {

  Teuchos::RCP<DefaultNode::DefaultNodeType> DefaultNode::getDefaultNode() {
    if (node_ == Teuchos::null) {
      Teuchos::ParameterList pl;
#ifdef HAVE_KOKKOS_THREADPOOL
      pl.set<int>("Num Threads",1);
      node_ = Teuchos::rcp<TPINode>(new TPINode(pl));
#else
#  ifdef HAVE_KOKKOS_TBB
      pl.set<int>("Num Threads",0);
      node_ = Teuchos::rcp<TBBNode>(new TBBNode(pl));
#  else
      node_ = Teuchos::rcp<SerialNode>(new SerialNode(pl));
#  endif
#endif
    }
    return node_;
  }

}
