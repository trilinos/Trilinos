#ifndef FAKEKOKKOS_DEFAULT_NODE_HPP_
#define FAKEKOKKOS_DEFAULT_NODE_HPP_

#include <Teuchos_RCP.hpp>

namespace Kokkos {

  class DefaultNode {
    public:

      typedef int DefaultNodeType;

    static Teuchos::RCP<DefaultNodeType> getDefaultNode() { return Teuchos::null; }
  };

}

#endif
