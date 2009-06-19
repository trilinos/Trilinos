#ifndef KOKKOS_DEFAULT_NODE_HPP_
#define KOKKOS_DEFAULT_NODE_HPP_

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif

#include <Teuchos_RCP.hpp>

namespace Kokkos {

  class DefaultNode {
    public:
#ifdef HAVE_KOKKOS_TBB
      typedef TBBNode DefaultNodeType;
#else
      typedef SerialNode DefaultNodeType;
#endif

      static DefaultNodeType & getDefaultNode();

    private:
      static Teuchos::RCP<DefaultNodeType> node_;
  };
}

#endif
