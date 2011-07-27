#ifndef KOKKOS_DEFAULT_NODE_HPP_
#define KOKKOS_DEFAULT_NODE_HPP_

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif

#include <Teuchos_RCP.hpp>

namespace Kokkos {

  /** \brief Class to specify %Kokkos default node type and instantiate the default node.
      \ingroup kokkos_node_api
    */
  class DefaultNode {
    public:
#if defined(HAVE_KOKKOS_THREADPOOL)
      typedef TPINode DefaultNodeType;
#elif defined(HAVE_KOKKOS_TBB)
      typedef TBBNode DefaultNodeType;
#else
      //! Typedef specifying the default node type.
      typedef SerialNode DefaultNodeType;
#endif

      //! \brief Return a pointer to the default node.
      static RCP<DefaultNodeType> getDefaultNode();

    private:
      static RCP<DefaultNodeType> node_;
  };

}

#endif
