#ifndef MUELU_ZOLTANINTERFACE_FWD_HPP
#define MUELU_ZOLTANINTERFACE_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)

namespace MueLu {
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class ZoltanInterface;
}

#ifndef MUELU_ZOLTANINTERFACE_SHORT
#define MUELU_ZOLTANINTERFACE_SHORT
#endif // MUELU_ZOLTANINTERFACE_SHORT

#endif // HAVE_MUELU_ZOLTAN && HAVE_MPI

#endif // MUELU_ZOLTANINTERFACE_FWD_HPP
