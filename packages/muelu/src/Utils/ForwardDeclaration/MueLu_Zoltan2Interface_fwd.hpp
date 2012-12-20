#ifndef MUELU_ZOLTAN2INTERFACE_FWD_HPP
#define MUELU_ZOLTAN2INTERFACE_FWD_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)

namespace MueLu {
  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class Zoltan2Interface;
}

#ifndef MUELU_ZOLTAN2INTERFACE_SHORT
#define MUELU_ZOLTAN2INTERFACE_SHORT
#endif // MUELU_ZOLTAN2INTERFACE_SHORT

#endif // HAVE_MUELU_ZOLTAN2 && HAVE_MPI

#endif // MUELU_ZOLTAN2INTERFACE_FWD_HPP
