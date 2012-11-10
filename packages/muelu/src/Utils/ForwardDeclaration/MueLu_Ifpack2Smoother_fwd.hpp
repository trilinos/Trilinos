#ifndef MUELU_IFPACK2SMOOTHER_FWDDECL_HPP
#define MUELU_IFPACK2SMOOTHER_FWDDECL_HPP

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MUELU_IFPACK2

namespace MueLu {
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class Ifpack2Smoother;
}

#ifndef MUELU_IFPACK2SMOOTHER_SHORT
#define MUELU_IFPACK2SMOOTHER_SHORT
#endif

#endif // HAVE_MUELU_IFPACK2
#endif // MUELU_IFPACK2SMOOTHER_FWDDECL_HPP
