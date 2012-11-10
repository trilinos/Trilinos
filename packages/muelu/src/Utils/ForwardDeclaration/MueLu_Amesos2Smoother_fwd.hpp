#ifndef MUELU_AMESOS2SMOOTHER_FWDDECL_HPP
#define MUELU_AMESOS2SMOOTHER_FWDDECL_HPP

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MUELU_AMESOS2

namespace MueLu {
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class Amesos2Smoother;
}

#ifndef MUELU_AMESOS2SMOOTHER_SHORT
#define MUELU_AMESOS2SMOOTHER_SHORT
#endif

#endif // HAVE_MUELU_AMESOS2
#endif // MUELU_AMESOS2SMOOTHER_FWDDECL_HPP
