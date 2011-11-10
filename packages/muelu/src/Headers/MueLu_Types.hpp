#ifndef MUELU_TYPES_HPP
#define MUELU_TYPES_HPP

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {
  enum CycleType {
    VCYCLE,
    WCYCLE
  };

  enum PreOrPost {
    PRE,
    POST,
    BOTH
  };
}

#endif //ifndef MUELU_TYPES_HPP
