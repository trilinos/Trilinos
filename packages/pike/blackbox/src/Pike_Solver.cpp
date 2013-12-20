#include "Pike_Solver.hpp"

namespace pike {

  std::ostream& operator<<(std::ostream& os, const pike::Solver& solver)
  {
    solver.describe(os,solver.getVerbLevel());
    return os;
  }

}
