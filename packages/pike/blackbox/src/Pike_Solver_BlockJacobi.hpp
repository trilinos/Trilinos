#ifndef PIKE_SOLVER_BLOCK_JACOBI_HPP
#define PIKE_SOLVER_BLOCK_JACOBI_HPP

#include "Pike_Solver_DefaultImplementation.hpp"

namespace pike {

  class BlockJacobi : public pike::SolverDefaultImplementation {
    
  public:

    void stepImplementation();
    
  };

}

#endif
