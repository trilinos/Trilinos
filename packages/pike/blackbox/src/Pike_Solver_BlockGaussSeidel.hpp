#ifndef PIKE_SOLVER_BLOCK_GAUSS_SEIDEL_HPP
#define PIKE_SOLVER_BLOCK_GAUSS_SEIDEL_HPP

#include "Pike_Solver_DefaultImplementation.hpp"

namespace pike {

  class BlockGaussSeidel : public pike::SolverDefaultImplementation {
    
  public:
    
    void stepImplementation();
    
  };

}

#endif
