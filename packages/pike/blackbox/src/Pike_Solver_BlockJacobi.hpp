#ifndef PIKE_SOLVER_BLOCK_JACOBI_HPP
#define PIKE_SOLVER_BLOCK_JACOBI_HPP

#include "Pike_Solver_DefaultBase.hpp"

namespace pike {

  class BlockJacobi : public pike::SolverDefaultBase {
    
  public:

    BlockJacobi();

    void stepImplementation();
    
  };

}

#endif
