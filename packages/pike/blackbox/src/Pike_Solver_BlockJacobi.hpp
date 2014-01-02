#ifndef PIKE_SOLVER_BLOCK_JACOBI_HPP
#define PIKE_SOLVER_BLOCK_JACOBI_HPP

#include "Pike_Solver_DefaultImpl.hpp"

namespace pike {

  class BlockJacobi : public pike::SolverDefaultImpl {
    
  public:

    void stepImplementation();
    
  };

}

#endif
