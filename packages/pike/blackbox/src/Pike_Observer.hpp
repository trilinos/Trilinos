#ifndef PIKE_OBSERVER_HPP
#define PIKE_OBSERVER_HPP

#include "Pike_BlackBox_config.hpp"

namespace pike {

  class Solver;

  class Observer {
    
  public:

    virtual ~Observer();

    virtual void observeBeginSolve(const pike::Solver& solver) = 0;

    virtual void observeEndSolve(const pike::Solver& solver) = 0;

    virtual void observeBeginStep(const pike::Solver& solver) = 0;

    virtual void observeEndStep(const pike::Solver& solver) = 0;

    virtual void observeConvergedSolve(const pike::Solver& solver) = 0;

    virtual void observeFailedSolve(const pike::Solver& solver) = 0;

  };

}

#endif
