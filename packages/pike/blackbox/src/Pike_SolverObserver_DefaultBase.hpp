#ifndef PIKE_OBSERVER_DEFAULT_BASE_HPP
#define PIKE_OBSERVER_DEFAULT_BASE_HPP

#include "Pike_SolverObserver.hpp"

namespace pike {

  class ObserverDefaultBase : public pike::SolverObserver {

  public:

    virtual ~ObserverDefaultBase() {}

    virtual void observeInitialization(const pike::Solver& solver);

    virtual void observeFinalization(const pike::Solver& solver);

    virtual void observeBeginSolve(const Solver& solver);

    virtual void observeEndSolve(const Solver& solver);

    virtual void observeBeginStep(const Solver& solver);

    virtual void observeEndStep(const Solver& solver);

    virtual void observeConvergedSolve(const Solver& solver);

    virtual void observeFailedSolve(const Solver& solver);

  };

}

#endif
