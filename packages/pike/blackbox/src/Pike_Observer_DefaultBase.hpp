#ifndef PIKE_OBSERVER_DEFAULT_BASE_HPP
#define PIKE_OBSERVER_DEFAULT_BASE_HPP

#include "Pike_Observer.hpp"

namespace pike {

  class ObserverDefaultBase : public pike::Observer {

  public:

    virtual ~ObserverDefaultBase() {}

    virtual void observeBeginSolve(const Solver& solver);

    virtual void observeEndSolve(const Solver& solver);

    virtual void observeBeginStep(const Solver& solver);

    virtual void observeEndStep(const Solver& solver);

    virtual void observeConvergedSolve(const Solver& solver);

    virtual void observeFailedSolve(const Solver& solver);

  };

}

#endif
