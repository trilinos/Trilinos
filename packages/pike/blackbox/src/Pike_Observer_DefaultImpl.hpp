#ifndef PIKE_OBSERVER_DEFAULT_IMPL_HPP
#define PIKE_OBSERVER_DEFAULT_IMPL_HPP

#include "Pike_Observer.hpp"

namespace pike {

  class ObserverDefaultImpl : public pike::Observer {

  public:

    virtual ~ObserverDefaultImpl() {}

    virtual void observeBeginSolve(const Solver& solver);

    virtual void observeEndSolve(const Solver& solver);

    virtual void observeBeginStep(const Solver& solver);

    virtual void observeEndStep(const Solver& solver);

    virtual void observeConvergedSolve(const Solver& solver);

    virtual void observeFailedSolve(const Solver& solver);

  };

}

#endif
