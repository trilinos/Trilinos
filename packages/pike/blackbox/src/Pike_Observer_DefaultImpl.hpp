#ifndef NOX_OBSERVER_DEFAULTIMPL_HPP
#define NOX_OBSERVER_DEFAULTIMPL_HPP

namespace nox {

  template<typename Solver>
  class ObserverDefaultImpl {

    virtual void observeBeginSolve(const Solver& solver);

    virtual void observeEndSolve(const Solver& solver);

    virtual void observeBeginStep(const Solver& solver);

    virtual void observeEndStep(const Solver& solver);

    virtual void observeConvergedSolve(const Solver& solver);

    virtual void observeFailedSolve(const Solver& solver);

  };

}

#endif
