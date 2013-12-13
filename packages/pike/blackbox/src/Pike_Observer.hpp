#ifndef NOX_OBSERVER_HPP
#define NOX_OBSERVER_HPP

namespace nox {

  template<typename Solver>
  class Observer {

    virtual ~Observer();

    virtual void observeBeginSolve(const Solver& solver) = 0;

    virtual void observeEndSolve(const Solver& solver) = 0;

    virtual void observeBeginStep(const Solver& solver) = 0;

    virtual void observeEndStep(const Solver& solver) = 0;

    virtual void observeConvergedSolve(const Solver& solver) = 0;

    virtual void observeFailedSolve(const Solver& solver) = 0;

  };

}

#endif
