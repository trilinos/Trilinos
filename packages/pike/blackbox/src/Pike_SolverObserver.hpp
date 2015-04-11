#ifndef PIKE_SOLVER_OBSERVER_HPP
#define PIKE_SOLVER_OBSERVER_HPP

namespace pike {

  class Solver;

  /** \brief Observer design pattern for monitoring pike::Solver methods */ 
  class SolverObserver {
    
  public:

    virtual ~SolverObserver() {}

    virtual void observeInitialization(const pike::Solver& solver) = 0;

    virtual void observeFinalization(const pike::Solver& solver) = 0;

    virtual void observeBeginSolve(const pike::Solver& solver) = 0;

    virtual void observeEndSolve(const pike::Solver& solver) = 0;

    virtual void observeBeginStep(const pike::Solver& solver) = 0;

    virtual void observeEndStep(const pike::Solver& solver) = 0;

    virtual void observeConvergedSolve(const pike::Solver& solver) = 0;

    virtual void observeFailedSolve(const pike::Solver& solver) = 0;

  };

}

#endif
