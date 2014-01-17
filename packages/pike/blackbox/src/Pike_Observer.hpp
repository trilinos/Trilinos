#ifndef PIKE_OBSERVER_HPP
#define PIKE_OBSERVER_HPP

namespace pike {

  class Solver;

  /** \brief Observer design pattern for monitoring pike::Solver methods */ 
  class Observer {
    
  public:

    virtual ~Observer() {}

    virtual void observeBeginSolve(const pike::Solver& solver) = 0;

    virtual void observeEndSolve(const pike::Solver& solver) = 0;

    virtual void observeBeginStep(const pike::Solver& solver) = 0;

    virtual void observeEndStep(const pike::Solver& solver) = 0;

    virtual void observeConvergedSolve(const pike::Solver& solver) = 0;

    virtual void observeFailedSolve(const pike::Solver& solver) = 0;

  };

}

#endif
