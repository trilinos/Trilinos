#include "Pike_SolverObserver_DefaultBase.hpp"

namespace pike {

  void ObserverDefaultBase::observeInitialization(const pike::Solver& solver) {}
  
  void ObserverDefaultBase::observeFinalization(const pike::Solver& solver) {}

  void ObserverDefaultBase::observeBeginSolve(const Solver& solver) {}

  void ObserverDefaultBase::observeEndSolve(const Solver& solver) {}
  
  void ObserverDefaultBase::observeBeginStep(const Solver& solver) {}
  
  void ObserverDefaultBase::observeEndStep(const Solver& solver) {}
  
  void ObserverDefaultBase::observeConvergedSolve(const Solver& solver) {}
  
  void ObserverDefaultBase::observeFailedSolve(const Solver& solver) {}

}
