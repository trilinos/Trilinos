#include "Pike_Observer_DefaultImpl.hpp"

namespace pike {

  void ObserverDefaultImpl::observeBeginSolve(const Solver& solver)
  { }

  void ObserverDefaultImpl::observeEndSolve(const Solver& solver)
  { }
  
  void ObserverDefaultImpl::observeBeginStep(const Solver& solver)
  { }
  
  void ObserverDefaultImpl::observeEndStep(const Solver& solver)
  { }
  
  void ObserverDefaultImpl::observeConvergedSolve(const Solver& solver)
  { }
  
  void ObserverDefaultImpl::observeFailedSolve(const Solver& solver)
  { }

}
