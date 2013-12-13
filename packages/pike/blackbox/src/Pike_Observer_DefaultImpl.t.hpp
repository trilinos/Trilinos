#ifndef NOX_OBSERVER_DEFAULTIMPL_T_HPP
#define NOX_OBSERVER_DEFAULTIMPL_T_HPP

namespace nox {

  template<typename Solver>
  void ObserverDefaultImpl<Solver>::observeBeginSolve(const Solver& solver)
  { }

  template<typename Solver>
  void ObserverDefaultImpl<Solver>::observeEndSolve(const Solver& solver)
  { }
  
  template<typename Solver>
  void ObserverDefaultImpl<Solver>::observeBeginStep(const Solver& solver)
  { }
  
  template<typename Solver>
  void ObserverDefaultImpl<Solver>::observeEndStep(const Solver& solver)
  { }
  
  template<typename Solver>
  void ObserverDefaultImpl<Solver>::observeConvergedSolve(const Solver& solver)
  { }
  
  template<typename Solver>
  void ObserverDefaultImpl<Solver>::observeFailedSolve(const Solver& solver)
  { }

}

#endif
