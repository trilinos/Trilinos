#ifndef PIKE_OBSERVER_LOGGER_HPP
#define PIKE_OBSERVER_LOGGER_HPP

#include "Pike_SolverObserver.hpp"
#include "Teuchos_RCP.hpp"
#include <vector>
#include <string>

namespace pike {

  class LoggerObserver : public pike::SolverObserver {

  public:
    
    LoggerObserver();

    void setLog(const Teuchos::RCP<std::vector<std::string> >& log);

    Teuchos::RCP<const std::vector<std::string> > getLog() const;

    Teuchos::RCP<std::vector<std::string> > getNonConstLog() const;

    void observeInitialization(const pike::Solver& solver);

    void observeFinalization(const pike::Solver& solver);

    void observeBeginSolve(const Solver& solver);

    void observeEndSolve(const Solver& solver);

    void observeBeginStep(const Solver& solver);

    void observeEndStep(const Solver& solver);

    void observeConvergedSolve(const Solver& solver);

    void observeFailedSolve(const Solver& solver);

  private:

    Teuchos::RCP<std::vector<std::string> > log_;
  };

  /** \brief Non-member ctor
      \relates LoggerObserver
   */
  Teuchos::RCP<pike::LoggerObserver> loggerObserver();

}

#endif
