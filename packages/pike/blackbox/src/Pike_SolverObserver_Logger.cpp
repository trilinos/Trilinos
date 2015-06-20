#include "Pike_SolverObserver_Logger.hpp"
#include "Pike_Solver.hpp"

namespace pike {

  LoggerObserver::LoggerObserver()
  {
    log_ = Teuchos::rcp(new std::vector<std::string>);
  }

  void LoggerObserver::setLog(const Teuchos::RCP<std::vector<std::string> >& log)
  {
    log_ = log;
  }

  Teuchos::RCP<const std::vector<std::string> > LoggerObserver::getLog() const
  {
    return log_;
  }
  
  Teuchos::RCP<std::vector<std::string> > LoggerObserver::getNonConstLog() const
  {
    return log_;
  }
  
  void LoggerObserver::observeInitialization(const pike::Solver& solver)
  {
    log_->push_back(solver.name()+": observeInitialization");
  }

  void LoggerObserver::observeFinalization(const pike::Solver& solver)
  {
    log_->push_back(solver.name()+": observeFinalization");
  }

  void LoggerObserver::observeBeginSolve(const Solver& solver)
  { 
    log_->push_back(solver.name()+": observeBeginSolve");
  }

  void LoggerObserver::observeEndSolve(const Solver& solver)
  {
    log_->push_back(solver.name()+": observeEndSolve");
  }

  void LoggerObserver::observeBeginStep(const Solver& solver)
  {
    log_->push_back(solver.name()+": observeBeginStep");
  }

  void LoggerObserver::observeEndStep(const Solver& solver)
  {
    log_->push_back(solver.name()+": observeEndStep");
  }

  void LoggerObserver::observeConvergedSolve(const Solver& solver)
  {
    log_->push_back(solver.name()+": observeConvergedSolve");
  }

  void LoggerObserver::observeFailedSolve(const Solver& solver)
  {
    log_->push_back(solver.name()+": observeFailedSolve");
  }

  Teuchos::RCP<pike::LoggerObserver> loggerObserver()
  {
    return Teuchos::rcp(new pike::LoggerObserver);
  }
  
}
