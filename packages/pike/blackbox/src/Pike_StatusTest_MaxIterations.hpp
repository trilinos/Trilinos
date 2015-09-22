
#ifndef PIKE_STATUS_TESTS_MAX_ITERATIONS_HPP
#define PIKE_STATUS_TESTS_MAX_ITERATIONS_HPP

#include "Pike_StatusTest.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include <iostream>

namespace pike {

  class MaxIterations : 
    public pike::StatusTest,
    public Teuchos::ParameterListAcceptorDefaultBase {

  public:
    MaxIterations(const int maxIterations = 100);

    pike::SolveStatus checkStatus(const pike::Solver& solver, const CheckType checkType = pike::COMPLETE);
    
    pike::SolveStatus getStatus() const;
    
    void reset();

    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=verbLevel_default) const;

    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  private:
    int maximumIterations_;
    int currentIterations_;
    Teuchos::RCP<Teuchos::ParameterList> validParameters_;
    pike::SolveStatus status_;
  };

}

#endif
