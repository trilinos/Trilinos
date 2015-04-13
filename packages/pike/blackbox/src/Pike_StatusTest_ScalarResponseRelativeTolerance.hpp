#ifndef PIKE_STATUS_TESTS_SCALAR_RESPONSE_RELATIVE_TOLERANCE_HPP
#define PIKE_STATUS_TESTS_SCALAR_RESPONSE_RELATIVE_TOLERANCE_HPP

#include "Pike_StatusTest.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include <iostream>

namespace pike {

  class BlackBoxModelEvaluator;

  class ScalarResponseRelativeTolerance : 
    public pike::StatusTest,
    public Teuchos::ParameterListAcceptorDefaultBase {

  public:
    ScalarResponseRelativeTolerance();

    pike::SolveStatus checkStatus(const pike::Solver& solver, const CheckType checkType = pike::COMPLETE);
    
    pike::SolveStatus getStatus() const;
    
    void reset();
    
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=verbLevel_default) const;

    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  private:
    std::string applicationName_;
    std::string responseName_;
    int responseIndex_;
    double tolerance_;
    int previousIteration_;
    double previousValue_;
    int currentIteration_;
    double currentValue_;
    pike::SolveStatus status_;
    Teuchos::RCP<Teuchos::ParameterList> validParameters_;
    Teuchos::RCP<const pike::BlackBoxModelEvaluator> application_;
  };

}

#endif
