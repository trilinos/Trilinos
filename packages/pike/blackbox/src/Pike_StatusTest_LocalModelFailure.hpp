#ifndef PIKE_STATUS_TESTS_LOCAL_MODEL_FAILURE_HPP
#define PIKE_STATUS_TESTS_LOCAL_MODEL_FAILURE_HPP

#include "Pike_StatusTest.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include <iostream>

namespace pike {

  class BlackBoxModelEvaluator;

  /** \brief Failure test for local convergence of a model.
   */
  class LocalModelFailure : 
    public pike::StatusTest,
    public Teuchos::ParameterListAcceptorDefaultBase {

  public:

    LocalModelFailure();

    pike::SolveStatus checkStatus(const pike::Solver& solver, const CheckType checkType = pike::COMPLETE);
    
    pike::SolveStatus getStatus() const;
    
    void reset();
    
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=verbLevel_default) const;

    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  private:
    std::string applicationName_;
    bool isLocallyConverged_;
    pike::SolveStatus status_;
    Teuchos::RCP<Teuchos::ParameterList> validParameters_;
    Teuchos::RCP<const pike::BlackBoxModelEvaluator> application_;
  };

}

#endif
