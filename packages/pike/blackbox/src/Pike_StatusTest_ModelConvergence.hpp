#ifndef PIKE_STATUS_TESTS_MODEL_CONVERGENCE_HPP
#define PIKE_STATUS_TESTS_MODEL_CONVERGENCE_HPP

#include "Pike_StatusTest.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include <iostream>

namespace pike {

  class BlackBoxModelEvaluator;

  /** \brief Failure test for local convergence of a model or convergence test for global convergence of a model.
   */
  class ModelConvergence : 
    public pike::StatusTest,
    public Teuchos::ParameterListAcceptorDefaultBase {

  public:

    enum ConvergenceType {
      Local,
      Global
    };

    ModelConvergence();

    pike::SolveStatus checkStatus(const pike::Solver& solver, const CheckType checkType = pike::COMPLETE);
    
    pike::SolveStatus getStatus() const;
    
    void reset();
    
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=verbLevel_default) const;

    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  private:
    std::string applicationName_;
    ConvergenceType convergenceType_;
    bool isLocallyConverged_;
    bool isGloballyConverged_;
    pike::SolveStatus status_;
    Teuchos::RCP<Teuchos::ParameterList> validParameters_;
    Teuchos::RCP<const pike::BlackBoxModelEvaluator> application_;
  };

}

#endif
