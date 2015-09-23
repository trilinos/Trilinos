#include "Pike_StatusTest_LocalModelConvergence.hpp"
#include "Pike_Solver.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Assert.hpp"

namespace pike {

  LocalModelConvergence::LocalModelConvergence() :
    applicationName_("???"),
    status_(pike::UNCHECKED)
  {
    validParameters_ = Teuchos::parameterList("Valid Parameters: LocalModelConvergence");
    validParameters_->set("Type","Local Model Convergence","Type of object to build.");
    validParameters_->set("Model Name","???","Name of the BlackBoxModelEvaluator that should be checked.");
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  pike::SolveStatus LocalModelConvergence::checkStatus(const pike::Solver& solver, const CheckType checkType)
  {
    if (is_null(application_)) {
      application_ = solver.getModelEvaluator(applicationName_);
      TEUCHOS_ASSERT(nonnull(application_));
    }

    if (solver.getNumberOfIterations() == 0) {
      status_ = pike::UNCHECKED;
      return status_; 
    }

    // Triggers convergence if the local model converged.
    isLocallyConverged_ = application_->isLocallyConverged();
    status_ = isLocallyConverged_ ? pike::CONVERGED : pike::UNCONVERGED;

    return status_;
  }

  pike::SolveStatus LocalModelConvergence::getStatus() const
  { return status_; }
  
  void LocalModelConvergence::reset()
  {
    status_ = pike::UNCHECKED;
  }
  
  void LocalModelConvergence::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
  {
    out << pike::statusToString(status_)
	<< "Local Model Convergence for \"" << applicationName_ << "\"."
	<< std::endl;
  }

  void LocalModelConvergence::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setMyParamList(paramList);
    TEUCHOS_ASSERT(paramList->get<std::string>("Type") == "Local Model Convergence");
    applicationName_ = paramList->get<std::string>("Model Name");
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> LocalModelConvergence::getValidParameters() const
  {
    return validParameters_;
  }
}
