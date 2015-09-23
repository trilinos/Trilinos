#include "Pike_StatusTest_LocalModelFailure.hpp"
#include "Pike_Solver.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Assert.hpp"

namespace pike {

  LocalModelFailure::LocalModelFailure() :
    applicationName_("???"),
    status_(pike::UNCHECKED)
  {
    validParameters_ = Teuchos::parameterList("Valid Parameters: LocalModelFailure");
    validParameters_->set("Type","Local Model Failure","Type of object to build.");
    validParameters_->set("Model Name","???","Name of the BlackBoxModelEvaluator that should be checked.");
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  pike::SolveStatus LocalModelFailure::checkStatus(const pike::Solver& solver, const CheckType checkType)
  {
    if (is_null(application_)) {
      application_ = solver.getModelEvaluator(applicationName_);
      TEUCHOS_ASSERT(nonnull(application_));
    }

    if (solver.getNumberOfIterations() == 0) {
      status_ = pike::UNCHECKED;
      return status_; 
    }

    // Triggers a failure if the local model fails to converge.
    isLocallyConverged_ = application_->isLocallyConverged();
    status_ = isLocallyConverged_ ? pike::UNCONVERGED : pike::FAILED;

    return status_;
  }

  pike::SolveStatus LocalModelFailure::getStatus() const
  { return status_; }
  
  void LocalModelFailure::reset()
  {
    status_ = pike::UNCHECKED;
  }
  
  void LocalModelFailure::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
  {
    out << pike::statusToString(status_)
	<< "Local Model Failure for \"" << applicationName_ << "\"."
	<< std::endl;
  }

  void LocalModelFailure::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setMyParamList(paramList);
    TEUCHOS_ASSERT(paramList->get<std::string>("Type") == "Local Model Failure");
    applicationName_ = paramList->get<std::string>("Model Name");
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> LocalModelFailure::getValidParameters() const
  {
    return validParameters_;
  }
}
