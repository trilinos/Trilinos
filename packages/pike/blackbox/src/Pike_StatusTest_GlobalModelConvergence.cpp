#include "Pike_StatusTest_GlobalModelConvergence.hpp"
#include "Pike_Solver.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Assert.hpp"

namespace pike {

  GlobalModelConvergence::GlobalModelConvergence() :
    applicationName_("???"),
    status_(pike::UNCHECKED)
  {
    validParameters_ = Teuchos::parameterList("Valid Parameters: GlobalModelConvergence");
    validParameters_->set("Type","Global Model Convergence","Type of object to build.");
    validParameters_->set("Model Name","???","Name of the BlackBoxModelEvaluator that should be checked.");
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  pike::SolveStatus GlobalModelConvergence::checkStatus(const pike::Solver& solver, const CheckType checkType)
  {
    if (is_null(application_)) {
      application_ = solver.getModelEvaluator(applicationName_);
      TEUCHOS_ASSERT(nonnull(application_));
    }

    if (solver.getNumberOfIterations() == 0) {
      status_ = pike::UNCHECKED;
      return status_; 
    }

    // Triggers convergence if the model is satisfied with global
    // convergence of the coupled problem.
    isGloballyConverged_ = application_->isGloballyConverged();
    status_ = isGloballyConverged_ ? pike::CONVERGED : pike::UNCONVERGED;

    return status_;
  }

  pike::SolveStatus GlobalModelConvergence::getStatus() const
  { return status_; }
  
  void GlobalModelConvergence::reset()
  {
    status_ = pike::UNCHECKED;
  }
  
  void GlobalModelConvergence::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
  {
    out << pike::statusToString(status_)
	<< "Global Model Convergence for \"" << applicationName_ << "\"." 
	<< std::endl;
  }

  void GlobalModelConvergence::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setMyParamList(paramList);
    TEUCHOS_ASSERT(paramList->get<std::string>("Type") == "Global Model Convergence");
    applicationName_ = paramList->get<std::string>("Model Name");
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> GlobalModelConvergence::getValidParameters() const
  {
    return validParameters_;
  }
}
