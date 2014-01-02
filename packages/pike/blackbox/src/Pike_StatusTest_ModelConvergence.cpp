#include "Pike_StatusTest_ModelConvergence.hpp"
#include "Pike_Solver.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Assert.hpp"

namespace pike {

  ModelConvergence::ModelConvergence() :
    applicationName_("???"),
    convergenceType_(pike::ModelConvergence::Local),
    status_(pike::UNCHECKED)
  {
    validParameters_ = Teuchos::parameterList("Valid Parameters: ModelConvergence");
    validParameters_->set("Type","Model Convergence","Type of object to build.");
    validParameters_->set("Model Name","???","Name of the BlackBoxModelEvaluator that should be checked.");
    Teuchos::setStringToIntegralParameter<int>(
        "Check Type",
        "Local",
        "Determines whether to check the local convergence (status of last call to solve()) for this model or to check the global convergence of the coupled problem as assessed by this model.",
        Teuchos::tuple<std::string>("Local", "Global"),
        validParameters_.get()
        );
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  pike::SolveStatus ModelConvergence::checkStatus(const pike::Solver& solver, const CheckType checkType)
  {
    if (is_null(application_)) {
      application_ = solver.getModelEvaluator(applicationName_);
      TEUCHOS_ASSERT(nonnull(application_));
    }

    if (solver.getNumberOfIterations() == 0) {
      status_ = pike::UNCHECKED;
      return status_; 
    }

    if (convergenceType_ == pike::ModelConvergence::Local) {
      // Triggers a failure if the local model fails to converge.
      isLocallyConverged_ = application_->isConverged();
      status_ = isLocallyConverged_ ? pike::UNCONVERGED : pike::FAILED;
    }
    else { // Global
      // Triggers convergence if the model is satisfied with global
      // convergence of the coupled problem.
      isGloballyConverged_ = application_->isGloballyConverged();
      status_ = isGloballyConverged_ ? pike::CONVERGED : pike::UNCONVERGED;
    }

    return status_;
  }

  pike::SolveStatus ModelConvergence::getStatus() const
  { return status_; }
  
  void ModelConvergence::reset()
  {
    status_ = pike::UNCHECKED;
  }
  
  void ModelConvergence::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
  {
    out << pike::statusToString(status_)
	<< "Model Convergence for \"" << applicationName_ << "\": ";
    if (convergenceType_ == Local)
      out << "Local = " << isLocallyConverged_ << std::endl;
    else
      out << "Global = " << isGloballyConverged_ << std::endl;
  }

  void ModelConvergence::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setMyParamList(paramList);
    TEUCHOS_ASSERT(paramList->get<std::string>("Type") == "Model Convergence");
    applicationName_ = paramList->get<std::string>("Model Name");
    std::string checkTypeName = paramList->get<std::string>("Check Type");

    if (checkTypeName == "Local")
      convergenceType_ = pike::ModelConvergence::Local;
    else if (checkTypeName == "Global")
      convergenceType_ = pike::ModelConvergence::Global;
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> ModelConvergence::getValidParameters() const
  {
    return validParameters_;
  }
}
