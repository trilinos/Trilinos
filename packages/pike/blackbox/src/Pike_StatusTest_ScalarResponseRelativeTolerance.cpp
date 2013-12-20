#include "Pike_StatusTest_ScalarResponseRelativeTolerance.hpp"
#include "Pike_Solver.hpp"
#include "Pike_Response_Scalar.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include <cmath>

namespace pike {

  ScalarResponseRelativeTolerance::ScalarResponseRelativeTolerance() :
    applicationName_("???"),
    responseName_("???"),
    tolerance_(1.0e-4),
    previousIteration_(-1),
    previousValue_(0.0),
    currentIteration_(-1),
    currentValue_(0.0),
    status_(pike::UNCHECKED)
  {
    validParameters_ = Teuchos::parameterList("Valid Parameters: ScalarResponseRelativeTolerance");
    validParameters_->set("Application Name","???","Name of the BlackBoxModelEvaluator that contains the response");
    validParameters_->set("Response Name","???","Name of the ScalarResponse to check");
    validParameters_->set("Tolerance",1.0e-4,"Relative tolerance required for convergence");
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  pike::SolveStatus ScalarResponseRelativeTolerance::checkStatus(const pike::Solver& solver, const CheckType checkType)
  {
    // since we need previous iteration info, we can't do a NONE or
    // MINIMAL check.

    // Relative test requires difference between successive steps, so
    // first step is always unconverged.
    if (solver.getNumberOfIterations() == 0) {
      if (is_null(application_)) {
	application_ = solver.getModelEvaluator(applicationName_);
	responseIndex_ = application_->getResponseIndex(responseName_);
      }

      previousIteration_ = -1;
      currentIteration_ = 0;
      previousValue_ = 0.0;
      currentValue_ = pike::getScalarResponse<double>(*application_->getResponse(responseIndex_));
      status_ = pike::UNCONVERGED;
      return status_;
    }

    if (solver.getNumberOfIterations() == previousIteration_)
      return status_; // has already been checked this iteration

    // took a step and now need to check
    if (solver.getNumberOfIterations() > previousIteration_) {
      previousIteration_ = currentIteration_;
      currentIteration_ = solver.getNumberOfIterations();
      previousValue_ = currentValue_;
      currentValue_ = pike::getScalarResponse<double>(*application_->getResponse(responseIndex_));
      if (std::abs(currentValue_ - previousValue_) < tolerance_)
	status_ = pike::UNCONVERGED;
      else
	status_ = pike::UNCONVERGED;
    }
    return status_;
  }

  pike::SolveStatus ScalarResponseRelativeTolerance::getStatus() const
  { return status_; }
  
  void ScalarResponseRelativeTolerance::reset()
  {
    previousIteration_ = -1;
    currentIteration_ = 0;
    previousValue_ = 0.0;
    currentValue_ = pike::getScalarResponse<double>(*application_->getResponse(responseIndex_));
    status_ = pike::UNCHECKED;
  }
  
  void ScalarResponseRelativeTolerance::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setParameterList(paramList);
    applicationName_ = paramList->get<std::string>("Application Name");
    responseName_ = paramList->get<std::string>("Response Name");
    tolerance_ = paramList->get<double>("Tolerance");
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> ScalarResponseRelativeTolerance::getValidParameters() const
  {
    return validParameters_;
  }
}
