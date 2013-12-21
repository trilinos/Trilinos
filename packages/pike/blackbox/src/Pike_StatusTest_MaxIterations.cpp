#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_Solver.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include <cmath>

namespace pike {

  MaxIterations::MaxIterations() :
    maximumIterations_(-1),
    currentIterations_(-1),
    status_(pike::UNCHECKED)
  {
    validParameters_ = Teuchos::parameterList("Valid Parameters: MaxIterations");
    validParameters_->set("Maximum Iterations",-1,"Maximum number of iterations before the test returns a failed state.");
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  pike::SolveStatus MaxIterations::checkStatus(const pike::Solver& solver, const CheckType checkType)
  {
    currentIterations_ = solver.getNumberOfIterations();

    if (currentIterations_ < maximumIterations_)
      status_ = pike::UNCONVERGED;
    else
      status_ = pike::FAILED;

    return status_;
  }

  pike::SolveStatus MaxIterations::getStatus() const
  { return status_; }
  
  void MaxIterations::reset()
  { status_ = pike::UNCHECKED; }
  
  void MaxIterations::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setParameterList(paramList);
    maximumIterations_ = paramList->get<int>("Maximum Iterations");
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> MaxIterations::getValidParameters() const
  {
    return validParameters_;
  }
}
