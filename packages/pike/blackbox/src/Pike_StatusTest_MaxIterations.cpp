#include "Pike_StatusTest_MaxIterations.hpp"
#include "Pike_Solver.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include <cmath>

namespace pike {

  MaxIterations::MaxIterations(const int maxIterations) :
    maximumIterations_(maxIterations),
    currentIterations_(-1),
    status_(pike::UNCHECKED)
  {
    validParameters_ = Teuchos::parameterList("Valid Parameters: MaxIterations");
    validParameters_->set("Type","Maximum Iterations","Type of test.");
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
  
  void MaxIterations::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel) const
  {
    out << pike::statusToString(status_) 
	<< "Num Iterations = " << currentIterations_ 
	<< ", limited to " << maximumIterations_ << std::endl;
  }

  void MaxIterations::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setMyParamList(paramList);
    TEUCHOS_ASSERT(paramList->get<std::string>("Type") == "Maximum Iterations");
    maximumIterations_ = paramList->get<int>("Maximum Iterations");
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> MaxIterations::getValidParameters() const
  {
    return validParameters_;
  }
}
