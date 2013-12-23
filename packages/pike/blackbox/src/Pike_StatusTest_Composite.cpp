#include "Pike_StatusTest_Composite.hpp"
#include "Pike_Solver.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include <cmath>

namespace pike {

  Composite::Composite(const pike::Composite::CompositeType type) :
    status_(pike::UNCHECKED)
  {
    validParameters_ = Teuchos::parameterList("Valid Parameters: Composite");
    Teuchos::setStringToIntegralParameter<int>(
        "Type",
        "AND",
        "Determines the form of the DCO_M model in the Momentum equation",
        Teuchos::tuple<std::string>("AND", "OR"),
        validParameters_.get()
        );
    validParameters_->disableRecursiveValidation();
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  pike::SolveStatus Composite::checkStatus(const pike::Solver& solver, const CheckType checkType)
  {
    if (type_ == pike::Composite::AND)
      this->checkAnd(solver,checkType);
    else
      this->checkOr(solver,checkType);
    
    return status_;
  }
  
  void Composite::checkAnd(const pike::Solver& solver, const CheckType checkType)
  {
    if (checkType == pike::NONE)
      status_ = pike::UNCHECKED;

    //for (

  }
  
  void Composite::checkOr(const pike::Solver& solver, const CheckType checkType)
  {
    if (checkType == pike::NONE)
      status_ = pike::UNCHECKED;

    //for (

  }

  pike::SolveStatus Composite::getStatus() const
  { return status_; }
  
  void Composite::reset()
  { status_ = pike::UNCHECKED; }
  
  void Composite::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setMyParamList(paramList);
    //maximumIterations_ = paramList->get<int>("Maximum Iterations");
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> Composite::getValidParameters() const
  {
    return validParameters_;
  }
}
