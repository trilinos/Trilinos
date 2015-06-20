#include "Pike_StatusTest_Composite.hpp"
#include "Pike_Solver.hpp"
#include "Pike_StatusTest_AbstractFactory.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ParameterList.hpp"

#include <cmath>

namespace pike {

  Composite::Composite(const pike::Composite::CompositeType type) :
    type_(type),
    status_(pike::UNCHECKED)
  {
    validParameters_ = Teuchos::parameterList("Valid Parameters: Composite");
    Teuchos::setStringToIntegralParameter<int>(
        "Type",
        "Composite OR",
        "Determines the form of the DCO_M model in the Momentum equation",
        Teuchos::tuple<std::string>("Composite AND", "Composite OR"),
        validParameters_.get()
        );
    validParameters_->disableRecursiveValidation();
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  void Composite::addTest(const Teuchos::RCP<pike::StatusTest>& t)
  {
    tests_.push_back(t);
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

    bool isUnconverged = false;
    pike::CheckType subCheckType = checkType;
    
    for (TestIterator test = tests_.begin(); test != tests_.end(); ++test) {
      pike::SolveStatus testStatus = (*test)->checkStatus(solver,subCheckType);
      
      if (testStatus == pike::UNCONVERGED) {
	isUnconverged = true;
	status_ = pike::UNCONVERGED;
	
	// If any are unconverged, "AND" means all are unconverged, so
	// we can disable the rest of the checks (but still need to
	// call check status in case they need to store solver
	// history.
	if (checkType == pike::MINIMAL)
	  subCheckType = pike::NONE;
      }

      // If we are not unconverged and 
      if ( (!isUnconverged) && (status_ == pike::UNCONVERGED) ) {
	status_ = testStatus;
      }
    }

  }
  
  void Composite::checkOr(const pike::Solver& solver, const CheckType checkType)
  {
    if (checkType == pike::NONE)
      status_ = pike::UNCHECKED;
    else
      status_ = pike::UNCONVERGED;

    pike::CheckType subCheckType = checkType;

    for (TestIterator test = tests_.begin(); test != tests_.end(); ++test) {
      pike::SolveStatus testStatus = (*test)->checkStatus(solver,subCheckType);

      if ( (status_ == pike::UNCONVERGED) && (testStatus != pike::UNCONVERGED) ) {
	status_ = testStatus;

	if (checkType == pike::MINIMAL)
	  subCheckType = pike::NONE;
      }
    }

  }

  pike::SolveStatus Composite::getStatus() const
  { return status_; }
  
  void Composite::reset()
  {    
    for (TestIterator test = tests_.begin(); test != tests_.end(); ++test)
      (*test)->reset();

    status_ = pike::UNCHECKED;
  }
    
  void Composite::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
  {
    out << pike::statusToString(status_);
    if (type_ == AND)
      out << "AND";
    else
      out << "OR";
    
    out << " Composite (" << tests_.size() << " subtests):" << std::endl;

    out.pushTab(defaultIndentation);

    for (TestConstIterator t = tests_.begin(); t != tests_.end(); ++t)
      (*t)->describe(out,verbLevel);

    out.popTab();
  }

  void Composite::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& p,
				   const pike::StatusTestAbstractFactory& factory)
  {
    // Don't call validation.  Sublist names can be arbitrary (or nonexistent)
    //paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    myParameters_ = p;

    // Loop over sublists for different tests
    for (Teuchos::ParameterList::ConstIterator sublistEntry = p->begin();
	 sublistEntry != p->end(); ++sublistEntry) {

      if ( (sublistEntry->first == "Type") && (sublistEntry->second.isType<std::string>()) ) {
	std::string stringType = "";
	if (sublistEntry->second.getValue(&stringType) == "Composite AND")
	  type_ = pike::Composite::AND;
	else if (sublistEntry->second.getValue(&stringType) == "Composite OR")
	  type_ = pike::Composite::OR;
	else {
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
				     "The parameter with key \"" << sublistEntry->first 
				     << "\" and value \"" << sublistEntry->second.getValue(&stringType)
				     << "\" is not valid for Composite object construction!");
	}
      }
      else if (sublistEntry->second.isList()) {
	Teuchos::RCP<Teuchos::ParameterList> sublist = Teuchos::sublist(p,sublistEntry->first,true);
	Teuchos::RCP<StatusTest> subtest = factory.buildStatusTests(sublist);
	this->addTest(subtest);
      }
      else {
	TEUCHOS_TEST_FOR_EXCEPTION(sublistEntry->second.isList(),
				   std::logic_error,
				   "The parameter sublist key \"" << sublistEntry->first << "\" must be a sublist or a string that determines the type of Composite test!"); 
      }
    }
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> Composite::getValidParameters() const
  {
    return validParameters_;
  }

  // nonmember ctor
  Teuchos::RCP<pike::Composite> composite(const pike::Composite::CompositeType type)
  {
    return Teuchos::rcp(new pike::Composite(type));
  }
  
}
