#include "Pike_Solver_DefaultImplementation.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_DataTransfer.hpp"
#include "Pike_Observer.hpp"
#include "Teuchos_Assert.hpp"

namespace pike {

  SolverDefaultImplementation::SolverDefaultImplementation() :
    numberOfIterations_(0),
    status_(pike::UNCHECKED),
    registrationComplete_(false)
  {
    validParameters_ = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }

  SolverDefaultImplementation::~SolverDefaultImplementation() {}

  void SolverDefaultImplementation::registerModelEvaluator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(registrationComplete_,
			       std::logic_error,
			       "Can NOT register model evaluators after registrationComplete() has been called!");
    models_.push_back(me);
  }
  
  void SolverDefaultImplementation::registerDataTransfer(const Teuchos::RCP<pike::DataTransfer>& dt)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(registrationComplete_,
			       std::logic_error,
			       "Can NOT register data transfers after registrationComplete() has been called!");
    transfers_.push_back(dt);
  }
  
  void SolverDefaultImplementation::completeRegistration()
  {
    registrationComplete_ = true;
  }

  Teuchos::RCP<const pike::BlackBoxModelEvaluator> SolverDefaultImplementation::getModelEvaluator(const std::string& name) const
  {
    typedef std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> >::const_iterator ModelIterator;
    for (ModelIterator m = models_.begin(); m != models_.end(); ++m)
      if ((*m)->name() == name)
	return *m;

    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Failed to find the DataTransfer named \"" << name << "\" in the solver.");
    return Teuchos::null;
  }

  Teuchos::RCP<const pike::DataTransfer> SolverDefaultImplementation::getDataTransfer(const std::string& name) const
  {
    typedef std::vector<Teuchos::RCP<pike::DataTransfer> >::const_iterator TransferIterator;
    for (TransferIterator t = transfers_.begin(); t != transfers_.end(); ++t)
      if ((*t)->name() == name)
	return *t;

    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Failed to find the DataTransfer named \"" << name << "\" in the solver.");
    return Teuchos::null;
  }

  pike::SolveStatus SolverDefaultImplementation::step()
  {
    for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
      (*observer)->observeBeginStep(*this);

    this->stepImplementation();
    ++numberOfIterations_;

    status_ = statusTests_->checkStatus(*this);
    
    for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
      (*observer)->observeEndStep(*this);

    if (status_ == CONVERGED)
      for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
	(*observer)->observeConvergedSolve(*this);

    if (status_ == FAILED)
      for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
	(*observer)->observeFailedSolve(*this);

    return status_;
  }
  
  pike::SolveStatus SolverDefaultImplementation::solve()
  {
    TEUCHOS_ASSERT(registrationComplete_);

    for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
      (*observer)->observeBeginSolve(*this);

    status_ = statusTests_->checkStatus(*this);

    while ( (status_ != CONVERGED) || (status_ != FAILED) )
      this->step();

    for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
      (*observer)->observeEndSolve(*this);

    return status_;
  }

  void SolverDefaultImplementation::reset()
  {
    numberOfIterations_ = 0;
    status_ = pike::UNCHECKED;
  }
  
  pike::SolveStatus SolverDefaultImplementation::getStatus() const
  { return status_; }

  int SolverDefaultImplementation::getNumberOfIterations() const
  { return numberOfIterations_; }
  
  void SolverDefaultImplementation::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    this->setParameterList(paramList);
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> SolverDefaultImplementation::getValidParameters() const
  { return validParameters_; }

  Teuchos::RCP<Teuchos::ParameterList> SolverDefaultImplementation::getNonconstValidParameters()
  { return validParameters_; }

  void SolverDefaultImplementation::addObserver(const Teuchos::RCP<pike::Observer>& observer)
  {
    observers_.push_back(observer);
  }

  std::vector<Teuchos::RCP<pike::Observer> > SolverDefaultImplementation::getObservers() const
  {
    return observers_;
  }

  void SolverDefaultImplementation::setStatusTests(const Teuchos::RCP<pike::StatusTest>& statusTests)
  {
    statusTests_ = statusTests;
  }

  Teuchos::RCP<const pike::StatusTest> SolverDefaultImplementation::getStatusTests() const
  {
    return statusTests_;
  }
}
