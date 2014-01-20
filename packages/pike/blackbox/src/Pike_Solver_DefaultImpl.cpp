#include "Pike_Solver_DefaultImpl.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_DataTransfer.hpp"
#include "Pike_Observer.hpp"
#include "Teuchos_Assert.hpp"

namespace pike {

  SolverDefaultImpl::SolverDefaultImpl() :
    numberOfIterations_(0),
    status_(pike::UNCHECKED),
    registrationComplete_(false)
  {
    validParameters_ = Teuchos::parameterList();
    validParameters_->set("Print Begin Solve Status",true, "If set to true the status tests will print current status at the beginning of the solve.");
    validParameters_->set("Print Step Status",true, "If set to true the status tests will print current status at the end of each step.");
    validParameters_->set("Print End Solve Status",true,"If set to true the status tests will print current status at the end of the solve.");
    validParameters_->set("Name","pike::Solver","A unique identifier chosen by the user for this solver. Used mainly for distinguishing nodes in a hierarchic problem.");
    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }

  SolverDefaultImpl::~SolverDefaultImpl() {}

  void SolverDefaultImpl::registerModelEvaluator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(registrationComplete_,
			       std::logic_error,
			       "Can NOT register model evaluators after registrationComplete() has been called!");
    models_.push_back(me);
  }
  
  void SolverDefaultImpl::registerDataTransfer(const Teuchos::RCP<pike::DataTransfer>& dt)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(registrationComplete_,
			       std::logic_error,
			       "Can NOT register data transfers after registrationComplete() has been called!");
    transfers_.push_back(dt);
  }
  
  void SolverDefaultImpl::completeRegistration()
  {
    // Set the defaults so the user doesn't have to set the parameter list
    if (is_null(this->getMyParamList())) {
      Teuchos::RCP<Teuchos::ParameterList> defaultParameters = Teuchos::parameterList();
      this->setParameterList(defaultParameters);
    }

    registrationComplete_ = true;
  }

  Teuchos::RCP<const pike::BlackBoxModelEvaluator> SolverDefaultImpl::getModelEvaluator(const std::string& name) const
  {
    for (ModelConstIterator m = models_.begin(); m != models_.end(); ++m)
      if ((*m)->name() == name)
	return *m;

    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Failed to find the ModelEvaluator named \"" << name << "\" in the solver.");
    return Teuchos::null;
  }

  const std::vector<Teuchos::RCP<const pike::BlackBoxModelEvaluator> > SolverDefaultImpl::getModelEvaluators() const
  {
    std::vector<Teuchos::RCP<const pike::BlackBoxModelEvaluator> > constModels;
    std::copy(models_.begin(),models_.end(),constModels.begin());
    return constModels;
  }

  Teuchos::RCP<const pike::DataTransfer> SolverDefaultImpl::getDataTransfer(const std::string& name) const
  {
    for (TransferConstIterator t = transfers_.begin(); t != transfers_.end(); ++t)
      if ((*t)->name() == name)
	return *t;

    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Failed to find the DataTransfer named \"" << name << "\" in the solver.");
    return Teuchos::null;
  }

  const std::vector<Teuchos::RCP<const pike::DataTransfer> > SolverDefaultImpl::getDataTransfers() const
  { std::vector<Teuchos::RCP<const pike::DataTransfer> > constTransfers;
    std::copy(transfers_.begin(),transfers_.end(),constTransfers.begin());
    return constTransfers;
  }

  pike::SolveStatus SolverDefaultImpl::step()
  {
    for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
      (*observer)->observeBeginStep(*this);

    this->stepImplementation();
    ++numberOfIterations_;

    status_ = statusTests_->checkStatus(*this);
    
    if (printStepStatus_) {
      Teuchos::RCP<Teuchos::FancyOStream> os = this->getOStream();
      *os << "\n** Step " << this->getNumberOfIterations() << " Status **" << std::endl;
      os->pushTab(defaultIndentation);
      *os << *statusTests_;
      os->popTab();
    }

    for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
      (*observer)->observeEndStep(*this);

    return status_;
  }
  
  pike::SolveStatus SolverDefaultImpl::solve()
  {
    TEUCHOS_ASSERT(registrationComplete_);

    for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
      (*observer)->observeBeginSolve(*this);
    
    if (printBeginSolveStatus_) {
      Teuchos::RCP<Teuchos::FancyOStream> os = this->getOStream();
      *os << "\n** Begin Solve Status **" << std::endl;
      os->pushTab(defaultIndentation);
      *os << *statusTests_;
      os->popTab();
    }

    status_ = statusTests_->checkStatus(*this);

    while ( (status_ != CONVERGED) && (status_ != FAILED) )
      this->step();
    
    if (printEndSolveStatus_) {
      Teuchos::RCP<Teuchos::FancyOStream> os = this->getOStream();
      *os << "\n** End Solve Status **" << std::endl;
      os->pushTab(defaultIndentation);
      *os << *statusTests_;
      os->popTab();
    }

    for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
      (*observer)->observeEndSolve(*this);

    if (status_ == CONVERGED)
      for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
	(*observer)->observeConvergedSolve(*this);

    if (status_ == FAILED)
      for (ObserverIterator observer = observers_.begin(); observer != observers_.end(); ++observer)
	(*observer)->observeFailedSolve(*this);

    return status_;
  }

  void SolverDefaultImpl::reset()
  {
    numberOfIterations_ = 0;
    status_ = pike::UNCHECKED;
  }
  
  pike::SolveStatus SolverDefaultImpl::getStatus() const
  { return status_; }

  int SolverDefaultImpl::getNumberOfIterations() const
  { return numberOfIterations_; }
  
  void SolverDefaultImpl::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    printBeginSolveStatus_ = paramList->get<bool>("Print Begin Solve Status");
    printStepStatus_ = paramList->get<bool>("Print Step Status");
    printEndSolveStatus_ = paramList->get<bool>("Print End Solve Status");
    name_ = paramList->get<std::string>("Name");
    this->setMyParamList(paramList);
  }
  
  Teuchos::RCP<const Teuchos::ParameterList> SolverDefaultImpl::getValidParameters() const
  { return validParameters_; }

  Teuchos::RCP<Teuchos::ParameterList> SolverDefaultImpl::getNonconstValidParameters()
  { return validParameters_; }

  void SolverDefaultImpl::addObserver(const Teuchos::RCP<pike::Observer>& observer)
  {
    observers_.push_back(observer);
  }

  std::vector<Teuchos::RCP<pike::Observer> > SolverDefaultImpl::getObservers() const
  {
    return observers_;
  }

  void SolverDefaultImpl::setStatusTests(const Teuchos::RCP<pike::StatusTest>& statusTests)
  {
    statusTests_ = statusTests;
  }

  Teuchos::RCP<const pike::StatusTest> SolverDefaultImpl::getStatusTests() const
  {
    return statusTests_;
  }
  
  std::string SolverDefaultImpl::name() const
  {
    return name_;
  }
}
