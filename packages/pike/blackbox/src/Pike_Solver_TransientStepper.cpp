#include "Pike_Solver_TransientStepper.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_Comm.hpp"
#include "Pike_BlackBoxModelEvaluator.hpp"

namespace pike {

  TransientStepper::TransientStepper() :
    currentTimeStep_(0),
    maxTimeSteps_(-1),
    beginTime_(-1.0),
    endTime_(-1.0),
    initialStepSize_(-1.0),
    minStepSize_(-1.0),
    maxStepSize_(-1.0),
    stepGrowthFactor_(2.0),
    stepDecreaseFactor_(0.5)
  {
    validParameters_ = Teuchos::parameterList("pike::TransientStepper::validParameters");

    this->getNonconstValidParameters()->set("Type","Transient Stepper");
    this->getNonconstValidParameters()->set("Internal Solver Sublist","");
    this->getNonconstValidParameters()->set("Maximum Number of Time Steps",-1);
    this->getNonconstValidParameters()->set("Begin Time",-1.0);
    this->getNonconstValidParameters()->set("End Time",-1.0);
    this->getNonconstValidParameters()->set("Initial Time Step Size",-1.0);
    this->getNonconstValidParameters()->set("Minimum Time Step Size",-1.0);
    this->getNonconstValidParameters()->set("Maximum Time Step Size",-1.0);
    this->getNonconstValidParameters()->set("Time Step Size Growth Factor", 2.0);
    this->getNonconstValidParameters()->set("Time Step Size Decrease Factor",0.5);
    this->getNonconstValidParameters()->set("Check Points",checkPoints_,"Time step size will be adjusted to hit these points exactly.");

    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  TransientStepper::~TransientStepper() {}

  void TransientStepper::setSolver(const Teuchos::RCP<pike::Solver>& solver)
  {
    solver_ = solver;
  }

  void TransientStepper::registerComm(const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    TEUCHOS_ASSERT(nonnull(solver_));
    solver_->registerComm(comm);
  }

  void TransientStepper::registerModelEvaluator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me)
  {
    TEUCHOS_ASSERT(nonnull(solver_));
    solver_->registerModelEvaluator(me);
  }

  void TransientStepper::registerDataTransfer(const Teuchos::RCP<pike::DataTransfer>& dt)
  {
    TEUCHOS_ASSERT(nonnull(solver_));
    solver_->registerDataTransfer(dt);
  }

  void TransientStepper::completeRegistration()
  {
    // Set the defaults so the user doesn't have to set the parameter list
    if (is_null(this->getMyParamList())) {
      Teuchos::RCP<Teuchos::ParameterList> defaultParameters = Teuchos::parameterList();
      this->setParameterList(defaultParameters);
    }

    solver_->completeRegistration();

    registrationComplete_ = true;
  }

  Teuchos::RCP<const pike::BlackBoxModelEvaluator> 
  TransientStepper::getModelEvaluator(const std::string& in_name) const
  { return solver_->getModelEvaluator(in_name); }

  const std::vector<Teuchos::RCP<const pike::BlackBoxModelEvaluator> > 
  TransientStepper::getModelEvaluators() const
  { return solver_->getModelEvaluators(); }

  Teuchos::RCP<const pike::DataTransfer> TransientStepper::getDataTransfer(const std::string& in_name) const
  { return solver_->getDataTransfer(in_name); }

  const std::vector<Teuchos::RCP<const pike::DataTransfer> > 
  TransientStepper::getDataTransfers() const
  {return solver_->getDataTransfers(); }

  void TransientStepper::initialize()
  { solver_->initialize(); }

  pike::SolveStatus TransientStepper::step()
  { 
    // NOT implemented yet!
    TEUCHOS_ASSERT(false);
  }

  pike::SolveStatus TransientStepper::solve()
  {
    // Not implemented yet!
    TEUCHOS_ASSERT(false);
  }

  void TransientStepper::finalize()
  { solver_->finalize(); }
  
  void TransientStepper::reset()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
			       "Error in pike::TransientStepper::reset() - The reset() method has been called on a Transient Solver.  This is currently not allowed!");
  }

  pike::SolveStatus TransientStepper::getStatus() const
  { return solver_->getStatus(); }

  int TransientStepper::getNumberOfIterations() const
  { return solver_->getNumberOfIterations(); }

  void TransientStepper::addObserver(const Teuchos::RCP<pike::SolverObserver>& observer)
  { solver_->addObserver(observer); }

  std::vector<Teuchos::RCP<pike::SolverObserver> > TransientStepper::getObservers() const
  { return solver_->getObservers(); }

  void TransientStepper::setStatusTests(const Teuchos::RCP<pike::StatusTest>& statusTests)
  { solver_->setStatusTests(statusTests); }

  Teuchos::RCP<const pike::StatusTest> 
  TransientStepper::getStatusTests() const
  {return solver_->getStatusTests(); }

  std::string TransientStepper::name() const
  { return solver_->name(); }

  void TransientStepper::setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList)
  {
    paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
    
    maxTimeSteps_ = paramList->get<double>("Maximum Number of Time Steps");
    beginTime_ = paramList->get<double>("Begin Time");
    endTime_ = paramList->get<double>("End Time",-1.0);
    initialStepSize_ = paramList->get<double>("Initial Time Step Size",-1.0);
    minStepSize_ = paramList->get<double>("Minimum Time Step Size",-1.0);
    maxStepSize_ = paramList->get<double>("Maximum Time Step Size",-1.0);
    stepGrowthFactor_ = paramList->get<double>("Time Step Size Growth Factor", 2.0);
    stepDecreaseFactor_ = paramList->get<double>("Time Step Size Decrease Factor",0.5);
    checkPoints_ = paramList->get<Teuchos::Array<double> >("Check Points");

    TEUCHOS_ASSERT(beginTime_ < endTime_);
    TEUCHOS_ASSERT(minStepSize_ < maxStepSize_);
    TEUCHOS_ASSERT(initialStepSize_ >= minStepSize_);
    TEUCHOS_ASSERT(initialStepSize_ <= maxStepSize_);
    TEUCHOS_ASSERT(stepGrowthFactor_ <= 1.0);
    TEUCHOS_ASSERT(stepDecreaseFactor_ <= 1.0);

    std::sort(checkPoints_.begin(),checkPoints_.end());
    for (Teuchos::Array<double>::size_type i=0; i < checkPoints_.size(); ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(checkPoints_[i] < beginTime_,
				 std::runtime_error,
				 "");
      TEUCHOS_TEST_FOR_EXCEPTION(checkPoints_[i] > endTime_,
				 std::runtime_error,
				 "");
    }

    TEUCHOS_TEST_FOR_EXCEPTION(paramList->get<std::string>("Internal Solver Sublist") == "",
			       std::runtime_error,
			       "Error in pike::TransientStepper::setParameterList(): The \"Internal Solver Sublist\" must be set to a valid sublist!");

    this->setMyParamList(paramList);
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  TransientStepper::getValidParameters() const
  { return validParameters_; }

  Teuchos::RCP<Teuchos::ParameterList>
  TransientStepper::getNonconstValidParameters()
  { return validParameters_; }


}
