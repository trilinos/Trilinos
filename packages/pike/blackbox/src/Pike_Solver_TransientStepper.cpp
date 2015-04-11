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
    currentTime_(beginTime_),
    initialStepSize_(-1.0),
    currentStepSize_(initialStepSize_),
    minStepSize_(-1.0),
    maxStepSize_(-1.0),
    stepGrowthFactor_(2.0),
    stepDecreaseFactor_(0.5),
    printTimeStepSummary_(true),
    printTimeStepDetails_(true),
    numConvergedTimeStepsBeforeGrowth_(3),
    overallStatus_(UNCHECKED),
    timeStepStatus_(UNCHECKED),
    totalNumFailedSteps_(0),
    numConsecutiveFailedTimeSteps_(0),
    numConsecutiveConvergedTimeSteps_(0),
    registrationComplete_(false)
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
    this->getNonconstValidParameters()->set("Check Points",checkPointsVec_,"Time step size will be adjusted to hit these points exactly.");
    this->getNonconstValidParameters()->set("Print Time Step Summary",true,"Prints time step summary information to ostream.");
    this->getNonconstValidParameters()->set("Print Time Step Details",true,"Prints details of time step to ostream.");
    this->getNonconstValidParameters()->set("Number Converged Time Steps for Growth",3,"Delays growing a time step size towards the maximum until a specified number of consecutive time steps have converged.  This helps prevent oscillation between cutting and increasing on alternate steps.");

    Teuchos::setupVerboseObjectSublist(validParameters_.get());
  }
  
  TransientStepper::~TransientStepper() {}

  void TransientStepper::setSolver(const Teuchos::RCP<pike::Solver>& solver)
  {
    solver_ = solver;
  }

  void TransientStepper::registerComm(const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    TEUCHOS_ASSERT(!registrationComplete_);
    TEUCHOS_ASSERT(nonnull(solver_));
    solver_->registerComm(comm);
  }

  void TransientStepper::registerModelEvaluator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& me)
  {
    TEUCHOS_ASSERT(!registrationComplete_);
    TEUCHOS_ASSERT(nonnull(solver_));
    solver_->registerModelEvaluator(me);

    // Store off the transient model evaluators.  We do suppor the
    // case where there can be a mix of steady-state and transient
    // model evaluators, so we only consider the transient ones for
    // determinig time step and accepting a completed time step solve.
    if (me->isTransient())
      transientModels_.push_back(me);
  }

  void TransientStepper::registerDataTransfer(const Teuchos::RCP<pike::DataTransfer>& dt)
  {
    TEUCHOS_ASSERT(!registrationComplete_);
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
  {
    solver_->initialize();
  }

  pike::SolveStatus TransientStepper::step()
  { 
    TEUCHOS_ASSERT(registrationComplete_);
    Teuchos::RCP<Teuchos::FancyOStream> rcpOs = this->getOStream();
    Teuchos::FancyOStream& os = *rcpOs;

    ++currentTimeStep_;

    if (printTimeStepSummary_) {
      os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      os << "Beginning time step " << currentTimeStep_ << ": time=" << currentTime_ << std::endl;
    }

    // Loop over a time step until it either converges or fails from
    // max number of time steps or minimum time step.
    bool hitMinTimeStep = false;
    bool achievedFinalTime = false;
    numConsecutiveFailedTimeSteps_ = 0;
    double previousStepSize = currentStepSize_;
    timeStepStatus_ = UNCONVERGED;

    while ( (timeStepStatus_ != CONVERGED) && (timeStepStatus_ != FAILED) ) {

      previousStepSize = currentStepSize_;

      if (printTimeStepDetails_)
	os << "\n  previous step size=" << previousStepSize << std::endl;

      // Increase or decrease current time step size, cap for max step
      // size
      if (currentTimeStep_ > 0) {
	if (numConsecutiveFailedTimeSteps_ > 0) {
	  currentStepSize_ *= stepDecreaseFactor_;
	  
	  if (printTimeStepDetails_)
	    os << "  Decreasing time step:"
	       << "\n    numConsecutiveFailedTimeSteps_ = " << numConsecutiveFailedTimeSteps_
	       << "\n    stepDecreaseFactor_ = " << stepDecreaseFactor_
	       << "\n    new step size = " << currentStepSize_ << std::endl;
	}
	else if (numConsecutiveConvergedTimeSteps_ >= numConvergedTimeStepsBeforeGrowth_) {
	  currentStepSize_ *= stepGrowthFactor_;
	  currentStepSize_ = std::min(currentStepSize_, maxStepSize_);

	  if (printTimeStepDetails_)
	    os << "  Increasing time step:"
	       << "\n    numConsecutiveConvergedTimeSteps_ = " << numConsecutiveConvergedTimeSteps_ 
	       << "\n    numConvergedTimeStepsBeforeGrowth_ = " << numConvergedTimeStepsBeforeGrowth_
	       << "\n    stepGrowthFactor_ = " << stepGrowthFactor_
	       << "\n    new step size = " << currentStepSize_ << std::endl;
	}
      }

      // Limit time step based on application numerical method
      // requirements
      double modelRequestedStepSize = maxStepSize_;
      for (std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> >::const_iterator m = transientModels_.begin();
	   m != transientModels_.end(); ++m) {
	double modelMaxStepSize = (*m)->getMaxTimeStepSize();
	double modelDesiredStepSize = (*m)->getDesiredTimeStepSize();

	// Do not allow simulations that violate application max step
	// size, but allow for violations of desired step size.
	TEUCHOS_TEST_FOR_EXCEPTION(modelMaxStepSize < minStepSize_, std::runtime_error,
				   "Error the application \"" << (*m)->name() 
				   << "\" requested a max step size of " 
				   << modelMaxStepSize 
				   << " but this is less than the minimum step size of " 
				   << minStepSize_ 
				   << " requested by the user.  Terminating run!");

	modelRequestedStepSize = std::min(modelRequestedStepSize, modelMaxStepSize);
	modelRequestedStepSize = std::min(modelRequestedStepSize, modelDesiredStepSize);
      }
      currentStepSize_ = std::min(currentStepSize_, modelRequestedStepSize);

      if (printTimeStepDetails_) {
	os << "  Step size adjustments for applications:" << std::endl;
	for (std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> >::const_iterator m = transientModels_.begin();
	     m != transientModels_.end(); ++m) {
	  os << "    model \"" << (*m)->name() << "\" max step size = " << (*m)->getMaxTimeStepSize() << std::endl;
	  os << "    model \"" << (*m)->name() << "\" desired step size = " << (*m)->getDesiredTimeStepSize() << std::endl;
	}
	os << "    step size after model adjustments = " << currentStepSize_ << std::endl;
      }

      // Clip against minimum step size (this could increase over
      // application desired step size, but still enforces application
      // max step size due to code above)
      if (currentStepSize_ <= minStepSize_) {
	currentStepSize_  = minStepSize_;
	hitMinTimeStep =  true;
      }

      if (printTimeStepDetails_) {
	os << "  Step size adjustments for minStepSize_ = " << minStepSize_ << std::endl;
	os << "    step size after min setp size adjustment = " << currentStepSize_ << std::endl;
      }

      // Cut the time step to hit checkpoints and end solve time
      // exactly (here we are allowed to go smaller than the minimum
      // time step requested by user so this block must occur after
      // the minimum step size adjustment).
      bool hitCheckPoint = false;
      {
	double nextTime = currentTime_ + currentStepSize_;

	if (checkPoints_.size() > 0) {
	  if (nextTime > *checkPoints_.begin()) {
	    currentStepSize_ = *checkPoints_.begin() - currentTime_;
	    hitCheckPoint = true;
	  }

	  if (printTimeStepDetails_) {
	    os << "  Step size adjustments for checkPoints_ = " << *checkPoints_.begin() << std::endl;
	    os << "    step size after checkPoints_ adjustment = " << currentStepSize_ << std::endl;
	  }

	}

	if (nextTime > endTime_) {
	  currentStepSize_ = endTime_ - currentTime_;
	  achievedFinalTime = true;

	  if (printTimeStepDetails_) {
	    os << "  Step size adjustments for endTime_ = " << endTime_ << std::endl;
	    os << "    step size after endTime_ adjustment = " << currentStepSize_ << std::endl;
	  }

	}  
	// Must have this in case we tentatively hit final time, but
	// the step fails and is cut back in this loop.
	else 
	  achievedFinalTime = false;
      }
      
      // Log whether we are at the minimum step size, so that if the
      // solve fails, we know to exit
      if (currentStepSize_ <= minStepSize_)
	hitMinTimeStep =  true;

      // Set the time steps
      for (std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> >::const_iterator m = transientModels_.begin();
	   m != transientModels_.end(); ++m)
	(*m)->setNextTimeStepSize(currentStepSize_);
      
      // Reset the solver status tests for a new solve
      solver_->reset();
  
      // solve step
      if (printTimeStepSummary_) {
	os << "\n  Starting inner solve: step size=" << currentStepSize_ 
	   << ", target time=" << currentTime_+currentStepSize_ << std::endl;
      }
      os.pushTab(defaultIndentation);
      pike::SolveStatus innerSolverStatus = solver_->solve();
      os.popTab();

      // Check for time step status change
      if (innerSolverStatus == CONVERGED) {
	timeStepStatus_ = CONVERGED;

	//make sure to pop checkpoint
	if (hitCheckPoint)
	  checkPoints_.pop_front();

	// log converged and failed time steps
	numConsecutiveFailedTimeSteps_ = 0;
	++numConsecutiveConvergedTimeSteps_;
	
	// time
	currentTime_ += currentStepSize_;

	// Accept the time step in the transient model evaluators
	for (std::vector<Teuchos::RCP<pike::BlackBoxModelEvaluator> >::const_iterator m = transientModels_.begin();
	     m != transientModels_.end(); ++m)
	  (*m)->acceptTimeStep();

	if (printTimeStepSummary_) {
	  os << "\nEnd time step " << currentTimeStep_ << ": status=" << "CONVERGED"
	     << ", time=" << currentTime_ 
	     << ", step size=" << currentStepSize_ << std::endl;
	}
      }
      else if (innerSolverStatus == FAILED) {
	if (hitMinTimeStep) {
	  timeStepStatus_ = FAILED;

	  if (printTimeStepSummary_) {
	    os << "\nEnd time step " << currentTimeStep_ << ": status=" << "FAILED"
	       << ", time=" << currentTime_ + currentStepSize_
	       << ", step size=" << currentStepSize_ << std::endl;
	    os << "Failed at minimum step size, ending simulation!" << std::endl;
	  }
	}
	// log converged and failed time steps
	++totalNumFailedSteps_;
	++numConsecutiveFailedTimeSteps_;
	numConsecutiveConvergedTimeSteps_ = 0;

	if (printTimeStepSummary_) {
	  os << "\nEnd time step " << currentTimeStep_ << ": status=" << "FAILED"
	     << ", time=" << currentTime_ + currentStepSize_
	     << ", step size=" << currentStepSize_ << std::endl;
	    os << "Reducing time step and retrying solve!" << std::endl;
	}
      }

    }


    // Check for overall status
    if ( (timeStepStatus_ == CONVERGED) && (achievedFinalTime) )
      overallStatus_ = CONVERGED;
    else if (timeStepStatus_ == FAILED)
      overallStatus_ = FAILED;
    else if (currentTimeStep_ >= maxTimeSteps_)
      overallStatus_ = FAILED;
    else {
      // Continue taking more steps
      overallStatus_ = UNCONVERGED;
    }

    return overallStatus_;
  }

  pike::SolveStatus TransientStepper::solve()
  {
    TEUCHOS_ASSERT(registrationComplete_);


    if (printTimeStepSummary_) {
      Teuchos::RCP<Teuchos::FancyOStream> os = this->getOStream();
      *os << "****************************************" << std::endl;
      *os << " Starting Transient Solve!" << std::endl;
      *os << "****************************************" << std::endl;
    }

    overallStatus_ = pike::UNCONVERGED;

    while ( (overallStatus_ != CONVERGED) && (overallStatus_ != FAILED) )
      overallStatus_ = this->step();

    if (printTimeStepSummary_) {
      Teuchos::RCP<Teuchos::FancyOStream> os = this->getOStream();
      *os << "****************************************" << std::endl;
      *os << " Completed Transient Solve: ";
      if (overallStatus_ == CONVERGED)
	*os << "Converged!";
      else
	*os << "Failed!";
      *os << std::endl;
      *os << "****************************************" << std::endl;
    }
    
    return overallStatus_;
  }

  void TransientStepper::finalize()
  { solver_->finalize(); }
  
  void TransientStepper::reset()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
			       "Error in pike::TransientStepper::reset() - The reset() method has been called on a Transient Solver.  This is currently not allowed!");
  }

  pike::SolveStatus TransientStepper::getStatus() const
  { return overallStatus_; }

  int TransientStepper::getNumberOfIterations() const
  { return currentTimeStep_; }

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
    
    maxTimeSteps_ = paramList->get<int>("Maximum Number of Time Steps");
    beginTime_ = paramList->get<double>("Begin Time");
    endTime_ = paramList->get<double>("End Time");
    currentTime_ = beginTime_;
    initialStepSize_ = paramList->get<double>("Initial Time Step Size");
    currentStepSize_ = initialStepSize_;
    minStepSize_ = paramList->get<double>("Minimum Time Step Size");
    maxStepSize_ = paramList->get<double>("Maximum Time Step Size");
    stepGrowthFactor_ = paramList->get<double>("Time Step Size Growth Factor");
    stepDecreaseFactor_ = paramList->get<double>("Time Step Size Decrease Factor");
    checkPointsVec_ = paramList->get<Teuchos::Array<double> >("Check Points");
    printTimeStepSummary_ = paramList->get<bool>("Print Time Step Summary");
    printTimeStepDetails_ = paramList->get<bool>("Print Time Step Details");
    numConvergedTimeStepsBeforeGrowth_ = paramList->get<int>("Number Converged Time Steps for Growth");

    TEUCHOS_ASSERT(beginTime_ < endTime_);
    TEUCHOS_ASSERT(minStepSize_ < maxStepSize_);
    TEUCHOS_ASSERT(initialStepSize_ >= minStepSize_);
    TEUCHOS_ASSERT(initialStepSize_ <= maxStepSize_);
    TEUCHOS_ASSERT(stepGrowthFactor_ >= 1.0);
    TEUCHOS_ASSERT(stepDecreaseFactor_ <= 1.0);

    std::sort(checkPointsVec_.begin(),checkPointsVec_.end());
    for (Teuchos::Array<double>::size_type i=0; i < checkPointsVec_.size(); ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(checkPointsVec_[i] < beginTime_,
				 std::runtime_error,
				 "");
      TEUCHOS_TEST_FOR_EXCEPTION(checkPointsVec_[i] > endTime_,
				 std::runtime_error,
				 "");
    }
    checkPoints_.assign(checkPointsVec_.begin(),checkPointsVec_.end());

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
