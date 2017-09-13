// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorBasic_impl_hpp
#define Tempus_IntegratorBasic_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <ctime>


namespace Tempus {

template<class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                inputPL,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model)
    : integratorObserver_(Teuchos::null),
      integratorStatus_(WORKING), isInitialized_(false)
{
  this->setTempusParameterList(inputPL);
  this->setStepper(model);
  this->initialize();
}


template<class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  std::string stepperType)
    : integratorObserver_(Teuchos::null),
      integratorStatus_(WORKING), isInitialized_(false)
{
  using Teuchos::RCP;
  RCP<StepperFactory<Scalar> > sf = Teuchos::rcp(new StepperFactory<Scalar>());
  RCP<Stepper<Scalar> > stepper = sf->createStepper(model, stepperType);

  this->setTempusParameterList(Teuchos::null);
  this->setStepperWStepper(stepper);
  this->initialize();
}


template<class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic()
  : integratorObserver_(Teuchos::null),
    integratorStatus_(WORKING), isInitialized_(false)
{
  this->setTempusParameterList(Teuchos::null);
}


template<class Scalar>
void IntegratorBasic<Scalar>::setStepper(
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  if (stepper_ == Teuchos::null) {
    // Construct from Integrator ParameterList
    RCP<StepperFactory<Scalar> > sf =Teuchos::rcp(new StepperFactory<Scalar>());
    std::string stepperName = integratorPL_->get<std::string>("Stepper Name");

    RCP<ParameterList> stepperPL = Teuchos::sublist(tempusPL_,stepperName,true);
    stepper_ = sf->createStepper(model, stepperPL);
  } else {
    stepper_->setModel(model);
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::setStepperWStepper(
  Teuchos::RCP<Stepper<Scalar> > newStepper)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Make integratorPL_ consistent with new stepper.
  RCP<ParameterList> newStepperPL = newStepper->getNonconstParameterList();
  integratorPL_->set("Stepper Name", newStepperPL->name());
  tempusPL_->set(newStepperPL->name(), *newStepperPL);
  stepper_ = newStepper;
}

/// This resets the SolutionHistory and sets the first SolutionState as the IC.
template<class Scalar>
void IntegratorBasic<Scalar>::
setInitialState(Teuchos::RCP<SolutionState<Scalar> >  state)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Construct from Integrator ParameterList
  RCP<ParameterList> shPL =
    Teuchos::sublist(integratorPL_, "Solution History", true);
  solutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));

  if (state == Teuchos::null) {
    // Construct default IC
    // Create initial condition metadata from TimeStepControl
    RCP<SolutionStateMetaData<Scalar> > md =
      rcp(new SolutionStateMetaData<Scalar> ());
    md->setTime (timeStepControl_->getInitTime());
    md->setIStep(timeStepControl_->getInitIndex());
    md->setDt   (timeStepControl_->getInitTimeStep());
    int orderTmp = timeStepControl_->getInitOrder();
    if (orderTmp == 0) orderTmp = stepper_->getOrder();
    md->setOrder(orderTmp);
    md->setSolutionStatus(Status::PASSED);  // ICs are considered passing.

    // Create initial condition from ModelEvaluator::getNominalValues()
    typedef Thyra::ModelEvaluatorBase MEB;
    Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgsIC =
      stepper_->getModel()->getNominalValues();
    RCP<Thyra::VectorBase<Scalar> > x = inArgsIC.get_x()->clone_v();
    RCP<Thyra::VectorBase<Scalar> > xdot;
    if (inArgsIC.supports(MEB::IN_ARG_x_dot)) {
      xdot = inArgsIC.get_x_dot()->clone_v();
    } else {
      xdot = x->clone_v();
    }
    RCP<Thyra::VectorBase<Scalar> > xdotdot;
    if (inArgsIC.supports(MEB::IN_ARG_x_dot_dot)) {
      xdotdot = inArgsIC.get_x_dot_dot()->clone_v();
    }
    else {
      xdotdot = x->clone_v();
    }
    RCP<SolutionState<Scalar> > newState = rcp(new SolutionState<Scalar>(
      md, x, xdot, xdotdot, stepper_->getDefaultStepperState()));

    solutionHistory_->addState(newState);
  } else {
    // Use state as IC
    solutionHistory_->addState(state);
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::
setInitialState(Scalar t0,
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x0,
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot0,
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xdotdot0)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Construct from Integrator ParameterList
  RCP<ParameterList> shPL =
    Teuchos::sublist(integratorPL_, "Solution History", true);
  solutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));

  // Create initial condition metadata from TimeStepControl
  RCP<SolutionStateMetaData<Scalar> > md =
    rcp(new SolutionStateMetaData<Scalar> ());
  md->setTime (timeStepControl_->getInitTime());
  md->setIStep(timeStepControl_->getInitIndex());
  md->setDt   (timeStepControl_->getInitTimeStep());
  int orderTmp = timeStepControl_->getInitOrder();
  if (orderTmp == 0) orderTmp = stepper_->getOrder();
  md->setOrder(orderTmp);
  md->setSolutionStatus(Status::PASSED);  // ICs are considered passing.

  // Create and set xdot and xdotdot.
  RCP<Thyra::VectorBase<Scalar> > xdot    = x0->clone_v();
  RCP<Thyra::VectorBase<Scalar> > xdotdot = x0->clone_v();
  if (xdot0 == Teuchos::null)
    Thyra::assign(xdot.ptr(),    Teuchos::ScalarTraits<Scalar>::zero());
  else
    Thyra::assign(xdot.ptr(),    *(xdot0));
  if (xdotdot0 == Teuchos::null)
    Thyra::assign(xdotdot.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  else
    Thyra::assign(xdotdot.ptr(), *(xdotdot0));

  RCP<SolutionState<Scalar> > newState = rcp(new SolutionState<Scalar>(
    md, x0, xdot, xdotdot, stepper_->getDefaultStepperState()));

  solutionHistory_->addState(newState);
}


template<class Scalar>
void IntegratorBasic<Scalar>::setSolutionHistory(
  Teuchos::RCP<SolutionHistory<Scalar> > sh)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  if (sh == Teuchos::null) {
    // Create default SolutionHistory, otherwise keep current history.
    if (solutionHistory_ == Teuchos::null) setInitialState();
  } else {

    TEUCHOS_TEST_FOR_EXCEPTION( sh->getNumStates() < 1,
      std::out_of_range,
         "Error - setSolutionHistory requires at least one SolutionState.\n"
      << "        Supplied SolutionHistory has only " << sh->getNumStates()
      << " SolutionStates.\n");

    // Make integratorPL_ consistent with new SolutionHistory.
    RCP<ParameterList> shPL = sh->getNonconstParameterList();
    integratorPL_->set("Solution History", shPL->name());
    integratorPL_->set(shPL->name(), shPL);

    solutionHistory_ = sh;
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::setTimeStepControl(
  Teuchos::RCP<TimeStepControl<Scalar> > tsc)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  if (tsc == Teuchos::null) {
    // Construct from Integrator ParameterList
    RCP<ParameterList> tscPL =
      Teuchos::sublist(integratorPL_,"Time Step Control",true);
    timeStepControl_ = rcp(new TimeStepControl<Scalar>(tscPL));

    if (timeStepControl_->getMinOrder() == 0)
      timeStepControl_->setMinOrder(stepper_->getOrderMin());
    if (timeStepControl_->getMaxOrder() == 0)
      timeStepControl_->setMaxOrder(stepper_->getOrderMax());
    if (timeStepControl_->getInitOrder() < timeStepControl_->getMinOrder())
      timeStepControl_->setInitOrder(timeStepControl_->getMinOrder());
    if (timeStepControl_->getInitOrder() > timeStepControl_->getMaxOrder())
      timeStepControl_->setInitOrder(timeStepControl_->getMaxOrder());

  } else {
    // Make integratorPL_ consistent with new TimeStepControl.
    RCP<ParameterList> tscPL = tsc->getNonconstParameterList();
    integratorPL_->set("Time Step Control", tscPL->name());
    integratorPL_->set(tscPL->name(), tscPL);

    timeStepControl_ = Teuchos::null;
    timeStepControl_ = tsc;
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::setObserver(
  Teuchos::RCP<IntegratorObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default IntegratorObserverBasic, otherwise keep current observer.
    if (integratorObserver_ == Teuchos::null) {
      integratorObserver_ =
        Teuchos::rcp(new IntegratorObserverBasic<Scalar>);
    }
  } else {
    integratorObserver_ = obs;
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::initialize()
{
  this->setTimeStepControl();
  this->parseScreenOutput();
  this->setSolutionHistory();
  this->setObserver();

  if (integratorTimer_ == Teuchos::null)
    integratorTimer_ = rcp(new Teuchos::Time("Integrator Timer"));
  if (stepperTimer_ == Teuchos::null)
    stepperTimer_    = rcp(new Teuchos::Time("Stepper Timer"));

  if (Teuchos::as<int>(this->getVerbLevel()) >=
      Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"IntegratorBasic::IntegratorBasic");
    *out << this->description() << std::endl;
  }

  isInitialized_ = true;
}


template<class Scalar>
std::string IntegratorBasic<Scalar>::description() const
{
  std::string name = "Tempus::IntegratorBasic";
  return(name);
}


template<class Scalar>
void IntegratorBasic<Scalar>::describe(
  Teuchos::FancyOStream          &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  out << description() << "::describe" << std::endl;
  out << "solutionHistory= " << solutionHistory_->description()<<std::endl;
  out << "timeStepControl= " << timeStepControl_->description()<<std::endl;
  out << "stepper        = " << stepper_        ->description()<<std::endl;

  if (Teuchos::as<int>(verbLevel) >=
              Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    out << "solutionHistory= " << std::endl;
    solutionHistory_->describe(out,verbLevel);
    out << "timeStepControl= " << std::endl;
    timeStepControl_->describe(out,verbLevel);
    out << "stepper        = " << std::endl;
    stepper_        ->describe(out,verbLevel);
  }
}


template <class Scalar>
bool IntegratorBasic<Scalar>::advanceTime(const Scalar timeFinal)
{
  if (timeStepControl_->timeInRange(timeFinal))
    timeStepControl_->setFinalTime(timeFinal);
  bool itgrStatus = advanceTime();
  return itgrStatus;
}


template <class Scalar>
void IntegratorBasic<Scalar>::startIntegrator()
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  if (isInitialized_ == false) {
    Teuchos::OSTab ostab(out,1,"StartIntegrator");
    *out << "Failure - IntegratorBasic is not initialized." << std::endl;
    integratorStatus_ = FAILED;
    return;
  }
  std::time_t begin = std::time(nullptr);
  integratorTimer_->start();
  Teuchos::OSTab ostab(out,0,"ScreenOutput");
  *out << "\nTempus - IntegratorBasic\n"
       << std::asctime(std::localtime(&begin)) << "\n"
       << "  Stepper = " << stepper_->description() << "\n"
       << "  Simulation Time Range  [" << timeStepControl_->getInitTime()
       << ", " << timeStepControl_->getFinalTime() << "]\n"
       << "  Simulation Index Range [" << timeStepControl_->getInitIndex()
       << ", " << timeStepControl_->getFinalIndex() << "]\n"
       << "============================================================================\n"
       << "  Step       Time         dt  Abs Error  Rel Error  Order  nFail  dCompTime"
       << std::endl;
  integratorStatus_ = WORKING;
}


template <class Scalar>
bool IntegratorBasic<Scalar>::advanceTime()
{
  TEMPUS_FUNC_TIME_MONITOR("Tempus::IntegratorBasic::advanceTime()");
  {
    startIntegrator();
    integratorObserver_->observeStartIntegrator(*this);

    while (integratorStatus_ == WORKING and
        timeStepControl_->timeInRange (solutionHistory_->getCurrentTime()) and
        timeStepControl_->indexInRange(solutionHistory_->getCurrentIndex())){

      stepperTimer_->reset();
      stepperTimer_->start();
      solutionHistory_->initWorkingState();

      startTimeStep();
      integratorObserver_->observeStartTimeStep(*this);

      timeStepControl_->getNextTimeStep(solutionHistory_, integratorStatus_);
      integratorObserver_->observeNextTimeStep(*this);

      if (integratorStatus_ == FAILED) break;
      solutionHistory_->getWorkingState()->getMetaData()->
        setSolutionStatus(WORKING);

      integratorObserver_->observeBeforeTakeStep(*this);

      stepper_->takeStep(solutionHistory_);

      integratorObserver_->observeAfterTakeStep(*this);

      stepperTimer_->stop();
      acceptTimeStep();
      integratorObserver_->observeAcceptedTimeStep(*this);
    }

    endIntegrator();
    integratorObserver_->observeEndIntegrator(*this);
  }

  return (integratorStatus_ == Status::PASSED);
}


template <class Scalar>
void IntegratorBasic<Scalar>::startTimeStep()
{
  Teuchos::RCP<SolutionStateMetaData<Scalar> > wsmd =
    solutionHistory_->getWorkingState()->getMetaData();

  // Check if we need to dump screen output this step
  std::vector<int>::const_iterator it =
    std::find(outputScreenIndices_.begin(),
              outputScreenIndices_.end(),
              wsmd->getIStep()+1);
  if (it == outputScreenIndices_.end())
    wsmd->setOutputScreen(false);
  else
    wsmd->setOutputScreen(true);
}


template <class Scalar>
void IntegratorBasic<Scalar>::acceptTimeStep()
{
  using Teuchos::RCP;
  RCP<SolutionStateMetaData<Scalar> > wsmd =
    solutionHistory_->getWorkingState()->getMetaData();

  // Too many failures
  if (wsmd->getNFailures() >= timeStepControl_->getMaxFailures()) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,2,"acceptTimeStep");
    *out << "Failure - Stepper has failed more than the maximum allowed.\n"
         << "  (nFailures = "<<wsmd->getNFailures()<< ") >= (nFailuresMax = "
         << timeStepControl_->getMaxFailures()<<")" << std::endl;
    integratorStatus_ = FAILED;
    return;
  }
  if (wsmd->getNConsecutiveFailures()
      >= timeStepControl_->getMaxConsecFailures()){
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"acceptTimeStep");
    *out << "Failure - Stepper has failed more than the maximum "
         << "consecutive allowed.\n"
         << "  (nConsecutiveFailures = "<<wsmd->getNConsecutiveFailures()
         << ") >= (nConsecutiveFailuresMax = "
         << timeStepControl_->getMaxConsecFailures()
         << ")" << std::endl;
    integratorStatus_ = FAILED;
    return;
  }

       // Stepper failure
  if ( solutionHistory_->getWorkingState()->getSolutionStatus() == FAILED or
       solutionHistory_->getWorkingState()->getStepperStatus() == FAILED or
       // Constant time step failure
       ((timeStepControl_->getStepType() == "Constant") and
        (wsmd->getDt() != timeStepControl_->getInitTimeStep()) and
        (wsmd->getOutput() != true) and
        (wsmd->getTime()+wsmd->getDt() != timeStepControl_->getFinalTime())
       )
     )
  {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,0,"acceptTimeStep");
    *out <<std::scientific
      <<std::setw( 6)<<std::setprecision(3)<<wsmd->getIStep()
      <<std::setw(11)<<std::setprecision(3)<<wsmd->getTime()
      <<std::setw(11)<<std::setprecision(3)<<wsmd->getDt()
      << "  STEP FAILURE!! - ";
    if ( solutionHistory_->getWorkingState()->getSolutionStatus() == FAILED) {
      *out << "Solution Status = "
           << toString(solutionHistory_->getWorkingState()->getSolutionStatus())
           << std::endl;
    } else if (solutionHistory_->getWorkingState()->getStepperStatus()==FAILED){
      *out << "Stepper Status = "
           << toString(solutionHistory_->getWorkingState()->getStepperStatus())
           << std::endl;
    } else if ((timeStepControl_->getStepType() == "Constant") and
               (wsmd->getDt() != timeStepControl_->getInitTimeStep())) {
      *out << "dt != Constant dt (="<<timeStepControl_->getInitTimeStep()<<")"
           << std::endl;
    }

    wsmd->setNFailures(wsmd->getNFailures()+1);
    wsmd->setNConsecutiveFailures(wsmd->getNConsecutiveFailures()+1);
    wsmd->setSolutionStatus(FAILED);
    return;
  }

  // =======================================================================
  // Made it here! Accept this time step

  solutionHistory_->promoteWorkingState();

  RCP<SolutionStateMetaData<Scalar> > csmd =
    solutionHistory_->getCurrentState()->getMetaData();

  csmd->setNFailures(std::max(csmd->getNFailures()-1,0));
  csmd->setNConsecutiveFailures(0);

  if ((csmd->getOutputScreen() == true) or
      (csmd->getOutput() == true) or
      (csmd->getTime() == timeStepControl_->getFinalTime())) {
    const double steppertime = stepperTimer_->totalElapsedTime();
    stepperTimer_->reset();
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,0,"ScreenOutput");
    *out
    <<std::scientific<<std::setw( 6)<<std::setprecision(3)<<csmd->getIStep()
                     <<std::setw(11)<<std::setprecision(3)<<csmd->getTime()
                     <<std::setw(11)<<std::setprecision(3)<<csmd->getDt()
                     <<std::setw(11)<<std::setprecision(3)<<csmd->getErrorAbs()
                     <<std::setw(11)<<std::setprecision(3)<<csmd->getErrorRel()
    <<std::fixed     <<std::setw( 7)<<std::setprecision(1)<<csmd->getOrder()
    <<std::scientific<<std::setw( 7)<<std::setprecision(3)<<csmd->getNFailures()
                     <<std::setw(11)<<std::setprecision(3)<<steppertime
    <<std::endl;
  }

  // Output and screen output
  if (csmd->getOutput() == true) {
    // Dump solution!
  }
}


template <class Scalar>
void IntegratorBasic<Scalar>::endIntegrator()
{
  std::string exitStatus;
  if (solutionHistory_->getCurrentState()->getSolutionStatus() ==
      Status::FAILED or integratorStatus_ == Status::FAILED) {
    exitStatus = "Time integration FAILURE!";
  } else {
    integratorStatus_ = Status::PASSED;
    exitStatus = "Time integration complete.";
  }

  integratorTimer_->stop();
  const double runtime = integratorTimer_->totalElapsedTime();
  std::time_t end = std::time(nullptr);
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,0,"ScreenOutput");
  *out << "============================================================================\n"
       << "  Total runtime = " << runtime << " sec = "
       << runtime/60.0 << " min\n"
       << std::asctime(std::localtime(&end))
       << exitStatus << "\n"
       << std::endl;
}


template <class Scalar>
void IntegratorBasic<Scalar>::parseScreenOutput()
{
  // This has been delayed until timeStepControl has been constructed.

  // Parse output indices
  outputScreenIndices_.clear();
  std::string str =
    integratorPL_->get<std::string>("Screen Output Index List", "");
  std::string delimiters(",");
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
    std::string token = str.substr(lastPos,pos-lastPos);
    outputScreenIndices_.push_back(int(std::stoi(token)));
    if(pos==std::string::npos)
      break;

    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }

  int outputScreenIndexInterval =
    integratorPL_->get<int>("Screen Output Index Interval", 1000000);
  int outputScreen_i   = timeStepControl_->getInitIndex();
  const int finalIStep = timeStepControl_->getFinalIndex();
  while (outputScreen_i <= finalIStep) {
    outputScreenIndices_.push_back(outputScreen_i);
    outputScreen_i += outputScreenIndexInterval;
  }

  // order output indices
  std::sort(outputScreenIndices_.begin(),outputScreenIndices_.end());

  return;
}


template <class Scalar>
void IntegratorBasic<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & inputPL)
{
  if (inputPL == Teuchos::null) {
    if (tempusPL_->isParameter("Integrator Name")) {
      // Set Integrator PL from Tempus PL
      std::string integratorName_ =
         tempusPL_->get<std::string>("Integrator Name");
      integratorPL_ = Teuchos::sublist(tempusPL_, integratorName_, true);
    } else {
      // Build default Integrator PL
      integratorPL_ = Teuchos::parameterList();
      integratorPL_->setName("Default Integrator");
      *integratorPL_ = *(this->getValidParameters());
      tempusPL_->set("Integrator Name", "Default Integrator");
      tempusPL_->set("Default Integrator", *integratorPL_);
    }
  } else {
    *integratorPL_ = *inputPL;
    tempusPL_->set("Integrator Name", integratorPL_->name());
    tempusPL_->set(integratorPL_->name(), *integratorPL_);
  }

  integratorPL_->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string integratorType =
    integratorPL_->get<std::string>("Integrator Type");
  TEUCHOS_TEST_FOR_EXCEPTION( integratorType != "Integrator Basic",
    std::logic_error,
    "Error - Inconsistent Integrator Type for IntegratorBasic\n"
    << "    Integrator Type = " << integratorType << "\n");

  return;
}


/** \brief Create valid IntegratorBasic ParameterList.
 */
template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
IntegratorBasic<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

  std::ostringstream tmp;
  tmp << "'Integrator Type' must be 'Integrator Basic'.";
  pl->set("Integrator Type", "Integrator Basic", tmp.str());

  tmp.clear();
  tmp << "Screen Output Index List.  Required to be in TimeStepControl range "
      << "['Minimum Time Step Index', 'Maximum Time Step Index']";
  pl->set("Screen Output Index List", "", tmp.str());
  pl->set("Screen Output Index Interval", 1000000,
    "Screen Output Index Interval (e.g., every 100 time steps)");

  pl->set("Stepper Name", "",
    "'Stepper Name' selects the Stepper block to construct (Required).");

  // Solution History
  pl->sublist("Solution History",false,"solutionHistory_docs")
      .disableRecursiveValidation();

  // Time Step Control
  pl->sublist("Time Step Control",false,"solutionHistory_docs")
      .disableRecursiveValidation();

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorBasic<Scalar>::getNonconstParameterList()
{
  return(integratorPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorBasic<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = integratorPL_;
  integratorPL_ = Teuchos::null;
  return(temp_param_list);
}

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                     pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model)
{
  Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorBasic<Scalar>(pList, model));
  return(integrator);
}

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integratorBasic(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model,
  std::string stepperType)
{
  Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorBasic<Scalar>(model, stepperType));
  return(integrator);
}

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integratorBasic()
{
  Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorBasic<Scalar>());
  return(integrator);
}

} // namespace Tempus
#endif // Tempus_IntegratorBasic_impl_hpp
