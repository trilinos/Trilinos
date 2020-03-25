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
  RCP<Stepper<Scalar> > stepper = sf->createStepper(stepperType, model);

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
IntegratorBasic<Scalar>::IntegratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                inputPL,
  std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models)
    : integratorObserver_(Teuchos::null),
      integratorStatus_(WORKING), isInitialized_(false)
{
  this->setTempusParameterList(inputPL);
  this->setStepper(models);
  this->initialize();
}


template<class Scalar>
void IntegratorBasic<Scalar>::setStepper(
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  if (stepper_ == Teuchos::null) {
    // Construct from Integrator ParameterList
    RCP<StepperFactory<Scalar> > sf =Teuchos::rcp(new StepperFactory<Scalar>());
    std::string stepperName = integratorPL_->get<std::string>("Stepper Name");

    RCP<ParameterList> stepperPL = Teuchos::sublist(tempusPL_,stepperName,true);
    stepper_ = sf->createStepper(stepperPL, model);
  } else {
    stepper_->setModel(model);
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::setStepper(
  std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  if (stepper_ == Teuchos::null) {
    // Construct from Integrator ParameterList
    RCP<StepperFactory<Scalar> > sf =Teuchos::rcp(new StepperFactory<Scalar>());
    std::string stepperName = integratorPL_->get<std::string>("Stepper Name");

    RCP<ParameterList> stepperPL = Teuchos::sublist(tempusPL_,stepperName,true);
    stepper_ = sf->createMultiSteppers(stepperPL, models);
  } else {
    stepper_->createSubSteppers(models);
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::setStepperWStepper(
  Teuchos::RCP<Stepper<Scalar> > newStepper)
{
  stepper_ = newStepper;
}

/// This resets the SolutionHistory and sets the first SolutionState as the IC.
template<class Scalar>
void IntegratorBasic<Scalar>::
initializeSolutionHistory(Teuchos::RCP<SolutionState<Scalar> > state)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Construct from Integrator ParameterList
  RCP<ParameterList> shPL =
    Teuchos::sublist(integratorPL_, "Solution History", true);
  solutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));

  if (state == Teuchos::null) {
    // Construct default IC from the application model
    RCP<SolutionState<Scalar> > newState = rcp(new SolutionState<Scalar>(
      stepper_->getModel(), stepper_->getDefaultStepperState(), Teuchos::null));

    // Set SolutionState from TimeStepControl
    newState->setTime    (timeStepControl_->getInitTime());
    newState->setIndex   (timeStepControl_->getInitIndex());
    newState->setTimeStep(timeStepControl_->getInitTimeStep());
    newState->setTolRel  (timeStepControl_->getMaxRelError());
    newState->setTolAbs  (timeStepControl_->getMaxAbsError());
    int order = timeStepControl_->getInitOrder();
    if (order == 0) order = stepper_->getOrder();
    if (order < stepper_->getOrderMin()) order = stepper_->getOrderMin();
    if (order > stepper_->getOrderMax()) order = stepper_->getOrderMax();
    newState->setOrder(order);
    newState->setSolutionStatus(Status::PASSED);  // ICs are considered passing.

    solutionHistory_->addState(newState);

  } else {
    // Use state as IC
    solutionHistory_->addState(state);
  }

  // Get IC from the application model via the stepper and ensure consistency.
  stepper_->setInitialConditions(solutionHistory_);
}


template<class Scalar>
void IntegratorBasic<Scalar>::
initializeSolutionHistory(Scalar t0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Construct from Integrator ParameterList
  RCP<ParameterList> shPL =
    Teuchos::sublist(integratorPL_, "Solution History", true);
  solutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));

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
    x0->clone_v(), xdot, xdotdot, stepper_->getDefaultStepperState()));

  // Set SolutionState from TimeStepControl
  newState->setTime    (t0);
  newState->setIndex   (timeStepControl_->getInitIndex());
  newState->setTimeStep(timeStepControl_->getInitTimeStep());
  int order = timeStepControl_->getInitOrder();
  if (order == 0) order = stepper_->getOrder();
  if (order < stepper_->getOrderMin()) order = stepper_->getOrderMin();
  if (order > stepper_->getOrderMax()) order = stepper_->getOrderMax();
  newState->setOrder(order);

  newState->setSolutionStatus(Status::PASSED);  // ICs are considered passing.

  solutionHistory_->addState(newState);

  // Get IC from the application model via the stepper and ensure consistency.
  stepper_->setInitialConditions(solutionHistory_);
}


template<class Scalar>
void IntegratorBasic<Scalar>::setSolutionHistory(
  Teuchos::RCP<SolutionHistory<Scalar> > sh)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  if (sh == Teuchos::null) {
    // Create default SolutionHistory, otherwise keep current history.
    if (solutionHistory_ == Teuchos::null) initializeSolutionHistory();
  } else {

    TEUCHOS_TEST_FOR_EXCEPTION( sh->getNumStates() < 1,
      std::out_of_range,
         "Error - setSolutionHistory requires at least one SolutionState.\n"
      << "        Supplied SolutionHistory has only " << sh->getNumStates()
      << " SolutionStates.\n");

    // Make integratorPL_ consistent with new SolutionHistory.
    RCP<ParameterList> shPL = sh->getNonconstParameterList();
    integratorPL_->set("Solution History", shPL->name());
    integratorPL_->set(shPL->name(), *shPL);

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
    // Create timeStepControl_ if null, otherwise keep current parameters.
    if (timeStepControl_ == Teuchos::null) {
      if (integratorPL_->isSublist("Time Step Control")) {
        // Construct from Integrator ParameterList
        RCP<ParameterList> tscPL =
          Teuchos::sublist(integratorPL_,"Time Step Control",true);
        timeStepControl_ = rcp(new TimeStepControl<Scalar>(tscPL));
      } else {
        // Construct default TimeStepControl
        timeStepControl_ = rcp(new TimeStepControl<Scalar>());
        RCP<ParameterList> tscPL = timeStepControl_->getNonconstParameterList();
        integratorPL_->set("Time Step Control", tscPL->name());
        integratorPL_->set(tscPL->name(), *tscPL);
      }
    }

  } else {
    // Make integratorPL_ consistent with new TimeStepControl.
    RCP<ParameterList> tscPL = tsc->getNonconstParameterList();
    integratorPL_->set("Time Step Control", tscPL->name());
    integratorPL_->set(tscPL->name(), *tscPL);
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
        Teuchos::rcp(new IntegratorObserverComposite<Scalar>);
      // Add basic observer to output time step info
      Teuchos::RCP<IntegratorObserverBasic<Scalar> > basicObs =
          Teuchos::rcp(new IntegratorObserverBasic<Scalar>);
      integratorObserver_->addObserver(basicObs);
    }
  } else {
    if (integratorObserver_ == Teuchos::null) {
      integratorObserver_ =
        Teuchos::rcp(new IntegratorObserverComposite<Scalar>);
    }
    integratorObserver_->addObserver(obs);
  }

}


template<class Scalar>
void IntegratorBasic<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( stepper_ == Teuchos::null, std::logic_error,
    "Error - Need to set the Stepper, setStepper(), before calling "
    "IntegratorBasic::initialize()\n");

  this->setTimeStepControl();
  this->parseScreenOutput();
  this->setSolutionHistory();
  this->setObserver();

  // Set initial conditions, make them consistent, and set stepper memory.
  stepper_->setInitialConditions(solutionHistory_);

  // Ensure TimeStepControl orders match the Stepper orders.
  if (timeStepControl_->getMinOrder() < stepper_->getOrderMin())
      timeStepControl_->setMinOrder(stepper_->getOrderMin());
  if (timeStepControl_->getMinOrder() > stepper_->getOrderMax())
      timeStepControl_->setMinOrder(stepper_->getOrderMax());

  if (timeStepControl_->getMaxOrder() == 0 ||
      timeStepControl_->getMaxOrder() > stepper_->getOrderMax())
      timeStepControl_->setMaxOrder(stepper_->getOrderMax());
  if (timeStepControl_->getMaxOrder() < timeStepControl_->getMinOrder())
      timeStepControl_->setMaxOrder(timeStepControl_->getMinOrder());

  if (timeStepControl_->getInitOrder() < timeStepControl_->getMinOrder())
      timeStepControl_->setInitOrder(timeStepControl_->getMinOrder());
  if (timeStepControl_->getInitOrder() > timeStepControl_->getMaxOrder())
      timeStepControl_->setInitOrder(timeStepControl_->getMaxOrder());

  TEUCHOS_TEST_FOR_EXCEPTION(
    timeStepControl_->getMinOrder() > timeStepControl_->getMaxOrder(),
    std::out_of_range,
       "Error - Invalid TimeStepControl min order greater than max order.\n"
    << "        Min order = " << timeStepControl_->getMinOrder() << "\n"
    << "        Max order = " << timeStepControl_->getMaxOrder() << "\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
    timeStepControl_->getInitOrder() < timeStepControl_->getMinOrder() ||
    timeStepControl_->getInitOrder() > timeStepControl_->getMaxOrder(),
    std::out_of_range,
       "Error - Initial TimeStepControl order is out of min/max range.\n"
    << "        Initial order = " << timeStepControl_->getInitOrder() << "\n"
    << "        Min order     = " << timeStepControl_->getMinOrder()  << "\n"
    << "        Max order     = " << timeStepControl_->getMaxOrder()  << "\n");

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
    integratorStatus_ = Status::FAILED;
    return;
  }
  integratorTimer_->start();
  // get optimal initial time step
  const Scalar initDt =
     std::min(timeStepControl_->getInitTimeStep(),
              stepper_->getInitTimeStep(solutionHistory_));
  // update initial time step
  timeStepControl_->setInitTimeStep(initDt);
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

      if (integratorStatus_ == Status::FAILED) break;
      solutionHistory_->getWorkingState()->setSolutionStatus(WORKING);

      integratorObserver_->observeBeforeTakeStep(*this);

      stepper_->takeStep(solutionHistory_);

      integratorObserver_->observeAfterTakeStep(*this);

      stepperTimer_->stop();
      checkTimeStep();
      integratorObserver_->observeAfterCheckTimeStep(*this);

      solutionHistory_->promoteWorkingState();
      integratorObserver_->observeEndTimeStep(*this);
    }

    endIntegrator();
    integratorObserver_->observeEndIntegrator(*this);
  }

  return (integratorStatus_ == Status::PASSED);
}


template <class Scalar>
void IntegratorBasic<Scalar>::startTimeStep()
{
  auto ws = solutionHistory_->getWorkingState();

  //set the Abs/Rel tolerance
  ws->setTolRel(timeStepControl_->getMaxRelError());
  ws->setTolAbs(timeStepControl_->getMaxAbsError());

  // Check if we need to dump screen output this step
  std::vector<int>::const_iterator it =
    std::find(outputScreenIndices_.begin(),
              outputScreenIndices_.end(),
              ws->getIndex());
  if (it == outputScreenIndices_.end())
    ws->setOutputScreen(false);
  else
    ws->setOutputScreen(true);

  const int initial = timeStepControl_->getInitIndex();
  const int interval = integratorPL_->get<int>("Screen Output Index Interval");
  if ( (ws->getIndex() - initial) % interval == 0)
    ws->setOutputScreen(true);
}


template <class Scalar>
void IntegratorBasic<Scalar>::checkTimeStep()
{
  using Teuchos::RCP;
  auto ws = solutionHistory_->getWorkingState();

  // Too many TimeStep failures, Integrator fails.
  if (ws->getNFailures() >= timeStepControl_->getMaxFailures()) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,2,"checkTimeStep");
    *out << "Failure - Stepper has failed more than the maximum allowed.\n"
         << "  (nFailures = "<<ws->getNFailures()<< ") >= (nFailuresMax = "
         << timeStepControl_->getMaxFailures()<<")" << std::endl;
    integratorStatus_ = Status::FAILED;
    return;
  }
  if (ws->getNConsecutiveFailures()
      >= timeStepControl_->getMaxConsecFailures()){
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"checkTimeStep");
    *out << "Failure - Stepper has failed more than the maximum "
         << "consecutive allowed.\n"
         << "  (nConsecutiveFailures = "<<ws->getNConsecutiveFailures()
         << ") >= (nConsecutiveFailuresMax = "
         << timeStepControl_->getMaxConsecFailures()
         << ")" << std::endl;
    integratorStatus_ = Status::FAILED;
    return;
  }

  // Check Stepper failure.
  if (ws->getSolutionStatus() == Status::FAILED or
       // Constant time step failure
       ((timeStepControl_->getStepType() == "Constant") and
        (ws->getTimeStep() != timeStepControl_->getInitTimeStep()) and
        (ws->getOutput() != true) and
        (ws->getTime() != timeStepControl_->getFinalTime())
       )
     )
  {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,0,"checkTimeStep");
    *out <<std::scientific
      <<std::setw( 6)<<std::setprecision(3)<<ws->getIndex()
      <<std::setw(11)<<std::setprecision(3)<<ws->getTime()
      <<std::setw(11)<<std::setprecision(3)<<ws->getTimeStep()
      << "  STEP FAILURE!! - ";
    if (ws->getSolutionStatus() == Status::FAILED) {
      *out << "Solution Status = " << toString(ws->getSolutionStatus())
           << std::endl;
    } else if ((timeStepControl_->getStepType() == "Constant") and
               (ws->getTimeStep() != timeStepControl_->getInitTimeStep())) {
      *out << "dt != Constant dt (="<<timeStepControl_->getInitTimeStep()<<")"
           << std::endl;
    }

    ws->setNFailures(ws->getNFailures()+1);
    ws->setNRunningFailures(ws->getNRunningFailures()+1);
    ws->setNConsecutiveFailures(ws->getNConsecutiveFailures()+1);
    ws->setSolutionStatus(Status::FAILED);
    return;
  }

  // TIME STEP PASSED basic tests!  Ensure it is set as such.
  ws->setSolutionStatus(Status::PASSED);

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
  runtime_ = integratorTimer_->totalElapsedTime();
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

  // order output indices and remove duplicates.
  std::sort(outputScreenIndices_.begin(),outputScreenIndices_.end());
  outputScreenIndices_.erase(std::unique(outputScreenIndices_.begin(),
                                         outputScreenIndices_.end()   ),
                                         outputScreenIndices_.end()     );
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

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                     pList,
  std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models)
{
  Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorBasic<Scalar>(pList, models));
  return(integrator);
}

} // namespace Tempus
#endif // Tempus_IntegratorBasic_impl_hpp
