//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorBasic_impl_hpp
#define Tempus_IntegratorBasic_impl_hpp

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_NumericalUtils.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperForwardEuler.hpp"

namespace Tempus {

template <class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic()
  : outputScreenIndices_(std::vector<int>()),
    outputScreenInterval_(1000000),
    integratorStatus_(WORKING),
    isInitialized_(false)
{
  setIntegratorName("Integrator Basic");
  setIntegratorType("Integrator Basic");
  setStepper(Teuchos::null);
  setSolutionHistory(Teuchos::null);
  setTimeStepControl(Teuchos::null);
  setObserver(Teuchos::null);

  integratorTimer_ = rcp(new Teuchos::Time("Integrator Timer"));
  stepperTimer_    = rcp(new Teuchos::Time("Stepper Timer"));

  // integrator is not initialized.  Still requires calls to setModel
  // and setSolutionHistory for initial conditions before calling
  // initialize() to be fully constructed.
}

template <class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic(
    Teuchos::RCP<Stepper<Scalar> > stepper,
    Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory,
    Teuchos::RCP<TimeStepControl<Scalar> > timeStepControl,
    Teuchos::RCP<IntegratorObserver<Scalar> > integratorObserver,
    std::vector<int> outputScreenIndices, int outputScreenInterval)
  : outputScreenIndices_(outputScreenIndices),
    outputScreenInterval_(outputScreenInterval),
    integratorStatus_(WORKING),
    isInitialized_(false)
{
  setIntegratorName("Integrator Basic");
  setIntegratorType("Integrator Basic");
  setStepper(stepper);
  setSolutionHistory(solutionHistory);
  setTimeStepControl(timeStepControl);
  setObserver(integratorObserver);

  integratorTimer_ = rcp(new Teuchos::Time("Integrator Timer"));
  stepperTimer_    = rcp(new Teuchos::Time("Stepper Timer"));

  initialize();
}

template <class Scalar>
void IntegratorBasic<Scalar>::copy(Teuchos::RCP<IntegratorBasic<Scalar> > iB)
{
  this->setIntegratorType(iB->getIntegratorType());
  this->setIntegratorName(iB->getIntegratorName());
  this->setStepper(iB->getStepper());
  this->setSolutionHistory(iB->getNonConstSolutionHistory());
  this->setTimeStepControl(iB->getNonConstTimeStepControl());
  this->setObserver(iB->getObserver());
  this->setScreenOutputIndexList(iB->getScreenOutputIndexList());
  this->setScreenOutputIndexInterval(iB->getScreenOutputIndexInterval());
  this->setStatus(iB->getStatus());
  integratorTimer_ = iB->getIntegratorTimer();
  stepperTimer_    = iB->getStepperTimer();
}

template <class Scalar>
void IntegratorBasic<Scalar>::setIntegratorType(std::string i)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      i != "Integrator Basic", std::logic_error,
      "Error - Integrator Type should be 'Integrator Basic'\n");

  this->integratorType_ = i;
}

template <class Scalar>
void IntegratorBasic<Scalar>::setModel(
    Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepper_ == Teuchos::null, std::logic_error,
      "Error - setModel(), need to set stepper first!\n");

  stepper_->setModel(model);
}

template <class Scalar>
void IntegratorBasic<Scalar>::setStepper(Teuchos::RCP<Stepper<Scalar> > stepper)
{
  if (stepper == Teuchos::null)
    stepper_ = Teuchos::rcp(new StepperForwardEuler<Scalar>());
  else
    stepper_ = stepper;
}

/// This resets the SolutionHistory and sets the first SolutionState as the IC.
template <class Scalar>
void IntegratorBasic<Scalar>::initializeSolutionHistory(
    Teuchos::RCP<SolutionState<Scalar> > state)
{
  using Teuchos::RCP;

  if (solutionHistory_ == Teuchos::null) {
    solutionHistory_ = rcp(new SolutionHistory<Scalar>());
  }
  else {
    solutionHistory_->clear();
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      stepper_ == Teuchos::null, std::logic_error,
      "Error - initializeSolutionHistory(), need to set stepper first!\n");

  if (state == Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(stepper_->getModel() == Teuchos::null,
                               std::logic_error,
                               "Error - initializeSolutionHistory(), need to "
                               "set stepper's model first!\n");
    // Construct default IC from the application model
    state = createSolutionStateME(stepper_->getModel(),
                                  stepper_->getDefaultStepperState());

    if (timeStepControl_ != Teuchos::null) {
      // Set SolutionState from TimeStepControl
      state->setTime(timeStepControl_->getInitTime());
      state->setIndex(timeStepControl_->getInitIndex());
      state->setTimeStep(timeStepControl_->getInitTimeStep());
      state->setTolRel(timeStepControl_->getMaxRelError());
      state->setTolAbs(timeStepControl_->getMaxAbsError());
    }
    state->setOrder(stepper_->getOrder());
    state->setSolutionStatus(Status::PASSED);  // ICs are considered passing.
  }

  solutionHistory_->addState(state);

  stepper_->setInitialConditions(solutionHistory_);
}

template <class Scalar>
void IntegratorBasic<Scalar>::initializeSolutionHistory(
    Scalar t0, Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0)
{
  using Teuchos::RCP;

  // Create and set xdot and xdotdot.
  RCP<Thyra::VectorBase<Scalar> > xdot    = x0->clone_v();
  RCP<Thyra::VectorBase<Scalar> > xdotdot = x0->clone_v();
  if (xdot0 == Teuchos::null)
    Thyra::assign(xdot.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  else
    Thyra::assign(xdot.ptr(), *(xdot0));
  if (xdotdot0 == Teuchos::null)
    Thyra::assign(xdotdot.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
  else
    Thyra::assign(xdotdot.ptr(), *(xdotdot0));

  TEUCHOS_TEST_FOR_EXCEPTION(
      stepper_ == Teuchos::null, std::logic_error,
      "Error - initializeSolutionHistory(), need to set stepper first!\n");

  auto state = createSolutionStateX(x0->clone_v(), xdot, xdotdot);
  state->setStepperState(stepper_->getDefaultStepperState());

  state->setTime(t0);
  if (timeStepControl_ != Teuchos::null) {
    // Set SolutionState from TimeStepControl
    state->setIndex(timeStepControl_->getInitIndex());
    state->setTimeStep(timeStepControl_->getInitTimeStep());
    state->setTolRel(timeStepControl_->getMaxRelError());
    state->setTolAbs(timeStepControl_->getMaxAbsError());
  }
  state->setOrder(stepper_->getOrder());
  state->setSolutionStatus(Status::PASSED);  // ICs are considered passing.

  initializeSolutionHistory(state);
}

template <class Scalar>
void IntegratorBasic<Scalar>::setSolutionHistory(
    Teuchos::RCP<SolutionHistory<Scalar> > sh)
{
  if (sh == Teuchos::null) {
    // Create default SolutionHistory, otherwise keep current history.
    if (solutionHistory_ == Teuchos::null)
      solutionHistory_ = rcp(new SolutionHistory<Scalar>());
  }
  else {
    solutionHistory_ = sh;
  }
}

template <class Scalar>
void IntegratorBasic<Scalar>::setTimeStepControl(
    Teuchos::RCP<TimeStepControl<Scalar> > tsc)
{
  if (tsc == Teuchos::null) {
    // Create timeStepControl_ if null, otherwise keep current parameters.
    if (timeStepControl_ == Teuchos::null) {
      // Construct default TimeStepControl
      timeStepControl_ = rcp(new TimeStepControl<Scalar>());
    }
  }
  else {
    timeStepControl_ = tsc;
  }
}

template <class Scalar>
void IntegratorBasic<Scalar>::setObserver(
    Teuchos::RCP<IntegratorObserver<Scalar> > obs)
{
  if (obs == Teuchos::null)
    integratorObserver_ = Teuchos::rcp(new IntegratorObserverBasic<Scalar>());
  else
    integratorObserver_ = obs;
}

template <class Scalar>
void IntegratorBasic<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepper_ == Teuchos::null, std::logic_error,
      "Error - Need to set the Stepper, setStepper(), before calling "
      "IntegratorBasic::initialize()\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
      solutionHistory_->getNumStates() < 1, std::out_of_range,
      "Error - SolutionHistory requires at least one SolutionState.\n"
          << "        Supplied SolutionHistory has only "
          << solutionHistory_->getNumStates() << " SolutionStates.\n");

  stepper_->initialize();
  solutionHistory_->initialize();
  timeStepControl_->initialize();

  isInitialized_ = true;
}

template <class Scalar>
std::string IntegratorBasic<Scalar>::description() const
{
  std::string name = "Tempus::IntegratorBasic";
  return (name);
}

template <class Scalar>
void IntegratorBasic<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << "\n--- " << this->description() << " ---" << std::endl;

  if (solutionHistory_ != Teuchos::null) {
    solutionHistory_->describe(*l_out, verbLevel);
  }
  else {
    *l_out << "solutionHistory = " << solutionHistory_ << std::endl;
  }

  if (timeStepControl_ != Teuchos::null) {
    timeStepControl_->describe(out, verbLevel);
  }
  else {
    *l_out << "timeStepControl = " << timeStepControl_ << std::endl;
  }

  if (stepper_ != Teuchos::null) {
    stepper_->describe(out, verbLevel);
  }
  else {
    *l_out << "stepper         = " << stepper_ << std::endl;
  }
  *l_out << std::string(this->description().length() + 8, '-') << std::endl;
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
  out->setOutputToRootOnly(0);
  if (isInitialized_ == false) {
    Teuchos::OSTab ostab(out, 1, "StartIntegrator");
    *out << "Failure - IntegratorBasic is not initialized." << std::endl;
    setStatus(Status::FAILED);
    return;
  }

  // set the Abs/Rel tolerance
  auto cs = solutionHistory_->getCurrentState();
  cs->setTolRel(timeStepControl_->getMaxRelError());
  cs->setTolAbs(timeStepControl_->getMaxAbsError());

  integratorTimer_->start();
  // get optimal initial time step
  const Scalar initDt = std::min(timeStepControl_->getInitTimeStep(),
                                 stepper_->getInitTimeStep(solutionHistory_));
  // update initial time step
  timeStepControl_->setInitTimeStep(initDt);
  timeStepControl_->initialize();
  setStatus(WORKING);
}

template <class Scalar>
bool IntegratorBasic<Scalar>::advanceTime()
{
  TEMPUS_FUNC_TIME_MONITOR("Tempus::IntegratorBasic::advanceTime()");
  {
    startIntegrator();
    integratorObserver_->observeStartIntegrator(*this);

    while (
        integratorStatus_ == WORKING &&
        timeStepControl_->timeInRange(solutionHistory_->getCurrentTime()) &&
        timeStepControl_->indexInRange(solutionHistory_->getCurrentIndex())) {
      stepperTimer_->reset();
      stepperTimer_->start();
      solutionHistory_->initWorkingState();

      startTimeStep();
      integratorObserver_->observeStartTimeStep(*this);

      timeStepControl_->setNextTimeStep(solutionHistory_, integratorStatus_);
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

  // set the Abs/Rel tolerance
  ws->setTolRel(timeStepControl_->getMaxRelError());
  ws->setTolAbs(timeStepControl_->getMaxAbsError());

  // Check if we need to dump screen output this step
  std::vector<int>::const_iterator it = std::find(
      outputScreenIndices_.begin(), outputScreenIndices_.end(), ws->getIndex());
  if (it == outputScreenIndices_.end())
    ws->setOutputScreen(false);
  else
    ws->setOutputScreen(true);

  const int initial = timeStepControl_->getInitIndex();
  if ((ws->getIndex() - initial) % outputScreenInterval_ == 0)
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
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 2, "checkTimeStep");
    *out << "Failure - Stepper has failed more than the maximum allowed.\n"
         << "  (nFailures = " << ws->getNFailures()
         << ") >= (nFailuresMax = " << timeStepControl_->getMaxFailures() << ")"
         << std::endl;
    setStatus(Status::FAILED);
    return;
  }
  if (ws->getNConsecutiveFailures() >=
      timeStepControl_->getMaxConsecFailures()) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1, "checkTimeStep");
    *out << "Failure - Stepper has failed more than the maximum "
         << "consecutive allowed.\n"
         << "  (nConsecutiveFailures = " << ws->getNConsecutiveFailures()
         << ") >= (nConsecutiveFailuresMax = "
         << timeStepControl_->getMaxConsecFailures() << ")" << std::endl;
    setStatus(Status::FAILED);
    return;
  }

  // Timestep size is at the minimum timestep size and the step failed.
  if (ws->getTimeStep() <= timeStepControl_->getMinTimeStep() &&
      ws->getSolutionStatus() == Status::FAILED) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1, "checkTimeStep");
    *out << "Failure - Stepper has failed and the time step size is "
         << "at the minimum.\n"
         << "  Solution Status = " << toString(ws->getSolutionStatus())
         << std::endl
         << "  (TimeStep = " << ws->getTimeStep()
         << ") <= (Minimum TimeStep = " << timeStepControl_->getMinTimeStep()
         << ")" << std::endl;
    setStatus(Status::FAILED);
    return;
  }

  // Check Stepper failure.
  if (ws->getSolutionStatus() == Status::FAILED ||
      // Constant time step failure
      ((timeStepControl_->getStepType() == "Constant") &&
       !approxEqual(ws->getTimeStep(), timeStepControl_->getInitTimeStep()))) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 0, "checkTimeStep");
    *out << std::scientific << std::setw(6) << std::setprecision(3)
         << ws->getIndex() << std::setw(11) << std::setprecision(3)
         << ws->getTime() << std::setw(11) << std::setprecision(3)
         << ws->getTimeStep() << "  STEP FAILURE!! - ";
    if (ws->getSolutionStatus() == Status::FAILED) {
      *out << "Solution Status = " << toString(ws->getSolutionStatus())
           << std::endl;
    }
    else if ((timeStepControl_->getStepType() == "Constant") &&
             (ws->getTimeStep() != timeStepControl_->getInitTimeStep())) {
      *out << "dt != Constant dt (=" << timeStepControl_->getInitTimeStep()
           << ")" << std::endl;
    }

    ws->setNFailures(ws->getNFailures() + 1);
    ws->setNRunningFailures(ws->getNRunningFailures() + 1);
    ws->setNConsecutiveFailures(ws->getNConsecutiveFailures() + 1);
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
          Status::FAILED ||
      integratorStatus_ == Status::FAILED) {
    exitStatus = "Time integration FAILURE!";
  }
  else {
    setStatus(Status::PASSED);
    exitStatus = "Time integration complete.";
  }

  integratorTimer_->stop();
}

template <class Scalar>
void IntegratorBasic<Scalar>::setScreenOutputIndexList(std::string str)
{
  // Parse output indices
  std::string delimiters(",");
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
    std::string token = str.substr(lastPos, pos - lastPos);
    outputScreenIndices_.push_back(int(std::stoi(token)));
    if (pos == std::string::npos) break;

    lastPos = str.find_first_not_of(delimiters, pos);
    pos     = str.find_first_of(delimiters, lastPos);
  }

  // order output indices and remove duplicates.
  std::sort(outputScreenIndices_.begin(), outputScreenIndices_.end());
  outputScreenIndices_.erase(
      std::unique(outputScreenIndices_.begin(), outputScreenIndices_.end()),
      outputScreenIndices_.end());
  return;
}

template <class Scalar>
std::string IntegratorBasic<Scalar>::getScreenOutputIndexListString() const
{
  std::stringstream ss;
  for (size_t i = 0; i < outputScreenIndices_.size(); ++i) {
    if (i != 0) ss << ", ";
    ss << outputScreenIndices_[i];
  }
  return ss.str();
}

/** \brief Create valid IntegratorBasic ParameterList.
 */
template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
IntegratorBasic<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList(getIntegratorName());

  pl->set("Integrator Type", getIntegratorType(),
          "'Integrator Type' must be 'Integrator Basic'.");

  pl->set("Screen Output Index List", getScreenOutputIndexListString(),
          "Screen Output Index List.  Required to be in TimeStepControl range "
          "['Minimum Time Step Index', 'Maximum Time Step Index']");

  pl->set("Screen Output Index Interval", getScreenOutputIndexInterval(),
          "Screen Output Index Interval (e.g., every 100 time steps)");

  pl->set("Stepper Name", stepper_->getStepperName(),
          "'Stepper Name' selects the Stepper block to construct (Required).");

  pl->set("Solution History", *solutionHistory_->getValidParameters());
  pl->set("Time Step Control", *timeStepControl_->getValidParameters());

  Teuchos::RCP<Teuchos::ParameterList> tempusPL =
      Teuchos::parameterList("Tempus");

  tempusPL->set("Integrator Name", pl->name());
  tempusPL->set(pl->name(), *pl);
  tempusPL->set(stepper_->getStepperName(), *stepper_->getValidParameters());

  return tempusPL;
}

// Nonmember constructor
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList> tempusPL, bool runInitialize)
{
  auto integrator = Teuchos::rcp(new IntegratorBasic<Scalar>());
  if (tempusPL == Teuchos::null || tempusPL->numParams() == 0)
    return integrator;  // integrator is not initialized (missing model and IC).

  auto integratorName = tempusPL->get<std::string>("Integrator Name");
  auto integratorPL   = Teuchos::sublist(tempusPL, integratorName, true);

  std::string integratorType =
      integratorPL->get<std::string>("Integrator Type");
  TEUCHOS_TEST_FOR_EXCEPTION(
      integratorType != "Integrator Basic", std::logic_error,
      "Error - For IntegratorBasic, 'Integrator Type' should be "
          << "'Integrator Basic'.\n"
          << "    Integrator Type = " << integratorType << "\n");

  integrator->setIntegratorName(integratorName);

  // Validate the Integrator ParameterList
  auto validPL = Teuchos::rcp_const_cast<Teuchos::ParameterList>(
      integrator->getValidParameters());
  auto vIntegratorName = validPL->template get<std::string>("Integrator Name");
  auto vIntegratorPL   = Teuchos::sublist(validPL, vIntegratorName, true);
  integratorPL->validateParametersAndSetDefaults(*vIntegratorPL, 1);

  // Set Stepper
  if (integratorPL->isParameter("Stepper Name")) {
    // Construct from Integrator ParameterList
    auto stepperName = integratorPL->get<std::string>("Stepper Name");
    auto stepperPL   = Teuchos::sublist(tempusPL, stepperName, true);
    stepperPL->setName(stepperName);
    auto sf = Teuchos::rcp(new StepperFactory<Scalar>());
    integrator->setStepper(sf->createStepper(stepperPL));
  }
  else {
    // Construct default Stepper
    auto stepper = Teuchos::rcp(new StepperForwardEuler<Scalar>());
    integrator->setStepper(stepper);
  }

  // Set TimeStepControl
  if (integratorPL->isSublist("Time Step Control")) {
    // Construct from Integrator ParameterList
    auto tscPL = Teuchos::sublist(integratorPL, "Time Step Control", true);
    integrator->setTimeStepControl(
        createTimeStepControl<Scalar>(tscPL, runInitialize));
  }
  else {
    // Construct default TimeStepControl
    integrator->setTimeStepControl(rcp(new TimeStepControl<Scalar>()));
  }

  // Set SolutionHistory
  if (integratorPL->isSublist("Solution History")) {
    // Construct from Integrator ParameterList
    auto shPL = Teuchos::sublist(integratorPL, "Solution History", true);
    auto sh   = createSolutionHistoryPL<Scalar>(shPL);
    integrator->setSolutionHistory(sh);
  }
  else {
    // Construct default SolutionHistory
    integrator->setSolutionHistory(createSolutionHistory<Scalar>());
  }

  // Set Observer to default.
  integrator->setObserver(Teuchos::null);

  // Set screen output interval.
  integrator->setScreenOutputIndexInterval(
      integratorPL->get<int>("Screen Output Index Interval",
                             integrator->getScreenOutputIndexInterval()));

  // Parse screen output indices
  auto str = integratorPL->get<std::string>("Screen Output Index List", "");
  integrator->setScreenOutputIndexList(str);

  return integrator;  // integrator is not initialized (missing model and IC).
}

// Nonmember constructor
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList> tempusPL,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    bool runInitialize)
{
  auto integrator = createIntegratorBasic<Scalar>(tempusPL, runInitialize);
  if (model == Teuchos::null) return integrator;

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > constModel = model;
  integrator->setModel(constModel);

  // Construct default IC state from the application model and TimeStepControl
  auto newState =
      createSolutionStateME(integrator->getStepper()->getModel(),
                            integrator->getStepper()->getDefaultStepperState());
  newState->setTime(integrator->getTimeStepControl()->getInitTime());
  newState->setIndex(integrator->getTimeStepControl()->getInitIndex());
  newState->setTimeStep(integrator->getTimeStepControl()->getInitTimeStep());
  newState->setTolRel(integrator->getTimeStepControl()->getMaxRelError());
  newState->setTolAbs(integrator->getTimeStepControl()->getMaxAbsError());
  newState->setOrder(integrator->getStepper()->getOrder());
  newState->setSolutionStatus(Status::PASSED);  // ICs are considered passing.

  // Set SolutionHistory IC
  auto sh = integrator->getNonConstSolutionHistory();
  sh->addState(newState);
  integrator->getStepper()->setInitialConditions(sh);

  if (runInitialize) integrator->initialize();

  return integrator;
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType)
{
  auto integrator = Teuchos::rcp(new IntegratorBasic<Scalar>());

  auto sf      = Teuchos::rcp(new StepperFactory<Scalar>());
  auto stepper = sf->createStepper(stepperType, model);
  integrator->setStepper(stepper);
  integrator->initializeSolutionHistory();
  integrator->initialize();

  return integrator;
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic()
{
  return Teuchos::rcp(new IntegratorBasic<Scalar>());
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorBasic<Scalar> > createIntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList> tempusPL,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models,
    bool runInitialize)
{
  auto integratorName = tempusPL->get<std::string>("Integrator Name");
  auto integratorPL   = Teuchos::sublist(tempusPL, integratorName, true);

  std::string integratorType =
      integratorPL->get<std::string>("Integrator Type");
  TEUCHOS_TEST_FOR_EXCEPTION(
      integratorType != "Integrator Basic", std::logic_error,
      "Error - For IntegratorBasic, 'Integrator Type' should be "
          << "'Integrator Basic'.\n"
          << "    Integrator Type = " << integratorType << "\n");

  auto integrator = Teuchos::rcp(new IntegratorBasic<Scalar>());
  integrator->setIntegratorName(integratorName);

  TEUCHOS_TEST_FOR_EXCEPTION(
      !integratorPL->isParameter("Stepper Name"), std::logic_error,
      "Error - Need to set the 'Stepper Name' in 'Integrator Basic'.\n");

  auto stepperName = integratorPL->get<std::string>("Stepper Name");
  TEUCHOS_TEST_FOR_EXCEPTION(
      stepperName == "Operator Split", std::logic_error,
      "Error - 'Stepper Name' should be 'Operator Split'.\n");

  // Construct Steppers from Integrator ParameterList
  auto stepperPL = Teuchos::sublist(tempusPL, stepperName, true);
  stepperPL->setName(stepperName);
  auto sf = Teuchos::rcp(new StepperFactory<Scalar>());
  integrator->setStepper(sf->createStepper(stepperPL, models));

  // Set TimeStepControl
  if (integratorPL->isSublist("Time Step Control")) {
    // Construct from Integrator ParameterList
    auto tscPL = Teuchos::sublist(integratorPL, "Time Step Control", true);
    integrator->setTimeStepControl(
        createTimeStepControl<Scalar>(tscPL, runInitialize));
  }
  else {
    // Construct default TimeStepControl
    integrator->setTimeStepControl(rcp(new TimeStepControl<Scalar>()));
  }

  // Construct default IC state from the application model and TimeStepControl
  auto newState =
      createSolutionStateME(integrator->getStepper()->getModel(),
                            integrator->getStepper()->getDefaultStepperState());
  newState->setTime(integrator->getTimeStepControl()->getInitTime());
  newState->setIndex(integrator->getTimeStepControl()->getInitIndex());
  newState->setTimeStep(integrator->getTimeStepControl()->getInitTimeStep());
  newState->setTolRel(integrator->getTimeStepControl()->getMaxRelError());
  newState->setTolAbs(integrator->getTimeStepControl()->getMaxAbsError());
  newState->setOrder(integrator->getStepper()->getOrder());
  newState->setSolutionStatus(Status::PASSED);  // ICs are considered passing.

  // Set SolutionHistory
  auto shPL = Teuchos::sublist(integratorPL, "Solution History", true);
  auto sh   = createSolutionHistoryPL<Scalar>(shPL);
  sh->addState(newState);
  integrator->getStepper()->setInitialConditions(sh);
  integrator->setSolutionHistory(sh);

  // Set Observer to default.
  integrator->setObserver(Teuchos::null);

  // Set screen output interval.
  integrator->setScreenOutputIndexInterval(
      integratorPL->get<int>("Screen Output Index Interval",
                             integrator->getScreenOutputIndexInterval()));

  // Parse screen output indices
  auto str = integratorPL->get<std::string>("Screen Output Index List", "");
  integrator->setScreenOutputIndexList(str);

  auto validPL = Teuchos::rcp_const_cast<Teuchos::ParameterList>(
      integrator->getValidParameters());

  // Validate the Integrator ParameterList
  auto vIntegratorName = validPL->template get<std::string>("Integrator Name");
  auto vIntegratorPL   = Teuchos::sublist(validPL, vIntegratorName, true);
  integratorPL->validateParametersAndSetDefaults(*vIntegratorPL);

  // Validate the Stepper ParameterList
  auto vStepperName = vIntegratorPL->template get<std::string>("Stepper Name");
  auto vStepperPL   = Teuchos::sublist(validPL, vStepperName, true);
  stepperPL->validateParametersAndSetDefaults(*vStepperPL);

  integrator->initialize();

  return integrator;
}

}  // namespace Tempus
#endif  // Tempus_IntegratorBasic_impl_hpp
