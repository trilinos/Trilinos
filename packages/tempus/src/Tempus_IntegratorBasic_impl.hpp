#ifndef Tempus_IntegratorBasic_impl_hpp
#define Tempus_IntegratorBasic_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <ctime>


namespace {

  static std::string integratorName_name    = "Integrator Name";
  static std::string StepperName_name       = "Stepper Name";

  static std::string integratorType_name    = "Integrator Type";
  static std::string integratorType_default = "Integrator Basic";

  static std::string initTime_name          = "Initial Time";
  static double      initTime_default       = 0.0;
  static std::string initTimeIndex_name     = "Initial Time Index";
  static int         initTimeIndex_default  = 0;
  static std::string initTimeStep_name      = "Initial Time Step";
  static double      initTimeStep_default   = std::numeric_limits<double>::epsilon();
  static std::string initOrder_name         = "Initial Order";
  static int         initOrder_default      = 0;
  static std::string finalTime_name         = "Final Time";
  static double      finalTime_default      = std::numeric_limits<double>::max();
  static std::string finalTimeIndex_name    = "Final Time Index";
  static int         finalTimeIndex_default = std::numeric_limits<int>::max();

  static std::string outputScreenTimeList_name     = "Screen Output Time List";
  static std::string outputScreenTimeList_default  = "";
  static std::string outputScreenIndexList_name    = "Screen Output Index List";
  static std::string outputScreenIndexList_default = "";
  static std::string outputScreenTimeInterval_name     =
    "Screen Output Time Interval";
  //static double      outputScreenTimeInterval_default  = 100.0;
  static std::string outputScreenIndexInterval_name    =
    "Screen Output Index Interval";
  static int         outputScreenIndexInterval_default = 100;

  static std::string SolutionHistory_name     = "Solution History";
  static std::string TimeStepControl_name     = "Time Step Control";

} // namespace


namespace Tempus {

template<class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic(
  Teuchos::RCP<Teuchos::ParameterList>                inputPL,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model)
     : tempusPL_(inputPL), integratorStatus_(WORKING)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Get the name of the integrator to build.
  std::string integratorName_ =tempusPL_->get<std::string>(integratorName_name);

  this->setParameterList(Teuchos::sublist(tempusPL_, integratorName_, true));
  this->setStepper(model);
  this->initialize();

  if (Teuchos::as<int>(this->getVerbLevel()) >=
      Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"IntegratorBasic::IntegratorBasic");
    *out << this->description() << std::endl;
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::setStepper(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  if (stepper_ == Teuchos::null) {
    // Construct from Integrator ParameterList
    RCP<StepperFactory<Scalar> > sf =Teuchos::rcp(new StepperFactory<Scalar>());
    std::string stepperName = integratorPL_->get<std::string>(StepperName_name);
    RCP<ParameterList> stepperPL = Teuchos::sublist(tempusPL_,stepperName,true);
    stepper_ = sf->createStepper(stepperPL, model);
  } else {
    stepper_->setModel(model);
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::setStepper(
  Teuchos::RCP<Stepper<Scalar> > newStepper)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Make integratorPL_ consistent with new stepper.
  RCP<ParameterList> newStepperPL = newStepper->getNonconstParameterList();
  integratorPL_->set(StepperName_name, newStepperPL->name());
  integratorPL_->set(newStepperPL->name(), newStepperPL);

  stepper_ = Teuchos::null;
  stepper_ = newStepper;
}


template<class Scalar>
void IntegratorBasic<Scalar>::setSolutionHistory(
  Teuchos::RCP<SolutionHistory<Scalar> > sh)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  if (sh == Teuchos::null) {
    // Construct from Integrator ParameterList
    RCP<ParameterList> shPL =
      Teuchos::sublist(integratorPL_, SolutionHistory_name, true);
    solutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));

    // Create IC SolutionState
    // Create meta data
    RCP<SolutionStateMetaData<Scalar> > md =
                                     rcp(new SolutionStateMetaData<Scalar> ());
    md->setTime (
      integratorPL_->get<double>(initTime_name,     initTime_default));
    md->setIStep(
      integratorPL_->get<int>   (initTimeIndex_name,initTimeIndex_default));
    md->setDt   (
      integratorPL_->get<double>(initTimeStep_name, initTimeStep_default));
    int orderTmp = integratorPL_->get<int> (initOrder_name, initOrder_default);
    if (orderTmp == 0) orderTmp = stepper_->getOrderMin();
    md->setOrder(orderTmp);
    md->setSolutionStatus(Status::PASSED);  // ICs are considered passing.

    // Create initial condition solution state
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
    RCP<Thyra::VectorBase<Scalar> > xdotdot = Teuchos::null;
    RCP<SolutionState<Scalar> > newState = rcp(new SolutionState<Scalar>(
      md, x, xdot, xdotdot, stepper_->getDefaultStepperState()));

    solutionHistory_->addState(newState);

  } else {

    TEUCHOS_TEST_FOR_EXCEPTION( sh->getNumStates() < 1,
      std::out_of_range,
         "Error - setSolutionHistory requires at least one SolutionState.\n"
      << "        Supplied SolutionHistory has only " << sh->getNumStates()
      << " SolutionStates.\n");

    // Make integratorPL_ consistent with new SolutionHistory.
    RCP<ParameterList> shPL = sh->getNonconstParameterList();
    integratorPL_->set(SolutionHistory_name, shPL->name());
    integratorPL_->set(shPL->name(), shPL);

    solutionHistory_ = Teuchos::null;
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
      Teuchos::sublist(integratorPL_,TimeStepControl_name,true);
    Scalar dtTmp =
      integratorPL_->get<double>(initTimeStep_name, initTimeStep_default);
    timeStepControl_ = rcp(new TimeStepControl<Scalar>(tscPL, dtTmp));

    if (timeStepControl_->orderMin_ == 0)
      timeStepControl_->orderMin_ = stepper_->getOrderMin();
    if (timeStepControl_->orderMax_ == 0)
      timeStepControl_->orderMax_ = stepper_->getOrderMax();

  } else {
    // Make integratorPL_ consistent with new TimeStepControl.
    RCP<ParameterList> tscPL = tsc->getNonconstParameterList();
    integratorPL_->set(SolutionHistory_name, tscPL->name());
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
    // Create default IntegratorObserver
    integratorObserver_ =
      Teuchos::rcp(new IntegratorObserver<Scalar>(solutionHistory_,
                                                  timeStepControl_));
  } else {
    integratorObserver_ = obs;
  }
}


template<class Scalar>
void IntegratorBasic<Scalar>::initialize()
{
  // Check that ParameterList is valid.
  integratorPL_->validateParameters(*this->getValidParameters());
  this->setSolutionHistory();
  this->setTimeStepControl();
  this->setObserver();

  // Check values and ranges of Parameters.
  Scalar initTime = integratorPL_->get<double>(initTime_name, initTime_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (initTime<timeStepControl_->timeMin_|| initTime>timeStepControl_->timeMax_),
    std::out_of_range,
    "Error - Initial time is out of range.\n"
    << "    [timeMin, timeMax] = [" << timeStepControl_->timeMin_ << ", "
                                    << timeStepControl_->timeMax_ << "]\n"
    << "    initTime = " << initTime << "\n");

  int iStep =
    integratorPL_->get<int>(initTimeIndex_name, initTimeIndex_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (iStep < timeStepControl_->iStepMin_|| iStep > timeStepControl_->iStepMax_),
    std::out_of_range,
    "Error - Initial time index is out of range.\n"
    << "    [iStepMin, iStepMax] = [" << timeStepControl_->iStepMin_ << ", "
                                      << timeStepControl_->iStepMax_ << "]\n"
    << "    iStep = " << iStep << "\n");

  Scalar dt =
    integratorPL_->get<double>(initTimeStep_name, initTimeStep_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dt < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative time step.  dt = "<<dt<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dt < timeStepControl_->dtMin_ || dt > timeStepControl_->dtMax_ ),
    std::out_of_range,
    "Error - Initial time step is out of range.\n"
    << "    [dtMin, dtMax] = [" << timeStepControl_->dtMin_ << ", "
                                << timeStepControl_->dtMax_ << "]\n"
    << "    dt = " << dt << "\n");

  int order = integratorPL_->get<int>(initOrder_name, initOrder_default);
  if (order == 0) order = stepper_->getOrderMin();
  TEUCHOS_TEST_FOR_EXCEPTION(
    (order<timeStepControl_->orderMin_ || order>timeStepControl_->orderMax_),
    std::out_of_range,
    "Error - Initial order is out of range.\n"
    << "    [orderMin, orderMax] = [" << timeStepControl_->orderMin_ << ", "
                                      << timeStepControl_->orderMax_ << "]\n"
    << "    order = " << order << "\n");

  Scalar finalTime =
    integratorPL_->get<double>(finalTime_name, finalTime_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (finalTime < timeStepControl_->timeMin_ ||
     finalTime > timeStepControl_->timeMax_),
    std::out_of_range,
    "Error - Final time is out of range.\n"
    << "    [timeMin, timeMax] = [" << timeStepControl_->timeMin_ << ", "
                                    << timeStepControl_->timeMax_ << "]\n"
    << "    finalTime = " << finalTime << "\n");

  int fiStep =
    integratorPL_->get<int>(finalTimeIndex_name, finalTimeIndex_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (fiStep<timeStepControl_->iStepMin_ || fiStep>timeStepControl_->iStepMax_),
    std::out_of_range,
    "Error - Final time index is out of range.\n"
    << "    [iStepMin, iStepMax] = [" << timeStepControl_->iStepMin_ << ", "
                                      << timeStepControl_->iStepMax_ << "]\n"
    << "    iStep = " << fiStep << "\n");


  // Parse output indices
  {
    outputScreenIndices_.clear();
    std::string str =
      integratorPL_->get<std::string>(outputScreenIndexList_name,
                                              outputScreenIndexList_default);
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

    Scalar outputScreenIndexInterval =
      integratorPL_->get<int>(outputScreenIndexInterval_name,
                              outputScreenIndexInterval_default);
    Scalar outputScreen_i = timeStepControl_->iStepMin_;
    while (outputScreen_i <= timeStepControl_->iStepMax_) {
      outputScreenIndices_.push_back(outputScreen_i);
      outputScreen_i += outputScreenIndexInterval;
    }

    // order output indices
    std::sort(outputScreenIndices_.begin(),outputScreenIndices_.end());
  }

  if (integratorTimer_ == Teuchos::null)
    integratorTimer_ = rcp(new Teuchos::Time("Integrator Timer"));
  if (stepperTimer_ == Teuchos::null)
    stepperTimer_    = rcp(new Teuchos::Time("Stepper Timer"));
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
    timeStepControl_->timeMax_ = timeFinal;
  bool itgrStatus = advanceTime();
  return itgrStatus;
}


template <class Scalar>
void IntegratorBasic<Scalar>::startIntegrator()
{
  std::time_t begin = std::time(nullptr);
  integratorTimer_->start();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,0,"ScreenOutput");
  *out << "\nTempus - IntegratorBasic\n"
       << std::asctime(std::localtime(&begin)) << "\n"
       << "  Stepper = " << stepper_->description() << "\n"
       << "  Simulation Time Range  [" << timeStepControl_->timeMin_
       << ", " << timeStepControl_->timeMax_ << "]\n"
       << "  Simulation Index Range [" << timeStepControl_->iStepMin_
       << ", " << timeStepControl_->iStepMax_ << "]\n"
       << "============================================================================\n"
       << "  Step       Time         dt  Abs Error  Rel Error  Order  nFail  dCompTime"
       << std::endl;
  integratorStatus_ = WORKING;
}


template <class Scalar>
bool IntegratorBasic<Scalar>::advanceTime()
{
  startIntegrator();
  integratorObserver_->observeStartIntegrator();

  while (integratorStatus_ == WORKING and
         timeStepControl_->timeInRange (solutionHistory_->getCurrentTime()) and
         timeStepControl_->indexInRange(solutionHistory_->getCurrentIndex())){

    stepperTimer_->reset();
    stepperTimer_->start();
    solutionHistory_->initWorkingState();

    startTimeStep();
    integratorObserver_->observeStartTimeStep();

    timeStepControl_->getNextTimeStep(solutionHistory_, integratorStatus_);
    integratorObserver_->observeNextTimeStep(integratorStatus_);

    if (integratorStatus_ == FAILED) break;

    integratorObserver_->observeBeforeTakeStep();

    stepper_->takeStep(solutionHistory_);

    integratorObserver_->observeAfterTakeStep();

    stepperTimer_->stop();
    acceptTimeStep();
    integratorObserver_->observeAcceptedTimeStep(integratorStatus_);
  }

  endIntegrator();
  integratorObserver_->observeEndIntegrator(integratorStatus_);

  return (integratorStatus_ == Status::PASSED);
}


template <class Scalar>
void IntegratorBasic<Scalar>::startTimeStep()
{
  Teuchos::RCP<SolutionStateMetaData<Scalar> > wsmd =
    solutionHistory_->getWorkingState()->metaData_;

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
    solutionHistory_->getWorkingState()->metaData_;

       // Stepper failure
  if ( solutionHistory_->getWorkingState()->getSolutionStatus() == FAILED or
       solutionHistory_->getWorkingState()->getStepperStatus() == FAILED or
       // Constant time step failure
       ((timeStepControl_->stepType_ == CONSTANT_STEP_SIZE) and
       (wsmd->getDt() != timeStepControl_->dtConstant_))
     )
  {
    wsmd->setNFailures(wsmd->getNFailures()+1);
    wsmd->setNConsecutiveFailures(wsmd->getNConsecutiveFailures()+1);
    wsmd->setSolutionStatus(FAILED);
  }

  // Too many failures
  if (wsmd->getNFailures() >= timeStepControl_->nFailuresMax_) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"continueIntegration");
    *out << "Failure - Stepper has failed more than the maximum allowed.\n"
         << "  (nFailures = "<<wsmd->getNFailures()<< ") >= (nFailuresMax = "
         <<timeStepControl_->nFailuresMax_<<")" << std::endl;
    integratorStatus_ = FAILED;
    return;
  }
  if (wsmd->getNConsecutiveFailures()
      >= timeStepControl_->nConsecutiveFailuresMax_){
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"continueIntegration");
    *out << "Failure - Stepper has failed more than the maximum "
         << "consecutive allowed.\n"
         << "  (nConsecutiveFailures = "<<wsmd->getNConsecutiveFailures()
         << ") >= (nConsecutiveFailuresMax = "
         <<timeStepControl_->nConsecutiveFailuresMax_
         << ")" << std::endl;
    integratorStatus_ = FAILED;
    return;
  }

  // =======================================================================
  // Made it here! Accept this time step

  solutionHistory_->promoteWorkingState();

  RCP<SolutionStateMetaData<Scalar> > csmd =
    solutionHistory_->getCurrentState()->metaData_;

  csmd->setNFailures(std::max(csmd->getNFailures()-1,0));
  csmd->setNConsecutiveFailures(0);

  // Output and screen output
  if (csmd->getOutput() == true) {
    // Dump solution!
  }

  if (csmd->getOutputScreen() == true) {
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
void IntegratorBasic<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & tmpPL)
{
  if (tmpPL == Teuchos::null)
   integratorPL_->validateParametersAndSetDefaults(*this->getValidParameters());
  else {
    integratorPL_ = tmpPL;
    integratorPL_->validateParameters(*this->getValidParameters());
  }

  Teuchos::readVerboseObjectSublist(&*integratorPL_,this);

  std::string integratorType =
    integratorPL_->get<std::string>(integratorType_name);
  TEUCHOS_TEST_FOR_EXCEPTION( integratorType != integratorType_default,
    std::logic_error,
    "Error - Inconsistent Integrator Type for IntegratorBasic\n"
    << "    " << integratorType_name << " = " << integratorType << "\n");

  return;
}


/** \brief Create valid IntegratorBasic ParameterList.
 */
template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
IntegratorBasic<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    std::ostringstream tmp;
    tmp << "'" << integratorType_name << "' must be '" << integratorType_default << "'.";
    pl->set(integratorType_name, integratorType_default, tmp.str());

    tmp.clear();
    tmp << "Initial time.  Required to be in TimeStepControl range "
        << "['Minimum Simulation Time', 'Maximum Simulation Time']";
    pl->set(initTime_name, initTime_default, tmp.str());

    tmp << "Initial time index.  Required to be in TimeStepControl range "
        << "['Minimum Time Step Index', 'Maximum Time Step Index']";
    pl->set(initTimeIndex_name, initTimeIndex_default, tmp.str());

    tmp.clear();
    tmp << "Initial time step.  Required to be positive and in TimeStepControl "
        << "range ['Minimum Time Step', 'Maximum Time Step']";
    pl->set(initTimeStep_name, initTimeStep_default, tmp.str());

    tmp.clear();
    tmp << "Initial order.  Required to be in TimeStepControl range "
        << "['Minimum Time Integration Order',"
        << " 'Maximum Time Integration Order']"
        << "  If set to zero (default), it "
        << "will be reset to the Stepper minimum order.";
    pl->set(initOrder_name, initOrder_default, tmp.str());

    tmp.clear();
    tmp << "Final time.  Required to be in the TimeStepControl range "
        << "['Minimum Simulation Time', 'Maximum Simulation Time']";
    pl->set(finalTime_name, finalTime_default, tmp.str());

    tmp.clear();
    tmp << "Final time index.  Required to be in TimeStepControl range "
        << "['Minimum Time Step Index', 'Maximum Time Step Index']";
    pl->set(finalTimeIndex_name, finalTimeIndex_default, tmp.str());

    tmp.clear();
    tmp << "Screen Output Index List.  Required to be in TimeStepControl range "
        << "['Minimum Time Step Index', 'Maximum Time Step Index']";
    pl->set(outputScreenIndexList_name, outputScreenIndexList_default,
            tmp.str());
    pl->set(outputScreenIndexInterval_name, outputScreenIndexInterval_default,
      "Screen Output Index Interval (e.g., every 100 time steps");

    pl->set( StepperName_name, "",
      "'Stepper Name' selects the Stepper block to construct (Required).");

    // Solution History
    pl->sublist(SolutionHistory_name,false,"solutionHistory_docs")
        .disableRecursiveValidation();

    // Time Step Control
    pl->sublist(TimeStepControl_name,false,"solutionHistory_docs")
        .disableRecursiveValidation();

    validPL = pl;
  }
  return validPL;
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

} // namespace Tempus
#endif // Tempus_IntegratorBasic_impl_hpp
