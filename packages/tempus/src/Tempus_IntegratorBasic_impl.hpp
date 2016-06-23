#ifndef TEMPUS_INTEGRATORBASIC_IMPL_HPP
#define TEMPUS_INTEGRATORBASIC_IMPL_HPP

// Teuchos
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
// Tempus
#include "Tempus_StepperFactory.hpp"
#include "Tempus_TimeStepControl.hpp"

#include <ctime>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;

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
  static int         initOrder_default      = 1;
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
  static double      outputScreenTimeInterval_default  = 100.0;
  static std::string outputScreenIndexInterval_name    =
    "Screen Output Index Interval";
  static int         outputScreenIndexInterval_default = 100;

  static std::string SolutionHistory_name     = "Solution History";
  static std::string TimeStepControl_name     = "Time Step Control";

} // namespace


namespace Tempus {

template<class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic(
  RCP<ParameterList>                         tempusPL,
  const RCP<Thyra::ModelEvaluator<Scalar> >& model,
  const RCP<IntegratorObserver<Scalar> >&    integratorObserver)
     : integratorStatus (WORKING)
{
  // Get the name of the integrator to build.
  integratorName_ = tempusPL->get<std::string>(integratorName_name);

  // Create classes from nested blocks prior to setParameters call.
  RCP<ParameterList> tmpiPL = Teuchos::sublist(tempusPL, integratorName_, true);

  //    Create solution history
  RCP<ParameterList> shPL = Teuchos::sublist(tmpiPL, SolutionHistory_name,true);
  solutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));

  //    Create TimeStepControl
  RCP<ParameterList> tscPL = Teuchos::sublist(tmpiPL,TimeStepControl_name,true);
  Scalar dtTmp = tmpiPL->get<double>(initTimeStep_name, initTimeStep_default);
  timeStepControl_ = rcp(new TimeStepControl<Scalar>(tscPL, dtTmp));

  this->setParameterList(tempusPL);

  // Create Stepper
  RCP<StepperFactory<Scalar> > sf = rcp(new StepperFactory<Scalar>());
  std::string stepperName = pList_->get<std::string>(StepperName_name);
  RCP<ParameterList> s_pl = Teuchos::sublist(tempusPL, stepperName, true);
  stepper_ = sf->createStepper(s_pl, model);

  // Create meta data
  RCP<SolutionStateMetaData<Scalar> > md =
                                   rcp(new SolutionStateMetaData<Scalar> ());
  md->time  = pList_->get<double>(initTime_name,      initTime_default);
  md->iStep = pList_->get<int>   (initTimeIndex_name, initTimeIndex_default);
  md->dt    = pList_->get<double>(initTimeStep_name,  initTimeStep_default);
  md->order = pList_->get<int>   (initOrder_name,     initOrder_default);

  // Create initial condition solution state
  typedef Thyra::ModelEvaluatorBase MEB;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgsIC =model->getNominalValues();
  RCP<Thyra::VectorBase<Scalar> > x = inArgsIC.get_x()->clone_v();
  RCP<Thyra::VectorBase<Scalar> > xdot;
  if (inArgsIC.supports(MEB::IN_ARG_x_dot))
    xdot = inArgsIC.get_x_dot()->clone_v();
  else
    xdot = x->clone_v();
  RCP<Thyra::VectorBase<Scalar> > xdotdot = Teuchos::null;
  RCP<SolutionState<Scalar> > newState = rcp(new SolutionState<Scalar>(
    md, x, xdot, xdotdot, stepper_->getDefaultStepperState()));
  solutionHistory_->addState(newState);

  if (integratorObserver == Teuchos::null) {
    // Create default IntegratorObserver
    integratorObserver_ = rcp(new IntegratorObserver<Scalar>(solutionHistory_,
                                                             timeStepControl_));
  } else {
    integratorObserver_ = integratorObserver;
  }

  integratorTimer = rcp(new Teuchos::Time("Integrator Timer"));
  stepperTimer    = rcp(new Teuchos::Time("Stepper Timer"));

  if (Teuchos::as<int>(this->getVerbLevel()) >=
      Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"IntegratorBasic::IntegratorBasic");
    *out << this->description() << std::endl;
  }
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
    timeStepControl_->timeMax = timeFinal;
  bool itgrStatus = advanceTime();
  return itgrStatus;
}


template <class Scalar>
void IntegratorBasic<Scalar>::startIntegrator()
{
  std::time_t begin = std::time(nullptr);
  integratorTimer->start();
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,0,"ScreenOutput");
  *out << "\nTempus - IntegratorBasic\n"
       << std::asctime(std::localtime(&begin)) << "\n"
       << "  Stepper = " << stepper_->description() << "\n"
       << "  Simulation Time Range  [" << timeStepControl_->timeMin
       << ", " << timeStepControl_->timeMax << "]\n"
       << "  Simulation Index Range [" << timeStepControl_->iStepMin
       << ", " << timeStepControl_->iStepMax << "]\n"
       << "============================================================================\n"
       << "  Step       Time         dt  Abs Error  Rel Error  Order  nFail  dCompTime"
       << std::endl;
  integratorStatus = WORKING;
}


template <class Scalar>
bool IntegratorBasic<Scalar>::advanceTime()
{
  startIntegrator();
  integratorObserver_->observeStartIntegrator();

  while (integratorStatus == WORKING and
         timeStepControl_->timeInRange (solutionHistory_->getCurrentTime()) and
         timeStepControl_->indexInRange(solutionHistory_->getCurrentIndex())){

    stepperTimer->reset();
    stepperTimer->start();
    integratorObserver_->observeStartTimeStep();

    solutionHistory_->initWorkingState();

    timeStepControl_->getNextTimeStep(solutionHistory_, integratorStatus);

    integratorObserver_->observeNextTimeStep(integratorStatus);

    if (integratorStatus == FAILED) break;

    integratorObserver_->observeBeforeTakeStep();

    stepper_->takeStep(solutionHistory_);

    integratorObserver_->observeAfterTakeStep();

    stepperTimer->stop();
    acceptTimeStep();
    integratorObserver_->observeAcceptedTimeStep(integratorStatus);
  }

  endIntegrator();
  integratorObserver_->observeEndIntegrator(integratorStatus);

  return (integratorStatus == Status::PASSED);
}


template <class Scalar>
void IntegratorBasic<Scalar>::acceptTimeStep()
{
  RCP<SolutionStateMetaData<Scalar> > wsmd =
    solutionHistory_->getWorkingState()->metaData;

       // Stepper failure
  if ( solutionHistory_->getWorkingState()->getSolutionStatus() == FAILED or
       solutionHistory_->getWorkingState()->getStepperStatus() == FAILED or
       // Constant time step failure
       ((timeStepControl_->stepType == CONSTANT_STEP_SIZE) and
       (wsmd->dt != timeStepControl_->dtConstant))
     )
  {
    wsmd->nFailures++;
    wsmd->nConsecutiveFailures++;
    wsmd->solutionStatus = FAILED;
  }

  // Too many failures
  if (wsmd->nFailures >= timeStepControl_->nFailuresMax) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"continueIntegration");
    *out << "Failure - Stepper has failed more than the maximum allowed.\n"
         << "  (nFailures = "<<wsmd->nFailures << ") >= (nFailuresMax = "
         <<timeStepControl_->nFailuresMax<<")" << std::endl;
    integratorStatus = FAILED;
    return;
  }
  if (wsmd->nConsecutiveFailures >= timeStepControl_->nConsecutiveFailuresMax) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"continueIntegration");
    *out << "Failure - Stepper has failed more than the maximum "
         << "consecutive allowed.\n"
         << "  (nConsecutiveFailures = "<<wsmd->nConsecutiveFailures
         << ") >= (nConsecutiveFailuresMax = "
         <<timeStepControl_->nConsecutiveFailuresMax
         << ")" << std::endl;
    integratorStatus = FAILED;
    return;
  }

  // =======================================================================
  // Made it here! Accept this time step

  solutionHistory_->promoteWorkingState();

  RCP<SolutionStateMetaData<Scalar> > csmd =
    solutionHistory_->getCurrentState()->metaData;

  csmd->nFailures = std::max(csmd->nFailures-1,0);
  csmd->nConsecutiveFailures = 0;

  // Output and screen output
  std::vector<int>::const_iterator it;
  it = std::find(outputScreenIndices.begin(),
                 outputScreenIndices.end(), csmd->iStep);
  if (it != outputScreenIndices.end()) {
    const double steppertime = stepperTimer->totalElapsedTime();
    stepperTimer->reset();
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,0,"ScreenOutput");
    *out <<std::scientific<<std::setw( 6)<<std::setprecision(3)<< csmd->iStep
         <<std::scientific<<std::setw(11)<<std::setprecision(3)<< csmd->time
         <<std::scientific<<std::setw(11)<<std::setprecision(3)<< csmd->dt
         <<std::scientific<<std::setw(11)<<std::setprecision(3)<< csmd->errorAbs
         <<std::scientific<<std::setw(11)<<std::setprecision(3)<< csmd->errorRel
         <<std::scientific<<std::setw( 7)<<std::setprecision(3)<< csmd->order
         <<std::scientific<<std::setw( 7)<<std::setprecision(3)<< csmd->nFailures
         <<std::scientific<<std::setw(11)<<std::setprecision(3)<< steppertime
         <<std::endl;
  }

  if (csmd->output == true) {
    // Dump solution!
  }
}


template <class Scalar>
void IntegratorBasic<Scalar>::endIntegrator()
{
  std::string exitStatus;
  if (solutionHistory_->getCurrentState()->getSolutionStatus() ==
      Status::FAILED or integratorStatus == Status::FAILED) {
    exitStatus = "Time integration FAILURE!";
  } else {
    integratorStatus = Status::PASSED;
    exitStatus = "Time integration complete.";
  }

  integratorTimer->stop();
  const double runtime = integratorTimer->totalElapsedTime();
  std::time_t end = std::time(nullptr);
  RCP<Teuchos::FancyOStream> out = this->getOStream();
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
  const RCP<ParameterList> & tempusPL)
{
  if (tempusPL == Teuchos::null)
    pList_->validateParametersAndSetDefaults(*this->getValidParameters());
  else {
    pList_ = Teuchos::sublist(tempusPL, integratorName_, true);
    pList_->validateParameters(*this->getValidParameters());
  }

  Teuchos::readVerboseObjectSublist(&*pList_,this);

  std::string integratorType = pList_->get<std::string>(integratorType_name);
  TEUCHOS_TEST_FOR_EXCEPTION( integratorType != integratorType_default,
    std::logic_error,
    "Error - Inconsistent Integrator Type for IntegratorBasic\n"
    << "    " << integratorType_name << " = " << integratorType << "\n");

  Scalar initTime = pList_->get<double>(initTime_name, initTime_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (initTime<timeStepControl_->timeMin || initTime>timeStepControl_->timeMax),
    std::out_of_range,
    "Error - Initial time is out of range.\n"
    << "    [timeMin, timeMax] = [" << timeStepControl_->timeMin << ", "
                                    << timeStepControl_->timeMax << "]\n"
    << "    initTime = " << initTime << "\n");

  int iStep = pList_->get<int>(initTimeIndex_name, initTimeIndex_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (iStep < timeStepControl_->iStepMin || iStep > timeStepControl_->iStepMax),
    std::out_of_range,
    "Error - Initial time index is out of range.\n"
    << "    [iStepMin, iStepMax] = [" << timeStepControl_->iStepMin << ", "
                                      << timeStepControl_->iStepMax << "]\n"
    << "    iStep = " << iStep << "\n");

  Scalar dt = pList_->get<double>(initTimeStep_name, initTimeStep_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dt < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative time step.  dt = "<<dt<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dt < timeStepControl_->dtMin || dt > timeStepControl_->dtMax ),
    std::out_of_range,
    "Error - Initial time step is out of range.\n"
    << "    [dtMin, dtMax] = [" << timeStepControl_->dtMin << ", "
                                << timeStepControl_->dtMax << "]\n"
    << "    dt = " << dt << "\n");

  int order = pList_->get<int>(initOrder_name, initOrder_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (order < timeStepControl_->orderMin || order > timeStepControl_->orderMax),
    std::out_of_range,
    "Error - Initial order is out of range.\n"
    << "    [orderMin, orderMax] = [" << timeStepControl_->orderMin << ", "
                                      << timeStepControl_->orderMax << "]\n"
    << "    order = " << order << "\n");

  Scalar finalTime = pList_->get<double>(finalTime_name, finalTime_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (finalTime<timeStepControl_->timeMin || finalTime>timeStepControl_->timeMax),
    std::out_of_range,
    "Error - Final time is out of range.\n"
    << "    [timeMin, timeMax] = [" << timeStepControl_->timeMin << ", "
                                    << timeStepControl_->timeMax << "]\n"
    << "    finalTime = " << finalTime << "\n");

  int fiStep = pList_->get<int>(finalTimeIndex_name, finalTimeIndex_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (fiStep < timeStepControl_->iStepMin || fiStep > timeStepControl_->iStepMax),
    std::out_of_range,
    "Error - Final time index is out of range.\n"
    << "    [iStepMin, iStepMax] = [" << timeStepControl_->iStepMin << ", "
                                      << timeStepControl_->iStepMax << "]\n"
    << "    iStep = " << fiStep << "\n");

  // Parse output indices
  {
    outputScreenIndices.clear();
    std::string str = pList_->get<std::string>(outputScreenIndexList_name,
                                              outputScreenIndexList_default);
    std::string delimiters(",");
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
      std::string token = str.substr(lastPos,pos-lastPos);
      outputScreenIndices.push_back(int(std::stoi(token)));
      if(pos==std::string::npos)
        break;

      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
    }

    Scalar outputScreenIndexInterval =
      pList_->get<int>(outputScreenIndexInterval_name, outputScreenIndexInterval_default);
    Scalar outputScreen_i = timeStepControl_->iStepMin;
    while (outputScreen_i <= timeStepControl_->iStepMax) {
      outputScreenIndices.push_back(outputScreen_i);
      outputScreen_i += outputScreenIndexInterval;
    }

    // order output indices
    std::sort(outputScreenIndices.begin(),outputScreenIndices.end());
  }

  return;
}


template<class Scalar>
RCP<const ParameterList> IntegratorBasic<Scalar>::getValidParameters() const
{
  static RCP<ParameterList> validPL;

  if (is_null(validPL)) {

    RCP<ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    std::ostringstream tmp;
    tmp << "'" << integratorType_name << "' must be '" << integratorType_default << "'.";
    pl->set(integratorType_name, integratorType_default, tmp.str());

    tmp.clear();
    tmp << "Initial time.  Required to be in range ["
        << timeStepControl_->timeMin << ", "<< timeStepControl_->timeMax << "].";
    pl->set(initTime_name, initTime_default, tmp.str());

    tmp << "Initial time index.  Required to be range ["
        << timeStepControl_->iStepMin << ", "<< timeStepControl_->iStepMax <<"].";
    pl->set(initTimeIndex_name, initTimeIndex_default, tmp.str());

    tmp.clear();
    tmp << "Initial time step.  Required to be positive and in range ["
        << timeStepControl_->dtMin << ", "<< timeStepControl_->dtMax << "].";
    pl->set(initTimeStep_name, initTimeStep_default, tmp.str());

    tmp.clear();
    tmp << "Initial order.  Required to be range ["
        << timeStepControl_->orderMin << ", "<< timeStepControl_->orderMax <<"].";
    pl->set(initOrder_name, initOrder_default, tmp.str());

    tmp.clear();
    tmp << "Final time.  Required to be in range ["
        << timeStepControl_->timeMin << ", "<< timeStepControl_->timeMax << "].";
    pl->set(finalTime_name, finalTime_default, tmp.str());

    tmp.clear();
    tmp << "Final time index.  Required to be range ["
        << timeStepControl_->iStepMin << ", "<< timeStepControl_->iStepMax <<"].";
    pl->set(finalTimeIndex_name, finalTimeIndex_default, tmp.str());

    tmp.clear();
    tmp << "Screen Output Index List.  Required to be range ["
        << timeStepControl_->iStepMin << ", "<< timeStepControl_->iStepMax <<"].";
    pl->set(outputScreenIndexList_name, outputScreenIndexList_default, tmp.str());
    pl->set(outputScreenIndexInterval_name, outputScreenIndexInterval_default,
      "Screen Output Index Interval (e.g., every 100 time steps");

    pl->set( StepperName_name, "",
      "'Stepper Name' selects the Stepper block to construct (Required input).\n");

    // Solution History
    ParameterList& solutionHistoryPL =
      pl->sublist(SolutionHistory_name,false,"solutionHistory_docs")
         .disableRecursiveValidation();
    solutionHistoryPL.setParameters(*(solutionHistory_->getValidParameters()));

    // Time Step Control
    ParameterList& timeStepControlPL =
      pl->sublist(TimeStepControl_name,false,"solutionHistory_docs")
         .disableRecursiveValidation();
    timeStepControlPL.setParameters(*(timeStepControl_->getValidParameters()));


    validPL = pl;
  }
  return validPL;
}


template <class Scalar>
RCP<ParameterList>
IntegratorBasic<Scalar>::getNonconstParameterList()
{
  return(pList_);
}


template <class Scalar>
RCP<ParameterList> IntegratorBasic<Scalar>::unsetParameterList()
{
  RCP<ParameterList> temp_param_list = pList_;
  pList_ = Teuchos::null;
  return(temp_param_list);
}


} // namespace Tempus
#endif // TEMPUS_INTEGRATORBASIC_IMPL_HPP
