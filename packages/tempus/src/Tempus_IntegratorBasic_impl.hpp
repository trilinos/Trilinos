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

  static std::string SolutionHistory_name     = "Solution History";
  static std::string TimeStepControl_name     = "Time Step Control";

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

} // namespace


namespace Tempus {

template<class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic(
  RCP<ParameterList>                              pList_,
  const RCP<Thyra::ModelEvaluator<Scalar> >&      model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver)
     : integratorStatus (WORKING)
{

  // Create classes from nested blocks prior to setParameters call.
  RCP<ParameterList> tmpiPL = Teuchos::sublist(pList_,"Integrator Basic",true);

  //    Create solution history
  RCP<ParameterList> shPL = Teuchos::sublist(tmpiPL, SolutionHistory_name,true);
  solutionHistory = rcp(new SolutionHistory<Scalar>(shPL));

  //    Create TimeStepControl
  RCP<ParameterList> tscPL = Teuchos::sublist(tmpiPL,TimeStepControl_name,true);
  Scalar dtTmp = tmpiPL->get<double>(initTimeStep_name, initTimeStep_default);
  timeStepControl = rcp(new TimeStepControl<Scalar>(tscPL, dtTmp));

  this->setParameterList(pList_);

  // Create Stepper
  RCP<StepperFactory<Scalar> > sf = rcp(new StepperFactory<Scalar>());
  std::string s = pList->get<std::string>(getStepperName(),getStepperDefault());
  RCP<ParameterList> s_pl = Teuchos::sublist(pList, s);
  stepper = sf->createStepper(s, s_pl, model);

  // Create meta data
  RCP<SolutionStateMetaData<Scalar> > md =
                                   rcp(new SolutionStateMetaData<Scalar> ());
  md->time  = pList->get<double>(initTime_name,      initTime_default);
  md->iStep = pList->get<int>   (initTimeIndex_name, initTimeIndex_default);
  md->dt    = pList->get<double>(initTimeStep_name,  initTimeStep_default);
  md->order = pList->get<int>   (initOrder_name,     initOrder_default);

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
    md, x, xdot, xdotdot, stepper->getStepperState()));
  solutionHistory->addState(newState);

  // Create default IntegratorObserver
  integratorObserver = rcp(new IntegratorObserver<Scalar>(solutionHistory,
                                                          timeStepControl));

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
  out << "solutionHistory= " << solutionHistory->description()<<std::endl;
  out << "timeStepControl= " << timeStepControl->description()<<std::endl;
  out << "stepper        = " << stepper        ->description()<<std::endl;

  if (Teuchos::as<int>(verbLevel) >=
              Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    out << "solutionHistory= " << std::endl;
    solutionHistory->describe(out,verbLevel);
    out << "timeStepControl= " << std::endl;
    timeStepControl->describe(out,verbLevel);
    out << "stepper        = " << std::endl;
    stepper        ->describe(out,verbLevel);
  }
}


template <class Scalar>
bool IntegratorBasic<Scalar>::advanceTime(const Scalar timeFinal)
{
  if (timeStepControl->timeInRange(timeFinal))
    timeStepControl->timeMax = timeFinal;
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
       << "  Stepper = " << stepper->description() << "\n"
       << "  Simulation Time Range  [" << timeStepControl->timeMin
       << ", " << timeStepControl->timeMax << "]\n"
       << "  Simulation Index Range [" << timeStepControl->iStepMin
       << ", " << timeStepControl->iStepMax << "]\n"
       << "============================================================================\n"
       << "  Step       Time         dt  Abs Error  Rel Error  Order  nFail  dCompTime"
       << std::endl;
  integratorStatus = WORKING;
}


template <class Scalar>
bool IntegratorBasic<Scalar>::advanceTime()
{
  startIntegrator();
  integratorObserver->observeStartIntegrator();

  while (integratorStatus == WORKING and
         timeStepControl->timeInRange (solutionHistory->getCurrentTime()) and
         timeStepControl->indexInRange(solutionHistory->getCurrentIndex())){

    stepperTimer->reset();
    stepperTimer->start();
    integratorObserver->observeStartTimeStep();

    solutionHistory->initWorkingState();

    timeStepControl->getNextTimeStep(solutionHistory, integratorStatus);

    integratorObserver->observeNextTimeStep(integratorStatus);

    if (integratorStatus == FAILED) break;

    integratorObserver->observeBeforeTakeStep();

    stepper->takeStep(solutionHistory);

    integratorObserver->observeAfterTakeStep();

    stepperTimer->stop();
    acceptTimeStep();
    integratorObserver->observeAcceptedTimeStep(integratorStatus);
  }

  endIntegrator();
  integratorObserver->observeEndIntegrator(integratorStatus);

  return (integratorStatus == Status::PASSED);
}


template <class Scalar>
void IntegratorBasic<Scalar>::acceptTimeStep()
{
  RCP<SolutionStateMetaData<Scalar> > wsmd =
    solutionHistory->getWorkingState()->metaData;

       // Stepper failure
  if ( solutionHistory->getWorkingState()->getSolutionStatus() == FAILED or
       solutionHistory->getWorkingState()->getStepperStatus() == FAILED or
       // Constant time step failure
       ((timeStepControl->stepType == CONSTANT_STEP_SIZE) and
       (wsmd->dt != timeStepControl->dtConstant))
     )
  {
    wsmd->nFailures++;
    wsmd->nConsecutiveFailures++;
    wsmd->solutionStatus = FAILED;
  }

  // Too many failures
  if (wsmd->nFailures >= timeStepControl->nFailuresMax) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"continueIntegration");
    *out << "Failure - Stepper has failed more than the maximum allowed.\n"
         << "  (nFailures = "<<wsmd->nFailures << ") >= (nFailuresMax = "
         <<timeStepControl->nFailuresMax<<")" << std::endl;
    integratorStatus = FAILED;
    return;
  }
  if (wsmd->nConsecutiveFailures >= timeStepControl->nConsecutiveFailuresMax) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"continueIntegration");
    *out << "Failure - Stepper has failed more than the maximum "
         << "consecutive allowed.\n"
         << "  (nConsecutiveFailures = "<<wsmd->nConsecutiveFailures
         << ") >= (nConsecutiveFailuresMax = "
         <<timeStepControl->nConsecutiveFailuresMax
         << ")" << std::endl;
    integratorStatus = FAILED;
    return;
  }

  // =======================================================================
  // Made it here! Accept this time step

  solutionHistory->promoteWorkingState();

  RCP<SolutionStateMetaData<Scalar> > csmd =
    solutionHistory->getCurrentState()->metaData;

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
  if (solutionHistory->getCurrentState()->getSolutionStatus() ==
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
  const RCP<ParameterList> & pList_)
{
  if (pList_ == Teuchos::null)
    pList->validateParametersAndSetDefaults(*this->getValidParameters());
  else {
    pList = Teuchos::sublist(pList_,"Integrator Basic",true);
    pList->validateParameters(*this->getValidParameters());
  }

  Teuchos::readVerboseObjectSublist(&*pList,this);

  Scalar initTime = pList->get<double>(initTime_name, initTime_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (initTime<timeStepControl->timeMin || initTime>timeStepControl->timeMax),
    std::out_of_range,
    "Error - Initial time is out of range.\n"
    << "    [timeMin, timeMax] = [" << timeStepControl->timeMin << ", "
                                    << timeStepControl->timeMax << "]\n"
    << "    initTime = " << initTime << "\n");

  int iStep = pList->get<int>(initTimeIndex_name, initTimeIndex_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (iStep < timeStepControl->iStepMin || iStep > timeStepControl->iStepMax),
    std::out_of_range,
    "Error - Initial time index is out of range.\n"
    << "    [iStepMin, iStepMax] = [" << timeStepControl->iStepMin << ", "
                                      << timeStepControl->iStepMax << "]\n"
    << "    iStep = " << iStep << "\n");

  Scalar dt = pList->get<double>(initTimeStep_name, initTimeStep_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dt < Teuchos::ScalarTraits<Scalar>::zero() ), std::logic_error,
    "Error - Negative time step.  dt = "<<dt<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dt < timeStepControl->dtMin || dt > timeStepControl->dtMax ),
    std::out_of_range,
    "Error - Initial time step is out of range.\n"
    << "    [dtMin, dtMax] = [" << timeStepControl->dtMin << ", "
                                << timeStepControl->dtMax << "]\n"
    << "    dt = " << dt << "\n");

  int order = pList->get<int>(initOrder_name, initOrder_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (order < timeStepControl->orderMin || order > timeStepControl->orderMax),
    std::out_of_range,
    "Error - Initial order is out of range.\n"
    << "    [orderMin, orderMax] = [" << timeStepControl->orderMin << ", "
                                      << timeStepControl->orderMax << "]\n"
    << "    order = " << order << "\n");

  Scalar finalTime = pList->get<double>(finalTime_name, finalTime_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (finalTime<timeStepControl->timeMin || finalTime>timeStepControl->timeMax),
    std::out_of_range,
    "Error - Final time is out of range.\n"
    << "    [timeMin, timeMax] = [" << timeStepControl->timeMin << ", "
                                    << timeStepControl->timeMax << "]\n"
    << "    finalTime = " << finalTime << "\n");

  int fiStep = pList->get<int>(finalTimeIndex_name, finalTimeIndex_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (fiStep < timeStepControl->iStepMin || fiStep > timeStepControl->iStepMax),
    std::out_of_range,
    "Error - Final time index is out of range.\n"
    << "    [iStepMin, iStepMax] = [" << timeStepControl->iStepMin << ", "
                                      << timeStepControl->iStepMax << "]\n"
    << "    iStep = " << fiStep << "\n");

  // Parse output indices
  {
    outputScreenIndices.clear();
    std::string str = pList->get<std::string>(outputScreenIndexList_name,
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
      pList->get<int>(outputScreenIndexInterval_name, outputScreenIndexInterval_default);
    Scalar outputScreen_i = timeStepControl->iStepMin;
    while (outputScreen_i <= timeStepControl->iStepMax) {
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
    tmp << "Initial time.  Required to be in range ["
        << timeStepControl->timeMin << ", "<< timeStepControl->timeMax << "].";
    pl->set(initTime_name, initTime_default, tmp.str());

    tmp.clear();
    tmp << "Initial time index.  Required to be range ["
        << timeStepControl->iStepMin << ", "<< timeStepControl->iStepMax <<"].";
    pl->set(initTimeIndex_name, initTimeIndex_default, tmp.str());

    tmp.clear();
    tmp << "Initial time step.  Required to be positive and in range ["
        << timeStepControl->dtMin << ", "<< timeStepControl->dtMax << "].";
    pl->set(initTimeStep_name, initTimeStep_default, tmp.str());

    tmp.clear();
    tmp << "Initial order.  Required to be range ["
        << timeStepControl->orderMin << ", "<< timeStepControl->orderMax <<"].";
    pl->set(initOrder_name, initOrder_default, tmp.str());

    tmp.clear();
    tmp << "Final time.  Required to be in range ["
        << timeStepControl->timeMin << ", "<< timeStepControl->timeMax << "].";
    pl->set(finalTime_name, finalTime_default, tmp.str());

    tmp.clear();
    tmp << "Final time index.  Required to be range ["
        << timeStepControl->iStepMin << ", "<< timeStepControl->iStepMax <<"].";
    pl->set(finalTimeIndex_name, finalTimeIndex_default, tmp.str());

    tmp.clear();
    tmp << "Screen Output Index List.  Required to be range ["
        << timeStepControl->iStepMin << ", "<< timeStepControl->iStepMax <<"].";
    pl->set(outputScreenIndexList_name, outputScreenIndexList_default, tmp.str());
    pl->set(outputScreenIndexInterval_name, outputScreenIndexInterval_default,
      "Screen Output Index Interval (e.g., every 100 time steps");

    pl->set( getStepperName(), getStepperDefault(),
      "'Stepper' sets the Stepper.\n"
      "'Forward Euler' - performs classic first-order Forward Euler\n",
      StepperValidator);

    // Solution History
    ParameterList& solutionHistoryPL =
      pl->sublist(SolutionHistory_name,false,"solutionHistory_docs")
         .disableRecursiveValidation();
    solutionHistoryPL.setParameters(*(solutionHistory->getValidParameters()));

    // Time Step Control
    ParameterList& timeStepControlPL =
      pl->sublist(TimeStepControl_name,false,"solutionHistory_docs")
         .disableRecursiveValidation();
    timeStepControlPL.setParameters(*(timeStepControl->getValidParameters()));


    validPL = pl;
  }
  return validPL;
}


template <class Scalar>
RCP<ParameterList>
IntegratorBasic<Scalar>::getNonconstParameterList()
{
  return(pList);
}


template <class Scalar>
RCP<ParameterList> IntegratorBasic<Scalar>::unsetParameterList()
{
  RCP<ParameterList> temp_param_list = pList;
  pList = Teuchos::null;
  return(temp_param_list);
}


} // namespace Tempus
#endif // TEMPUS_INTEGRATORBASIC_IMPL_HPP
