#ifndef TEMPUS_INTEGRATORBASIC_IMPL_HPP
#define TEMPUS_INTEGRATORBASIC_IMPL_HPP

// Teuchos
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
// Tempus
#include "Tempus_StepperFactory.hpp"
#include "Tempus_TimeStepControl.hpp"

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

} // namespace


namespace Tempus {

template<class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic(
  RCP<ParameterList>                              pList_,
  const RCP<Thyra::ModelEvaluator<Scalar> >&      model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver)
{

  // Create classes from nested blocks prior to setParameters call.
  RCP<ParameterList> tmpiPL = Teuchos::sublist(pList_,"Integrator",true);

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
  RCP<SolutionState<Scalar> > currentState =
    rcp(new SolutionState<Scalar>(md, x, xdot, xdotdot,
                                  stepper->getStepperState()));
  solutionHistory->addState(currentState);

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
  out << "workingState   = " << workingState   ->description()<<std::endl;
  out << "timeStepControl= " << timeStepControl->description()<<std::endl;
  out << "stepper        = " << stepper        ->description()<<std::endl;

  if (Teuchos::as<int>(verbLevel) >=
              Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    out << "solutionHistory= " << std::endl;
    solutionHistory->describe(out,verbLevel);
    out << "workingState   = " << std::endl;
    workingState   ->describe(out,verbLevel);
    out << "timeStepControl= " << std::endl;
    timeStepControl->describe(out,verbLevel);
    out << "stepper        = " << std::endl;
    stepper        ->describe(out,verbLevel);
  }
}


template <class Scalar>
bool IntegratorBasic<Scalar>::advanceTime(const Scalar time_final)
{
  timeStepControl->timeMax = time_final;
  integratorObserver->observeStartIntegrator();

  bool integratorStatus = true;
  bool stepperStatus = true;

  while (integratorStatus == true and
         timeStepControl->timeInRange(currentState->getTime()) and
         timeStepControl->indexInRange(currentState->getIndex())){

    integratorObserver->observeStartTimeStep();

    workingState = solutionHistory->initWorkingState(stepperStatus);

    timeStepControl->getNextTimeStep(solutionHistory, stepperStatus,
                                     integratorStatus);

    integratorObserver->observeNextTimeStep(stepperStatus, integratorStatus);

    if (integratorStatus != true) break;

    integratorObserver->observeBeforeTakeStep();

    stepperStatus = stepper->takeStep(solutionHistory);

    integratorObserver->observeAfterTakeStep();

    acceptTimeStep(stepperStatus, integratorStatus);
    integratorObserver->observeAcceptedTimeStep(stepperStatus,integratorStatus);
  }

  integratorObserver->observeEndIntegrator(stepperStatus, integratorStatus);

  return integratorStatus;
}


template <class Scalar>
void IntegratorBasic<Scalar>::acceptTimeStep(
  bool & stepperStatus, bool & integratorStatus)
{
  RCP<SolutionStateMetaData<Scalar> > metaData = workingState->metaData;
  //const Scalar time = metaData->time;
  //const int iStep = metaData->iStep;
  //const Scalar errorAbs = metaData->errorAbs;
  //const Scalar errorRel = metaData->errorRel;
  //const int order = metaData->order;
  Scalar & dt = metaData->dt;
  int & nFailures = metaData->nFailures;
  int & nConsecutiveFailures = metaData->nConsecutiveFailures;

  // Stepper failure
  if (stepperStatus != true) {
    nFailures++;
    nConsecutiveFailures++;
  }

  // Constant time step failure
  if ((timeStepControl->stepType == CONSTANT_STEP_SIZE) and
      (dt != timeStepControl->dtConstant)) {
    nFailures++;
    nConsecutiveFailures++;
  }

  // Too many failures
  if (nFailures >= timeStepControl->nFailuresMax) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"continueIntegration");
    *out << "Failure - Stepper has failed more than the maximum allowed.\n"
         << "  (nFailures = "<<nFailures << ") >= (nFailuresMax = "
         <<timeStepControl->nFailuresMax<<")" << std::endl;
    integratorStatus = false;
    return;
  }
  if (nConsecutiveFailures >= timeStepControl->nConsecutiveFailuresMax) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"continueIntegration");
    *out << "Failure - Stepper has failed more than the maximum "
         << "consecutive allowed.\n"
         << "  (nConsecutiveFailures = "<<nConsecutiveFailures
         << ") >= (nConsecutiveFailuresMax = "
         <<timeStepControl->nConsecutiveFailuresMax
         << ")" << std::endl;
    integratorStatus = false;
    return;
  }

  // Made it here! Accept this time step

  nFailures = std::max(nFailures-1,0);
  nConsecutiveFailures = 0;

  solutionHistory->promoteWorkingState();

  return;
}


template <class Scalar>
void IntegratorBasic<Scalar>::setParameterList(
  const RCP<ParameterList> & pList_)
{
  if (pList_ == Teuchos::null)
    pList->validateParametersAndSetDefaults(*this->getValidParameters());
  else {
    pList = Teuchos::sublist(pList_,"Integrator",true);
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
