#ifndef TEMPUS_INTEGRATORBASIC_IMPL_HPP
#define TEMPUS_INTEGRATORBASIC_IMPL_HPP

#include "Tempus_IntegratorBasic.hpp"


namespace {

  static std::string SolutionHistory_name     = "Solution History";
  static std::string TimeStepControl_name     = "Time Step Control";

  static std::string initTime_name         = "Initial Time";
  static Scalar      initTime_default      = 0.0;
  static std::string initTimeStep_name     = "Initial Time Step";
  static Scalar      initTimeStep_default  = 0.0;
  static std::string initTimeIndex_name    = "Initial Time Index";
  static int         initTimeIndex_default = 0;
  static std::string initOrder_name        = "Initial Order";
  static int         initOrder_default     = 1;

  static std::string ForwardEuler_name = "Forward Euler";
  static std::string Stepper_name      = "Stepper";
  static std::string Stepper_default   = ForwardEuler_name;

  Teuchos::Array<std::string> Stepper_names = Teuchos::tuple<std::string>(
      ForwardEuler_name);

  const Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<tempus::StepType> >
    StepperValidator = Teuchos::rcp(
        new Teuchos::StringToIntegralParameterEntryValidator<tempus::StepType>(
          Stepper_names,
          Teuchos::tuple<tempus::Steppers>(
            tempus::FORWARD_EULER),
          Stepper_name));

} // namespace


namespace tempus {


template<class Scalar>
IntegratorBasic<Scalar>::IntegratorBasic(
  RCP<ParameterList>&                    pList_,
  const RCP<Thyra::VectorBase<Scalar> >& x,
  const RCP<Thyra::VectorBase<Scalar> >& xdot=Teuchos::null,
  const RCP<Thyra::VectorBase<Scalar> >& xdotdot=Teuchos::null )
{
  // Create TimeStepControl before ParameterList so we can check
  // ParameterList ranges.
  RCP<ParameterList> tsc_pl = Teuchos::sublist(pList, TimeStepControl_name);
  timeStepControl = rcp(new TimeStepControl<Scalar>(tsc_pl));

  if ( pList_ == Teuchos::null )
    pList = getValidParameters();
  else
    pList = pList_;

  this->setParameterList(pList);

  // Create meta data
  RCP<SolutionStateMetaData<Scalar> > md =
                                   rcp(new SolutionStateMetaData<Scalar> ());
  md->time = pList->get<double>(initTime_name, initTime_default);
  md->dt = pList->get<double>(initTimeStep_name, initTimeStep_default);
  md->iStep = pList->get<int>(initTimeIndex_name, initTimeIndex_default);
  md->order = pList->get<int>(initOrder_name, initOrder_default);

  // Create Stepper
  std::string s = pList->get<std::string>(Stepper_name, Stepper_default);
  RCP<ParameterList> s_pl = Teuchos::sublist(pList, s);
  stepper = rcp(new Stepper<Scalar>(s_pl));

  // Create working solution state
  workingState = rcp(new SolutionState<Scalar>(md, x, xdot, xdotdot,
                                               stepper->getStepperState());

  // Create solution history
  RCP<ParameterList> sh_pl = Teuchos::sublist(pList, SolutionHistory_name);
  solutionHistory = rcp(new SolutionHistory<Scalar>(sh_pl));
  solutionHistory->addState(workingState);

  if ( Teuchos::as<int>(this->getVerbLevel()) >=
       Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"IntegratorBasic::IntegratorBasic");
    *out << this->description() << std::endl;
  }
}


template<class Scalar>
std::string IntegratorBasic<Scalar>::description() const
{
  std::string name = "tempus::IntegratorBasic";
  return(name);
}


template<class Scalar>
void IntegratorBasic<Scalar>::describe(
  Teuchos::FancyOStream          &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  out << description() << "::describe" << std::endl;
  out << "solutionHistory    = " << solutionHistory->description()<<std::endl;
  out << "workingState       = " << workingState   ->description()<<std::endl;
  out << "timeStepControl    = " << timeStepControl->description()<<std::endl;
  out << "integratorObserver = " <<integratorObserver->description()<<std::endl;
  out << "stepper            = " << stepper        ->description()<<std::endl;

  if ( Teuchos::as<int>(verbLevel) >=
              Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    out << "solutionHistory    = " << solutionHistory   ->describe()<<std::endl;
    out << "workingState       = " << workingState      ->describe()<<std::endl;
    out << "timeStepControl    = " << timeStepControl   ->describe()<<std::endl;
    out << "integratorObserver = " << integratorObserver->describe()<<std::endl;
    out << "stepper            = " << stepper           ->describe()<<std::endl;
  }
}


template <class Scalar>
void IntegratorBasic<Scalar>::advanceTime(const Scalar time_final)
{
  timeStepControl->timeMax = time_final;
  integratorObserver->observeStartIntegrator();

  bool integratorStatus = true;
  bool stepperStatus = true;

  while ( timeStepControl->timeInRange(workingState->getTime()) and
          timeStepControl->indexInRange(workingState->getIndex()) ){

    integratorObserver->observeStartTimeStep();

    workingState = solutionHistory->setWorkingState();

    timeStepControl->getNextTimeStep(workingState->metaData,
                                     stepperStatus, integratorStatus);

    if (integratorStatus != true) {
      integratorObserver->observeFailedIntegrator();
      break;
    }

    integratorObserver->observeBeforeTakeStep();

    stepperStatus = stepper->takeStep(solutionHistory);

    if (stepperStatus != true) {
      integratorObserver->observeFailedTimeStep();
      continue;
    }

    acceptTimeStep();
    integratorObserver->observeAcceptedTimeStep();
  }

  integratorObserver->observeEndIntegrator();
}


template <class Scalar>
void IntegratorBasic<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& pList_)
{
  TEUCHOS_TEST_FOR_EXCEPT( is_null(pList_) );
  pList_->validateParameters(*this->getValidParameters());
  pList = pList_;

  Teuchos::readVerboseObjectSublist(&*pList,this);

  Scalar time = pList->get<double>(initTime_name, initTime_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (time < timeStepControl->timeMin || time > timeStepControl->timeMax ),
    std::out_of_range,
    "Error - Time is out of range.\n"
    << "    [timeMin, timeMax] = [" << timeStepControl->timeMin << ", "
                                    << timeStepControl->timeMax << "]\n"
    << "    time = " << time << "\n");

  Scalar dt = pList->get<double>(initTimeStep_name, initTimeStep_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dt < ST::zero() ), std::logic_error,
    "Error - Negative time step.  dt = "<<dt<<")\n");
  TEUCHOS_TEST_FOR_EXCEPTION(
    (dt < timeStepControl->dtMin || dt > timeStepControl->dtMax ),
    std::out_of_range,
    "Error - Time step is out of range.\n"
    << "    [dtMin, dtMax] = [" << timeStepControl->dtMin << ", "
                                << timeStepControl->dtMax << "]\n"
    << "    dt = " << dt << "\n");

  int iStep = pList->get<int>(initTimeIndex_name, initTimeIndex_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (iStep < timeStepControl->iStepMin || iStep > timeStepControl->iStepMax),
    std::out_of_range,
    "Error - Time step is out of range.\n"
    << "    [iStepMin, iStepMax] = [" << timeStepControl->iStepMin << ", "
                                      << timeStepControl->iStepMax << "]\n"
    << "    iStep = " << iStep << "\n");

  int order = pList->get<int>(initOrder_name, initOrder_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (order < timeStepControl->orderMin || order > timeStepControl->orderMax),
    std::out_of_range,
    "Error - Time step is out of range.\n"
    << "    [orderMin, orderMax] = [" << timeStepControl->orderMin << ", "
                                      << timeStepControl->orderMax << "]\n"
    << "    order = " << order << "\n");
}


template<class Scalar>
RCP<const Teuchos::ParameterList> IntegratorBasic<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    pl->set(initTime_name, initTime_default,
      "Initial time.  Required to be in range ["
      << timeStepControl->timeMin << ", "<< timeStepControl->timeMax << "].");
    pl->set(initTimeStep_name, initTimeStep_default,
      "Initial time step.  Required to be positive and in range ["
      << timeStepControl->dtMin << ", "<< timeStepControl->dtMax << "].");
    pl->set(initTimeIndex_name, initTimeIndex_default,
      "Initial time index.  Required to be range ["
      << timeStepControl->iStepMin << ", "<< timeStepControl->iStepMax << "].");
    pl->set(initOrder_name, initOrder_default,
      "Initial order.  Required to be range ["
      << timeStepControl->orderMin << ", "<< timeStepControl->orderMax << "].");

    pl->set(Stepper_name, Stepper_default,
      Stepper_name <<" sets the Stepper.\n"
      "'Forward Euler' - performs classic first-order Forward Euler\n",
      StepperValidator);

    validPL = pl;

  }
  return validPL;
}


template <class Scalar>
RCP<Teuchos::ParameterList>
IntegratorBasic<Scalar>::getNonconstParameterList()
{
  return(pList);
}


template <class Scalar>
RCP<Teuchos::ParameterList> IntegratorBasic<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = pList;
  pList = Teuchos::null;
  return(temp_param_list);
}


} // namespace tempus
#endif // TEMPUS_INTEGRATORBASIC_IMPL_HPP
