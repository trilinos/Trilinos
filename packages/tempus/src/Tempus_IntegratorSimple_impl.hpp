#ifndef TEMPUS_INTEGRATORSIMPLE_IMPL_HPP
#define TEMPUS_INTEGRATORSIMPLE_IMPL_HPP

#include "Tempus_IntegratorSimple.hpp"


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
IntegratorSimple<Scalar>::IntegratorSimple(
  RCP<ParameterList>              parameterList,
  RCP<Thyra::VectorBase<Scalar> > x,
  RCP<Thyra::VectorBase<Scalar> > xdot=Teuchos::null,
  RCP<Thyra::VectorBase<Scalar> > xdotdot=Teuchos::null )
{
  error    = 0.0;
  accuracy = 0.0;

  // Create TimeStepControl before ParameterList so we can check
  // ParameterList ranges.
  RCP<ParameterList> tsc_pl = Teuchos::sublist(paramList, TimeStepControl_name);
  timeStepControl = rcp(new TimeStepControl<Scalar>(tsc_pl));

  if ( paramList_ == Teuchos::null )
    paramList = getValidParameters();
  else
    paramList = paramList_;

  this->setParameterList(paramList);

  // Create initial condition solution state
  currentState = rcp(new SolutionState<Scalar>(x, xdot, xdotdot, time,
                                               timeStepControl->dtMin,
                                               timeStepControl->dtMax,
                                               iStep, 0, error, false,
                                               accuracy));

  // Create solution history, an array of solution states.
  RCP<ParameterList> sh_pl = Teuchos::sublist(paramList, SolutionHistory_name);
  solutionHistory = rcp(new SolutionHistory<Scalar>(sh_pl));
  addState(currentState);

  // Create Stepper
  std::string s = paramList->get<std::string>(Stepper_name, Stepper_default);
  RCP<ParameterList> s_pl = Teuchos::sublist(paramList, s);
  stepper = rcp(new Stepper<Scalar>(s_pl));

  if ( Teuchos::as<int>(this->getVerbLevel()) >=
       Teuchos::as<int>(Teuchos::VERB_HIGH) ) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"IntegratorSimple::IntegratorSimple");
    *out << this->description() << std::endl;
  }
}


template<class Scalar>
std::string IntegratorSimple<Scalar>::description() const
{
  std::string name = "tempus::IntegratorSimple";
  return(name);
}


template<class Scalar>
void IntegratorSimple<Scalar>::describe(
  Teuchos::FancyOStream          &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  out << description() << "::describe" << std::endl;
  out << "time   = " << time  << std::endl;
  out << "dt     = " << dt    << std::endl;
  out << "iStep  = " << iStep << std::endl;
  out << "error  = " << error  << std::endl;
  out << "order  = " << order  << std::endl;
  if ((Teuchos::as<int>(verbLevel)==Teuchos::as<int>(Teuchos::VERB_DEFAULT)) ||
      (Teuchos::as<int>(verbLevel)>=Teuchos::as<int>(Teuchos::VERB_LOW)    )  ){
    out << "solutionHistory    = " << solutionHistory->description()<<std::endl;
    out << "currentState       = " << currentState   ->description()<<std::endl;
    out << "timeStepControl    = " << timeStepControl->description()<<std::endl;
    out << "integratorObserver = " << integratorObserver->description()
        << std::endl;
    out << "stepper            = " << stepper        ->description()<<std::endl;
  } else if ( Teuchos::as<int>(verbLevel) >=
              Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    out << "solutionHistory    = " << solutionHistory   ->describe()<<std::endl;
    out << "currentState       = " << currentState      ->describe()<<std::endl;
    out << "timeStepControl    = " << timeStepControl   ->describe()<<std::endl;
    out << "integratorObserver = " << integratorObserver->describe()<<std::endl;
    out << "stepper            = " << stepper           ->describe()<<std::endl;
  }
}


template <class Scalar>
void IntegratorSimple<Scalar>::advanceTime(const Scalar time_final)
{
  integratorObserver->observeStartTime();

  while ( timeStepControl->timeInRange(time) and
          timeStepControl->indexInRange(iStep) ){

    integratorStatus = timeStepControl->getNextTimeStep();
    if (integratorStatus == FAILED) {
      integratorObserver->observeFailedIntegrator();
      break;
    }

    integratorObserver->observeStartTimeStep();

    stepperStatus = stepper->takeStep();

    if (stepperStatus == FAILED) {
      integratorObserver->observeFailedTimeStep();
      continue;
    }

    acceptTimeStep();
    integratorObserver->observeAcceptedTimeStep();

    printTimeStep();
  }

  integratorObserver->observeEndTime();
}


template <class Scalar>
void IntegratorSimple<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList_)
{
  TEUCHOS_TEST_FOR_EXCEPT( is_null(paramList_) );
  paramList_->validateParameters(*this->getValidParameters());
  paramList = paramList_;

  Teuchos::readVerboseObjectSublist(&*paramList,this);

  time = paramList->get<double>(initTime_name, initTime_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (time < timeStepControl->timeMin || time > timeStepControl->timeMax ),
    std::out_of_range,
    "Error - Time is out of range.\n"
    << "    [timeMin, timeMax] = [" << timeStepControl->timeMin << ", "
                                    << timeStepControl->timeMax << "]\n"
    << "    time = " << time << "\n");

  dt    = paramList->get<double>(initTimeStep_name, initTimeStep_default);
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

  iStep = paramList->get<int>(initTimeIndex_name, initTimeIndex_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (iStep < timeStepControl->iStepMin || iStep > timeStepControl->iStepMax),
    std::out_of_range,
    "Error - Time step is out of range.\n"
    << "    [iStepMin, iStepMax] = [" << timeStepControl->iStepMin << ", "
                                      << timeStepControl->iStepMax << "]\n"
    << "    iStep = " << iStep << "\n");

  order = paramList->get<int>(initOrder_name, initOrder_default);
  TEUCHOS_TEST_FOR_EXCEPTION(
    (order < timeStepControl->orderMin || order > timeStepControl->orderMax),
    std::out_of_range,
    "Error - Time step is out of range.\n"
    << "    [orderMin, orderMax] = [" << timeStepControl->orderMin << ", "
                                      << timeStepControl->orderMax << "]\n"
    << "    order = " << order << "\n");
}


template<class Scalar>
RCP<const Teuchos::ParameterList> IntegratorSimple<Scalar>::getValidParameters() const
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
IntegratorSimple<Scalar>::getNonconstParameterList()
{
  return(paramList);
}


template <class Scalar>
RCP<Teuchos::ParameterList> IntegratorSimple<Scalar>::unsetParameterList()
{
  RCP<Teuchos::ParameterList> temp_param_list = paramList;
  paramList = Teuchos::null;
  return(temp_param_list);
}


} // namespace tempus
#endif // TEMPUS_INTEGRATORSIMPLE_IMPL_HPP
