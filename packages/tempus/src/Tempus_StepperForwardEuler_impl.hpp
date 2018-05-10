// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperForwardEuler_impl_hpp
#define Tempus_StepperForwardEuler_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Tempus {

// StepperForwardEuler definitions:
template<class Scalar>
StepperForwardEuler<Scalar>::StepperForwardEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  // Set all the input parameters and call initialize
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->validExplicitODE(appModel);
  appModel_ = appModel;

  inArgs_  = appModel_->getNominalValues();
  outArgs_ = appModel_->createOutArgs();
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setSolver(std::string solverName)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperForwardEuler::setSolver()");
  *out << "Warning -- No solver to set for StepperForwardEuler "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperForwardEuler::setSolver()");
  *out << "Warning -- No solver to set for StepperForwardEuler "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperForwardEuler::setSolver()");
  *out << "Warning -- No solver to set for StepperForwardEuler "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperForwardEuler<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (stepperObserver_ == Teuchos::null) {
      stepperFEObserver_ =
        Teuchos::rcp(new StepperForwardEulerObserver<Scalar>());
      stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >(stepperFEObserver_);
    }
  } else {
    stepperObserver_ = obs;
    stepperFEObserver_ =
      Teuchos::rcp_dynamic_cast<StepperForwardEulerObserver<Scalar> >
        (stepperObserver_);
  }
}

template<class Scalar>
void StepperForwardEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperForwardEuler::takeStep()");
  {
    stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    typedef Thyra::ModelEvaluatorBase MEB;
    inArgs_.set_x(currentState->getX());
    if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(currentState->getTime());

    // For model evaluators whose state function f(x, x_dot, t) describes
    // an implicit ODE, and which accept an optional x_dot input argument,
    // make sure the latter is set to null in order to request the evaluation
    // of a state function corresponding to the explicit ODE formulation
    // x_dot = f(x, t)
    if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);
    RCP<Thyra::VectorBase<Scalar> > xDot = currentState->getXDot();
    if (xDot == Teuchos::null) xDot = getXDotTemp(currentState->getX());
    outArgs_.set_f(xDot);

    if (!Teuchos::is_null(stepperFEObserver_))
      stepperFEObserver_->observeBeforeExplicit(solutionHistory, *this);

    appModel_->evalModel(inArgs_,outArgs_);

    // Forward Euler update, x = x + dt*xdot
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    const Scalar dt = workingState->getTimeStep();
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
      *(currentState->getX()),dt,*(xDot));

    if (workingState->getXDot() != Teuchos::null)
      assign((workingState->getXDot()).ptr(),
        Teuchos::ScalarTraits<Scalar>::zero());
    workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    workingState->setOrder(this->getOrder());
    stepperObserver_->observeEndTakeStep(solutionHistory, *this);
  }
  return;
}


template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperForwardEuler<Scalar>::
getXDotTemp(Teuchos::RCP<Thyra::VectorBase<Scalar> > x)
{
  if (xDotTemp_ == Teuchos::null) {
    xDotTemp_ = x->clone_v();
    Thyra::assign(xDotTemp_.ptr(), Scalar(0.0));
  }
  return xDotTemp_;
}


/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> > StepperForwardEuler<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperForwardEuler<Scalar>::description() const
{
  std::string name = "Forward Euler";
  return(name);
}


template<class Scalar>
void StepperForwardEuler<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "appModel_ = " << appModel_->description() << std::endl;
}


template <class Scalar>
void StepperForwardEuler<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (stepperPL_ == Teuchos::null) stepperPL_ = this->getDefaultParameters();
  } else {
    stepperPL_ = pList;
  }
  stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string stepperType = stepperPL_->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Forward Euler",
    std::logic_error,
       "Error - Stepper Type is not 'Forward Euler'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperForwardEuler<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set("Stepper Type", "Forward Euler",
          "'Stepper Type' must be 'Forward Euler'.");

  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperForwardEuler<Scalar>::getDefaultParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  *pl = *(this->getValidParameters());
  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperForwardEuler<Scalar>::getNonconstParameterList()
{
  return(stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperForwardEuler<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperForwardEuler_impl_hpp
