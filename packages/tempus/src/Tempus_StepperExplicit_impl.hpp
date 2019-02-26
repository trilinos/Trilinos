// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicit_impl_hpp
#define Tempus_StepperExplicit_impl_hpp


namespace Tempus {


template<class Scalar>
void StepperExplicit<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->validExplicitODE(appModel);
  this->appModel_ = appModel;

  this->inArgs_  = this->appModel_->getNominalValues();
  this->outArgs_ = this->appModel_->createOutArgs();
}


template<class Scalar>
void StepperExplicit<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}

template<class Scalar>
void StepperExplicit<Scalar>::setSolver(std::string solverName)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperExplicit::setSolver()");
  *out << "Warning -- No solver to set for StepperExplicit "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperExplicit<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperExplicit::setSolver()");
  *out << "Warning -- No solver to set for StepperExplicit "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
void StepperExplicit<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperExplicit::setSolver()");
  *out << "Warning -- No solver to set for StepperExplicit "
       << "(i.e., explicit method).\n" << std::endl;
  return;
}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperExplicit<Scalar>::
getStepperX(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getX() != Teuchos::null) stepperX_ = state->getX();
  // Else use temporary storage stepperXp_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperX_ == Teuchos::null, std::logic_error,
    "Error - stepperX_ has not been set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperX_;
}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperExplicit<Scalar>::
getStepperXDot(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getXDot() != Teuchos::null) stepperXDot_ = state->getXDot();
  // Else use temporary storage stepperXDot_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperXDot_ == Teuchos::null, std::logic_error,
    "Error - stepperXDot_ has not set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperXDot_;
}

template<class Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
StepperExplicit<Scalar>::
getStepperXDotDot(Teuchos::RCP<SolutionState<Scalar> > state)
{
  if (state->getXDotDot() != Teuchos::null) stepperXDotDot_=state->getXDotDot();
  // Else use temporary storage stepperXDotDot_ which should have been set in
  // setInitialConditions().

  TEUCHOS_TEST_FOR_EXCEPTION( stepperXDotDot_==Teuchos::null, std::logic_error,
    "Error - stepperXDotDot_ has not set in setInitialConditions() or\n"
    "        can not be set from the state!\n");

  return stepperXDotDot_;
}

template<class Scalar>
void
StepperExplicit<Scalar>::
evaluateExplicitODE(Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
                    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
                    const Scalar time)
{
  typedef Thyra::ModelEvaluatorBase MEB;

  inArgs_.set_x(x);
  if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time);

  // For model evaluators whose state function f(x, xDot, t) describes
  // an implicit ODE, and which accept an optional xDot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // xDot = f(x, t).
  if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);

  outArgs_.set_f(xDot);

  appModel_->evalModel(inArgs_, outArgs_);
}

template<class Scalar>
void
StepperExplicit<Scalar>::
evaluateExplicitODE(Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot,
                    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
                    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDot,
                    const Scalar time)
{
  typedef Thyra::ModelEvaluatorBase MEB;

  inArgs_.set_x(x);
  if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(xDot);
  if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time);

  // For model evaluators whose state function f(x, xDot, xDotDot, t) describes
  // an implicit ODE, and which accept an optional xDotDot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // xDotDot = f(x, xDot, t).
  if (inArgs_.supports(MEB::IN_ARG_x_dot_dot))
    inArgs_.set_x_dot_dot(Teuchos::null);

  outArgs_.set_f(xDotDot);

  appModel_->evalModel(inArgs_, outArgs_);
}


} // namespace Tempus
#endif // Tempus_StepperExplicit_impl_hpp
