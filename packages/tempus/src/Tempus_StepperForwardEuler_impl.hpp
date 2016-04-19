#ifndef TEMPUS_STEPPERFORWARDEULER_IMPL_HPP
#define TEMPUS_STEPPERFORWARDEULER_IMPL_HPP

#include "Tempus_StepperForwardEuler.hpp"


namespace Tempus {

// StepperForwardEuler definitions:
template<class Scalar>
StepperForwardEuler<Scalar>::StepperForwardEuler(
  RCP<ParameterList>                               pList_,
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model_ )
 : model(model_)
{
  inArgs  = model->createInArgs();
  outArgs = model->createOutArgs();
  inArgs  = model->getNominalValues();

  stepperState = rcp(new StepperState<Scalar>(description()));
}

template<class Scalar>
bool takeStep(const Ptr<SolutionHistory<Scalar> >& solutionHistory);
{
  const Ptr<SolutionState<Scalar> > workingState =
    solutionHistory->getWorkingState();

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(workingState), std::logic_error,
    "Error - SolutionState, workingstate, is invalid!\n");

  typedef Thyra::ModelEvaluatorBase MEB;
  const Scalar time = workingState->getTime();
  const Scalar dt   = workingState->getTimeStep();

  inArgs.set_x(workingState->getX());
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(time);

  // For model evaluators whose state function f(x, x_dot, t) describes
  // an implicit ODE, and which accept an optional x_dot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // x_dot = f(x, t)
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(Teuchos::null);
  outArgs.set_f(workingState->getXDot());

  model->evalModel(inArgs,outArgs);

  // Forward Euler update, x = x + dt*xdot
  Thyra::Vp_StV(workingState->getX().ptr(),dt,*(workingState->getXDot));

  return (true);
}


template<class Scalar>
void StepperForwardEuler::setStepperState(
  const RCP<Tempus::StepperState<Scalar> >& stepperState_)
{
  TEUCHOS_TEST_FOR_EXCEPTION(stepperState_->stepperName != description(),
    std::logic_error,
    "Error - StepperState does not match Stepper!\n"
    "        stepperState = " << stepperState_->stepperName << "\n"
    "        stepper      = " << description() << std::endl);

  // ForwardEuler does not have any state information, but if it did we
  // would set the internal stepperState to the input stepperState_ here.
}


template<class Scalar>
RCP<Tempus::StepperState<Scalar> > getStepperState()
{
  return stepperState;
}


template<class Scalar>
std::string StepperForwardEuler<Scalar>::description() const
{
  std::string name = "Tempus::StepperForwardEuler";
  return(name);
}


template<class Scalar>
void StepperForwardEuler<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "model = " << model->description() << std::endl;
}


template <class Scalar>
void StepperForwardEuler<Scalar>::setParameterList(
  RCP<ParameterList> const& pList_)
{
  TEUCHOS_TEST_FOR_EXCEPT(!is_null(pList_));
  pList_->validateParameters(*this->getValidParameters());
  pList = pList_;

  Teuchos::readVerboseObjectSublist(&*pList,this);
}


template<class Scalar>
RCP<const ParameterList> StepperForwardEuler<Scalar>::getValidParameters() const
{
  static RCP<ParameterList> validPL;
  return validPL;
}


template <class Scalar>
RCP<ParameterList>
StepperForwardEuler<Scalar>::getNonconstParameterList()
{
  return(pList);
}


template <class Scalar>
RCP<ParameterList> StepperForwardEuler<Scalar>::unsetParameterList()
{
  RCP<ParameterList> temp_plist = pList;
  pList = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif TEMPUS_STEPPERFORWARDEULER_IMPL_HPP
