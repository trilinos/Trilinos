#ifndef TEMPUS_STEPPERFORWARDEULER_IMPL_HPP
#define TEMPUS_STEPPERFORWARDEULER_IMPL_HPP

// Teuchos
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Tempus {

// StepperForwardEuler definitions:
template<class Scalar>
StepperForwardEuler<Scalar>::StepperForwardEuler(
  Teuchos::RCP<Teuchos::ParameterList>                pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model )
  : model_(model)
{
  this->setParameterList(pList);

  inArgs_  = model_->createInArgs();
  outArgs_ = model_->createOutArgs();
  inArgs_  = model_->getNominalValues();
}

template<class Scalar>
void StepperForwardEuler<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  Teuchos::RCP<SolutionState<Scalar> > workingState = solutionHistory->getWorkingState();

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(workingState), std::logic_error,
    "Error - SolutionState, workingState, is invalid!\n");

  typedef Thyra::ModelEvaluatorBase MEB;
  const Scalar time = workingState->getTime();
  const Scalar dt   = workingState->getTimeStep();

  inArgs_.set_x(workingState->getX());
  if (inArgs_.supports(MEB::IN_ARG_t)) inArgs_.set_t(time);

  // For model evaluators whose state function f(x, x_dot, t) describes
  // an implicit ODE, and which accept an optional x_dot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // x_dot = f(x, t)
  if (inArgs_.supports(MEB::IN_ARG_x_dot)) inArgs_.set_x_dot(Teuchos::null);
  outArgs_.set_f(workingState->getXDot());

  model_->evalModel(inArgs_,outArgs_);

  // Forward Euler update, x = x + dt*xdot
  Thyra::Vp_StV(workingState->getX().ptr(),dt,*(workingState->getXDot().ptr()));

  workingState->stepperState_->stepperStatus_ = Status::PASSED;
  return;
}


/** \brief Provide a StepperState to the SolutionState.
 *  Forward Euler does not have any special state data,
 *  so just provide the base class StepperState with the
 *  ForwardEuler dsecription.  This can be checked to ensure
 *  that the input StepperState can be used by Forward Euler.
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
      << "model_ = " << model_->description() << std::endl;
}


template <class Scalar>
void StepperForwardEuler<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pList));
  //pList->validateParameters(*this->getValidParameters());
  pList_ = pList;

  std::string stepperType = pList_->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Forward Euler",
    std::logic_error,
       "Error - Stepper sublist is not 'Forward Euler'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");

  Teuchos::readVerboseObjectSublist(&*pList_,this);
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperForwardEuler<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<Teuchos::ParameterList> validPL;

  if (is_null(validPL)) {

    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    std::ostringstream tmp;
    tmp << "'Stepper Type' must be 'Forward Euler'.";
    pl->set("Stepper Type", "Forward Euler", tmp.str());

    validPL = pl;
  }
  return validPL;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperForwardEuler<Scalar>::getNonconstParameterList()
{
  return(pList_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperForwardEuler<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = pList_;
  pList_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // TEMPUS_STEPPERFORWARDEULER_IMPL_HPP
