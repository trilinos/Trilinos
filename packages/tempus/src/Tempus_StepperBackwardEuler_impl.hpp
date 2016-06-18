#ifndef TEMPUS_STEPPERBACKWARDEULER_IMPL_HPP
#define TEMPUS_STEPPERBACKWARDEULER_IMPL_HPP

// Teuchos
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;

namespace Tempus {


template <typename Scalar>
std::function<void (const Thyra::VectorBase<Scalar> &,
                          Thyra::VectorBase<Scalar> &)>
StepperBackwardEuler<Scalar>::xDotFunction(
  Scalar dt, Teuchos::RCP<const Thyra::VectorBase<Scalar> > x_old)
{
  return [=](const Thyra::VectorBase<Scalar> & x,
                   Thyra::VectorBase<Scalar> & x_dot)
    {
      // this is the Euler x dot vector
      Thyra::V_StVpStV(ptrFromRef(x_dot),1.0/dt,x,-1.0/dt,*x_old);
    };
}


// StepperBackwardEuler definitions:
template<class Scalar>
StepperBackwardEuler<Scalar>::StepperBackwardEuler(
  RCP<ParameterList>                         pList,
  const RCP<Thyra::ModelEvaluator<Scalar> >& transientModel )
{
  residualModel_ = rcp(new ResidualModelEvaluator<Scalar>(transientModel));

  inArgs_  = residualModel_->createInArgs();
  outArgs_ = residualModel_->createOutArgs();
  inArgs_  = residualModel_->getNominalValues();

  RCP<ParameterList> noxPL = Teuchos::sublist(pList,"NOX",true);
  solver_ = rcp(new Thyra::NOXNonlinearSolver());
  solver_->setParameterList(noxPL);
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::takeStep(
  const RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  RCP<SolutionState<Scalar> > workingState = solutionHistory->getWorkingState();

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

  residualModel_->evalModel(inArgs_,outArgs_);

  // Backward Euler update, x = x + dt*xdot
  Thyra::Vp_StV(workingState->getX().ptr(),dt,*(workingState->getXDot().ptr()));

  workingState->stepperState->stepperStatus = Status::PASSED;
  return;
}

/** \brief Provide a StepperState to the SolutionState.
 *  Backward Euler does not have any special state data,
 *  so just provide the base class StepperState with the
 *  BackwardEuler dsecription.  This can be checked to ensure
 *  that the input StepperState can be used by Backward Euler.
 */
template<class Scalar>
RCP<Tempus::StepperState<Scalar> > StepperBackwardEuler<Scalar>::
getDefaultStepperState()
{
  RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperBackwardEuler<Scalar>::description() const
{
  std::string name = "Backward Euler";
  return(name);
}


template<class Scalar>
void StepperBackwardEuler<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "residualModel_ = " << residualModel_->description() << std::endl;
}


template <class Scalar>
void StepperBackwardEuler<Scalar>::setParameterList(
  RCP<ParameterList> const& pList)
{
  TEUCHOS_TEST_FOR_EXCEPT(!is_null(pList));
  pList->validateParameters(*this->getValidParameters());
  pList_ = pList;

  Teuchos::readVerboseObjectSublist(&*pList_,this);
}


template<class Scalar>
RCP<const ParameterList> StepperBackwardEuler<Scalar>::getValidParameters() const
{
  static RCP<ParameterList> validPL;
  return validPL;
}


template <class Scalar>
RCP<ParameterList>
StepperBackwardEuler<Scalar>::getNonconstParameterList()
{
  return(pList_);
}


template <class Scalar>
RCP<ParameterList> StepperBackwardEuler<Scalar>::unsetParameterList()
{
  RCP<ParameterList> temp_plist = pList_;
  pList_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // TEMPUS_STEPPERBACKWARDEULER_IMPL_HPP
