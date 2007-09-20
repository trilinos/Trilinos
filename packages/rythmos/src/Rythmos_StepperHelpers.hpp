
#ifndef RYTHMOS_STEPPER_HELPERS_HPP
#define RYTHMOS_STEPPER_HELPERS_HPP


#include "Rythmos_StepperBase.hpp"


namespace Rythmos {


/** \brief Restart a time stepper.
 *
 * This simple helper function just grabs the state out of the
 * <tt>*stepper</tt> object and then resets it on itself as an initial
 * condition.  This is generally used to restart a stepper when passing over a
 * breakpoint where the model is expected to be discontinuous in some way.
 */
template<class Scalar>
void restart( StepperBase<Scalar> *stepper )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(0==stepper);
#endif // TEUCHOS_DEBUG
  typedef Thyra::ModelEvaluatorBase MEB;
  const Rythmos::StepStatus<double>
    stepStatus = stepper->getStepStatus();
  const RCP<const Thyra::ModelEvaluator<Scalar> >
    model = stepper->getModel();
  // First, copy all of the model's state, including parameter values etc.
  MEB::InArgs<double> initialCondition = model->createInArgs();
  initialCondition.setArgs(model->getNominalValues());
  // Set the current values of the state and time
  RCP<const Thyra::VectorBase<double> > x, x_dot;
  Rythmos::get_x_and_x_dot(*stepper,stepStatus.time,&x,&x_dot);
  initialCondition.set_x(x);
  initialCondition.set_x_dot(x_dot);
  initialCondition.set_t(stepStatus.time);
  // Set the new initial condition back on the stepper.  This will effectively
  // reset the stepper to think that it is starting over again (which it is).
  stepper->setInitialCondition(initialCondition);
}



} // namespace Rythmos


#endif // RYTHMOS_STEPPER_HELPERS_HPP
