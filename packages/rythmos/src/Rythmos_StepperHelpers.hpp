
#ifndef RYTHMOS_STEPPER_HELPERS_HPP
#define RYTHMOS_STEPPER_HELPERS_HPP


#include "Rythmos_StepperBase.hpp"
#include "Teuchos_Assert.hpp"
#include "Thyra_AssertOp.hpp"


namespace Rythmos {


/** \brief Assert valid transient model for StepperBase.
 *
 */
template<class Scalar>
void assertValidModel(
  const StepperBase<Scalar>& stepper,
  const Thyra::ModelEvaluator<Scalar>& model
  )
{

  typedef Thyra::ModelEvaluatorBase MEB;

  TEUCHOS_ASSERT(stepper.acceptsModel());

  const MEB::InArgs<Scalar> inArgs = model.createInArgs();
  const MEB::OutArgs<Scalar> outArgs = model.createOutArgs();

  //TEUCHOS_ASSERT(inArgs.supports(MEB::IN_ARG_t));
  TEUCHOS_ASSERT(inArgs.supports(MEB::IN_ARG_x));
  TEUCHOS_ASSERT(outArgs.supports(MEB::OUT_ARG_f));
  
  if (stepper.isImplicit()) { // implicit stepper
    TEUCHOS_ASSERT( inArgs.supports(MEB::IN_ARG_x_dot) );
    TEUCHOS_ASSERT( inArgs.supports(MEB::IN_ARG_alpha) );
    TEUCHOS_ASSERT( inArgs.supports(MEB::IN_ARG_beta) );
    TEUCHOS_ASSERT( outArgs.supports(MEB::OUT_ARG_W) );
  } 
  //else { // explicit stepper
  //  TEUCHOS_ASSERT( !inArgs.supports(MEB::IN_ARG_x_dot) );
  //  TEUCHOS_ASSERT( !inArgs.supports(MEB::IN_ARG_alpha) );
  //  TEUCHOS_ASSERT( !inArgs.supports(MEB::IN_ARG_beta) );
  //  TEUCHOS_ASSERT( !outArgs.supports(MEB::OUT_ARG_W) );
  //}

}


/** \brief Set an initial condition on a stepper from a model if the stepper
 * does not already have an initial condition.
 *
 * \returns Returns <tt>true</tt> if the stepper->setInitialCondition(...) was
 * called and returns <tt>false</tt> otherwise.
 */
template<class Scalar>
bool setDefaultInitialConditionFromNominalValues(
  const Thyra::ModelEvaluator<Scalar>& model,
  const Ptr<StepperBase<Scalar> >& stepper
  )
{

  typedef ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;

  if (isInitialized(*stepper))
    return false;  // Already has an initial condition
  
  MEB::InArgs<Scalar> initCond = model.getNominalValues();

  if (!is_null(initCond.get_x())) {
    // IC has x, we will assume that initCont.get_t() is the valid start time.
    // Therefore, we just need to check that x_dot is also set or we will
    // create a zero x_dot
#ifdef HAVE_RYTHMOS_DEBUG
    THYRA_ASSERT_VEC_SPACES( "setInitialConditionIfExists(...)", 
      *model.get_x_space(), *initCond.get_x()->space() );
#endif
    if (initCond.supports(MEB::IN_ARG_x_dot)) {
      if (is_null(initCond.get_x_dot())) {
        const RCP<Thyra::VectorBase<Scalar> > x_dot =
          createMember(model.get_x_space());
        assign(x_dot.ptr(), ST::zero());
      }
      else {
#ifdef HAVE_RYTHMOS_DEBUG
        THYRA_ASSERT_VEC_SPACES( "setInitialConditionIfExists(...)", 
          *model.get_x_space(), *initCond.get_x_dot()->space() );
#endif
      }
    }
    stepper->setInitialCondition(initCond);
    return true;
  }

  // The model has not nominal values for which to set the initial
  // conditions so wo don't do anything!  The stepper will still have not
  return false;

}


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
#ifdef HAVE_RYTHMOS_DEBUG
  TEST_FOR_EXCEPT(0==stepper);
#endif // HAVE_RYTHMOS_DEBUG
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
