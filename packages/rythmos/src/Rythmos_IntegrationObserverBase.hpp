
#ifndef RYTHMOS_INTEGRATION_OBSERVER_BASE_HPP
#define RYTHMOS_INTEGRATION_OBSERVER_BASE_HPP

#include "Rythmos_Types.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"

namespace Rythmos {


// Forwards
template<class Scalar> class TimeRange;
template<class Scalar> class StepControlInfo;
template<class Scalar> class StepperBase;


/** \brief Base class for strategy objects that observe and time integration
 * by observing the stepper object.
 *
 * ToDo: Finish Implementation!
 */
template<class Scalar>
class IntegrationObserverBase
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<IntegrationObserverBase<Scalar> >
{
public:

  /** \brief Clone this integration observer if supported .
   *
   * Here, the cloned object just has to have the ouptut strategy copied, not
   * the complete state of the object mid way through an integration.
   */
  virtual RCP<IntegrationObserverBase<Scalar> >
  cloneIntegrationObserver() const = 0;

  /** \brief Reset the observer to prepair it to observe another integration.
   *
   * \param integrationTimeDomain [in] The time domain over which the
   * integration will be defined.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>integrationTimeDomain.length() > 0.0</tt>
   * </ul>
   *
   * \todo Add initial guess as an argument
   */
  virtual void resetIntegrationObserver(
    const TimeRange<Scalar> &integrationTimeDomain
    // ToDo: Pass in the initial condition to the observer
    ) = 0;

  /** \brief Observe the beginning of a time integration loop.
   *
   * \param stepper [in] The stepper object.
   *
   * Warning! This function is *NOT* stateless.  It should be called once and
   * only once at the beginning of getFwdPoints().
   *
   * NOTE: The function <tt>resetIntegrationControlStrategy()</tt> must be
   * called prior to even the first call to function.
   *
   * NOTE: This method should be pure virtual but has been given a
   * default implementation for backwards compatibility.  We will make
   * this pure virtual in the future.
   */
  virtual void observeStartTimeIntegration(
    const StepperBase<Scalar> &stepper);

  /** \brief Observe the end of a time integration loop.
   *
   * \param stepper [in] The stepper object.
   *
   * Warning! This function is *NOT* stateless.  It should be called once and
   * only once at the end of getFwdPoints().
   *
   * NOTE: The function <tt>resetIntegrationControlStrategy()</tt> must be
   * called prior to even the first call to function.
   *
   * NOTE: This method should be pure virtual but has been given a
   * default implementation for backwards compatibility.  We will make
   * this pure virtual in the future.
   */
  virtual void observeEndTimeIntegration(
    const StepperBase<Scalar> &stepper);

  // ToDo: add observeStartTimeStep(stepper, ...)
  /** \brief Observer the beginning of an integration step.
   *
   * \param stepper [in] The stepper object.
   *
   * \param stepCtrlInfo [in] The info for the time step about to be
   * taken.
   *
   * \param timeStepIter [in] The time step iteration counter.  In the first
   * call to this function, this should be <tt>timeStepIter==0</tt> and it
   * should be incremented on each call only once.  While the concrete
   * implementation of <tt>*this</tt> could keep track of the this counter,
   * putting it in the argument list helps to simplify logic and helps to
   * validate correct usage.
   *
   * Warning! This function is *NOT* stateless.  It should be called once and
   * only once at the beginning of each time step.
   *
   * NOTE: The function <tt>resetIntegrationControlStrategy()</tt> must be
   * called prior to even the first call to function.
   *
   * NOTE: This method should be pure virtual but has been given a
   * default implementation for backwards compatibility.  We will make
   * this pure virtual in the future.
   */
  virtual void observeStartTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  /** \brief Observe a successfully completed integration step.
   *
   * \param stepper [in] The stepper object that was just stepped forward once
   * to integrate the transient ODE/DAE equations.  On the very first call and
   * every other call, this stepper should have a finite time range for a
   * successfull step.
   *
   * \param stepCtrlInfo [in] The info for the actual time step that was just
   * completed.
   *
   * \param timeStepIter [in] The time step iteration counter.  In the first
   * call to this function, this should be <tt>timeStepIter==0</tt> and it
   * should be incremented on each call only once.  While the concrete
   * implementation of <tt>*this</tt> could keep track of the this counter,
   * putting it in the argument list helps to simplify logic and helps to
   * validate correct usage.
   *
   * Warning! This function is *NOT* stateless.  It should be called once and
   * only once per time step iteration.
   *
   * NOTE: The function <tt>resetIntegrationControlStrategy()</tt> must be
   * called prior to even the first call to function.
   *
   * NOTE: If <tt>isInitialTimeStep(stepper->getTimeRange(), fullTimeRange) ==
   * true</tt> then this is the first time step (where <tt>fullTimeRange</tt>
   * was passed into <tt>resetIntegrationObserver()</tt>.
   *
   * NOTE: If <tt>isFinalTimeStep(stepper->getTimeRange(), fullTimeRange) ==
   * true</tt> then this is the last time step (where <tt>fullTimeRange</tt>
   * was passed into <tt>resetIntegrationObserver()</tt>.
   */
  virtual void observeCompletedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    ) = 0;

  /** \brief Observer a failed integration step.
   *
   * \param stepper [in] The stepper object that was just stepped forward once
   * to integrate the transient ODE/DAE equations.  On the very first call and
   * every other call, this stepper should have a finite time range for a
   * successfull step.
   *
   * \param stepCtrlInfo [in] The info for the actual time step that was just
   * attempted.
   *
   * \param timeStepIter [in] The time step iteration counter.  In the first
   * call to this function, this should be <tt>timeStepIter==0</tt> and it
   * should be incremented on each call only once.  While the concrete
   * implementation of <tt>*this</tt> could keep track of the this counter,
   * putting it in the argument list helps to simplify logic and helps to
   * validate correct usage.
   *
   * Warning! This function is *NOT* stateless.  It should be called once and
   * only once per failed time step iteration.
   *
   * NOTE: The function <tt>resetIntegrationControlStrategy()</tt> must be
   * called prior to even the first call to function.
   *
   * NOTE: If <tt>isInitialTimeStep(stepper->getTimeRange(), fullTimeRange) ==
   * true</tt> then this is the first time step (where <tt>fullTimeRange</tt>
   * was passed into <tt>resetIntegrationObserver()</tt>.
   *
   * NOTE: If <tt>isFinalTimeStep(stepper->getTimeRange(), fullTimeRange) ==
   * true</tt> then this is the last time step (where <tt>fullTimeRange</tt>
   * was passed into <tt>resetIntegrationObserver()</tt>.
   *
   * NOTE: This method should be pure virtual but has been given a
   * default implementation for backwards compatibility.  We will make
   * this pure virtual in the future.
   */
  virtual void observeFailedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

};


/** \brief Helper function to tell an IntegrationObserverBase object if an
 * observed time step is the first time step or not.
 *
 * \todo Unit test this function!
 */
template<class Scalar>
bool isInitialTimeStep(
  const TimeRange<Scalar> &currentTimeRange,
  const TimeRange<Scalar> &fullTimeRange
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return compareTimeValues(currentTimeRange.lower(), fullTimeRange.lower()) == ST::zero();
}


/** \brief Helper function to tell an IntegrationObserverBase object if an
 * observed time step is the last time step or not.
 *
 * \todo Unit test this function!
 */
template<class Scalar>
bool isFinalTimeStep(
  const TimeRange<Scalar> &currentTimeRange,
  const TimeRange<Scalar> &fullTimeRange
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return compareTimeValues(currentTimeRange.upper(), fullTimeRange.upper()) >= ST::zero();
}


// /////////////////////////////////////////////////////////////
/*  Default implementations for backwards compatibility

   Some of the observer methods were added after rythmos was released.
   Even though all methods should be pure virtual, we provide a
   default implementation here for the recently added methods to
   maintain backwards compatibility.  These should be removed in the
   future.
*/

template<class Scalar>
void IntegrationObserverBase<Scalar>::
observeStartTimeIntegration(const StepperBase<Scalar> &/* stepper */)
{

}    

template<class Scalar>
void IntegrationObserverBase<Scalar>::
observeEndTimeIntegration(const StepperBase<Scalar> &/* stepper */)
{

}    

template<class Scalar>
void IntegrationObserverBase<Scalar>::
observeStartTimeStep(
    const StepperBase<Scalar> &/* stepper */,
    const StepControlInfo<Scalar> &/* stepCtrlInfo */,
    const int /* timeStepIter */
    )
{

}    

template<class Scalar>
void IntegrationObserverBase<Scalar>::
observeFailedTimeStep(
    const StepperBase<Scalar> &/* stepper */,
    const StepControlInfo<Scalar> &/* stepCtrlInfo */,
    const int /* timeStepIter */
    )
{

}    


} // namespace Rythmos


#endif // RYTHMOS_INTEGRATION_OBSERVER_BASE_HPP
