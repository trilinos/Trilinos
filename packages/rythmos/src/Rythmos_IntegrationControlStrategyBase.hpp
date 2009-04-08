
#ifndef RYTHMOS_INTEGRATION_CONTROL_STRATEGY_BASE_HPP
#define RYTHMOS_INTEGRATION_CONTROL_STRATEGY_BASE_HPP

#include "Rythmos_Types.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"

namespace Rythmos {


// Forwards
template<class Scalar> class StepControlInfo;
template<class Scalar> class TimeRange;
template<class Scalar> class StepperBase;


/** \brief Base class for strategy objects that control integration by
 * selecting step sizes for a stepper.
 *
 * ToDo: Finish Implementation!
 */
template<class Scalar>
class IntegrationControlStrategyBase
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<IntegrationControlStrategyBase<Scalar> >,
    virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \brief Clone this integration control object if supported .
   *
   * Here, the cloned object just has to have the control information copied,
   * not the complete state of the object mid way through an integration.
   */
  virtual RCP<IntegrationControlStrategyBase<Scalar> >
  cloneIntegrationControlStrategy() const = 0;

  /** \brief Reset the control algorithm to the beginning to start a new
   * integration.
   *
   * \param integrationTimeDomain [in] The time domain over which the
   * integration will be defined.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>integrationTimeDomain.length() > 0.0</tt>
   * </ul>
   */
  virtual void resetIntegrationControlStrategy(
    const TimeRange<Scalar> &integrationTimeDomain
    ) = 0;

  /** \brief Select the next time step control info.
   *
   * \param stepper [in] The stepper object that is being stepped forward in
   * time to integrate the transient ODE/DAE equations.  On the very first
   * call, this stepper should just have the initial condition.
   *
   * \param stepCtrlInfoLast [in] The actual time step that was taken on the
   * last time step.
   *
   * \param timeStepIter [in] The (zero-based) time step iteration counter.
   * In the first call to this function, this should be passed as
   * <tt>timeStepIter==0</tt> and it should be incremented on each call only
   * once.  While the concrete implementation if <tt>*this</tt> could keep
   * track of the this counter, putting it in the argument list helps to
   * simplify logic and helps to validate correct usage.
   *
   * \returns Returns the control information about the next step, including
   * if <tt>returnVal.stepSize</tt> is limited by a breakpoint and if that
   * breakpoint requires a restart of the stepper.  If no more steps are to be
   * taken then a step size of <tt>returnVal.stepSize < 0.0</tt> will be
   * returned and the time integrator client should halt the integration
   * immediately!
   *
   * Warning! This function is *NOT* stateless.  It should be called once and
   * only once per time step iteration.
   *
   * NOTE: The function <tt>resetIntegrationControlStrategy()</tt> must be
   * called prior to even the first call to function.
   */
  virtual StepControlInfo<Scalar>
  getNextStepControlInfo(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfoLast,
    const int timeStepIter
    ) = 0;

};


} // namespace Rythmos


#endif // RYTHMOS_INTEGRATION_CONTROL_STRATEGY_BASE_HPP
