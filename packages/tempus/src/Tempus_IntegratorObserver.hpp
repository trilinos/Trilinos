#ifndef TEMPUS_INTEGRATOROBSERVER_HPP
#define TEMPUS_INTEGRATOROBSERVER_HPP

#include "Tempus_SolutionHistory.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace tempus {

/** \brief IntegratorObserver class for time integrators.
 *
 * This is a means for application developers to perform tasks
 * during the time integrator, e.g.,
 *   - Compute specific quantities
 *   - Output information
 *   - Adjust the time step (CFL stability)
 *   - "Massage" the working solution state
 *   - ...
 *
 * <b>Design Considerations</b>
 *   - IntegratorObserver should have access to the entire SolutionHistory,
 *     as application developers may have that need.
 *   - The needed IntegratorObserver functions are determined by the
 *     access needs in Integrator::advanceTime().
 *   - IntegratorObserver is not stateless!  Developers may touch the
 *     solution state!  Developers need to be careful not to break the
 *     restart (checkpoint) capability.
 *   - The functions in this base class are simple no-op functions, as
 *     all basic functionality should be handled through other methods.
 */
template<class Scalar>
class IntegratorObserver
{
public:

  /// Constructor
  virtual IntegratorObserver(
    const RCP<SolutionHistory<Scalar> >& solutionHistory_,
    const RCP<TimeStepControl<Scalar> >& timeStepControl_)
    : solutionHistory(solutionHistory), timeStepControl(timeStepControl_)
  {}

  /// Destructor
  virtual ~IntegratorObserver();

  /// \name Basic IntegratorObserver methods
  //@{
    /// Observe the beginning of the time integrator.
    virtual void observeStartIntegrator(){}

    /// Observe the beginning of the time step loop.
    virtual void observeStartTimeStep(){}

    /// Observe after Integrator has declared failure.
    virtual void observeFailedIntegrator(){}

    /// Observe before Stepper takes step.
    virtual void observeBeforeTakeStep(){}

    /// Observe after Stepper has declared failure.
    virtual void observeFailedTimeStep(){}

    /// Observe after accepting time step.
    virtual void observeAcceptedTimeStep(){}

    /// Observe the end of the time integrator.
    virtual void observeEndIntegrator(){}
  //@}

protected:

  RCP<SolutionHistory<Scalar> > solutionHistory;
  RCP<TimeStepControl<Scalar> > timeStepControl;

};
} // namespace tempus
#endif // TEMPUS_INTEGRATOROBSERVER_HPP
