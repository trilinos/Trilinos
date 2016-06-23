#ifndef TEMPUS_INTEGRATOROBSERVER_HPP
#define TEMPUS_INTEGRATOROBSERVER_HPP

#include "Tempus_TimeStepControl.hpp"

namespace Tempus {

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
  IntegratorObserver(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory,
    const Teuchos::RCP<TimeStepControl<Scalar> >& timeStepControl)
    : solutionHistory_(solutionHistory), timeStepControl_(timeStepControl)
  {}

  /// Destructor
  virtual ~IntegratorObserver(){}

  /// \name Basic IntegratorObserver methods
  //@{
    /// Observe the beginning of the time integrator.
    virtual void observeStartIntegrator(){}

    /// Observe the beginning of the time step loop.
    virtual void observeStartTimeStep(){}

    /// Observe after the next time step size is selected.
    virtual void observeNextTimeStep(Status & integratorStatus){}

    /// Observe before Stepper takes step.
    virtual void observeBeforeTakeStep(){}

    /// Observe after Stepper takes step.
    virtual void observeAfterTakeStep(){}

    /// Observe after accepting time step.
    virtual void observeAcceptedTimeStep(Status & integratorStatus){}

    /// Observe the end of the time integrator.
    virtual void observeEndIntegrator(const Status integratorStatus){}
  //@}

protected:

  Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory_;
  Teuchos::RCP<TimeStepControl<Scalar> > timeStepControl_;

};
} // namespace Tempus
#endif // TEMPUS_INTEGRATOROBSERVER_HPP
