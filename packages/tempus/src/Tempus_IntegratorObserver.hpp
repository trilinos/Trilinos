// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorObserver_hpp
#define Tempus_IntegratorObserver_hpp

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
 */
template<class Scalar>
class IntegratorObserver
{
public:

  /// \name Basic IntegratorObserver methods
  //@{
    /// Observe the beginning of the time integrator.
    virtual void observeStartIntegrator() = 0;

    /// Observe the beginning of the time step loop.
    virtual void observeStartTimeStep() = 0;

    /// Observe after the next time step size is selected.
    virtual void observeNextTimeStep(Status & integratorStatus) = 0;

    /// Observe before Stepper takes step.
    virtual void observeBeforeTakeStep() = 0;

    /// Observe after Stepper takes step.
    virtual void observeAfterTakeStep() = 0;

    /// Observe after accepting time step.
    virtual void observeAcceptedTimeStep(Status & integratorStatus) = 0;

    /// Observe the end of the time integrator.
    virtual void observeEndIntegrator(const Status integratorStatus) = 0;

    virtual void setSolutionHistory(
      Teuchos::RCP<SolutionHistory<Scalar> > sh) = 0;

    virtual void setTimeStepControl(
      Teuchos::RCP<TimeStepControl<Scalar> > tsc) = 0;
  //@}

};
} // namespace Tempus
#endif // Tempus_IntegratorObserver_hpp
