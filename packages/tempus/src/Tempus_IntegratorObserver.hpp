//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorObserver_hpp
#define Tempus_IntegratorObserver_hpp

#include "Tempus_config.hpp"
#include "Tempus_TimeStepControl.hpp"

// Forward declarations
namespace Tempus {
template <typename Scalar>
class Integrator;
}

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
template <class Scalar>
class IntegratorObserver {
 public:
  /// \name Basic IntegratorObserver methods
  //@{
  /// Observe the beginning of the time integrator.
  virtual void observeStartIntegrator(const Integrator<Scalar>& integrator) = 0;

  /// Observe the beginning of the time step loop.
  virtual void observeStartTimeStep(const Integrator<Scalar>& integrator) = 0;

  /// Observe after the next time step size is selected. The
  /// observer can choose to change the current integratorStatus.
  virtual void observeNextTimeStep(const Integrator<Scalar>& integrator) = 0;

  /// Observe before Stepper takes step.
  virtual void observeBeforeTakeStep(const Integrator<Scalar>& integrator) = 0;

  /// Observe after Stepper takes step.
  virtual void observeAfterTakeStep(const Integrator<Scalar>& integrator) = 0;

  /// Observe after checking time step. Observer can still fail the time step
  /// here.
  virtual void observeAfterCheckTimeStep(
      const Integrator<Scalar>& integrator) = 0;

  /// Observe the end of the time step loop.
  virtual void observeEndTimeStep(const Integrator<Scalar>& integrator) = 0;

  /// Observe the end of the time integrator.
  virtual void observeEndIntegrator(const Integrator<Scalar>& integrator) = 0;

  /// default destructor
  virtual ~IntegratorObserver() = default;
  //@}
};
}  // namespace Tempus
#endif  // Tempus_IntegratorObserver_hpp
