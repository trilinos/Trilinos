//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorObserverSubcycling_decl_hpp
#define Tempus_IntegratorObserverSubcycling_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_IntegratorObserver.hpp"
#include "Tempus_Integrator.hpp"
#include "Teuchos_Time.hpp"

namespace Tempus {

/** \brief IntegratorObserverSubcycling class for time integrators.
 *  This basic class has simple no-op functions, as all basic
 *  functionality should be handled through other methods.
 */
template <class Scalar>
class IntegratorObserverSubcycling
  : virtual public Tempus::IntegratorObserver<Scalar> {
 public:
  /// Constructor
  IntegratorObserverSubcycling();

  /// Destructor
  virtual ~IntegratorObserverSubcycling();

  /// \name Subcycling IntegratorObserver methods
  //@{
  /// Observe the beginning of the time integrator.
  virtual void observeStartIntegrator(
      const Integrator<Scalar>& integrator) override;

  /// Observe the beginning of the time step loop.
  virtual void observeStartTimeStep(
      const Integrator<Scalar>& integrator) override;

  /// Observe after the next time step size is selected.
  virtual void observeNextTimeStep(
      const Integrator<Scalar>& integrator) override;

  /// Observe before Stepper takes step.
  virtual void observeBeforeTakeStep(
      const Integrator<Scalar>& integrator) override;

  /// Observe after Stepper takes step.
  virtual void observeAfterTakeStep(
      const Integrator<Scalar>& integrator) override;

  /// Observe after checking time step.  Observer can still fail the time step
  /// here.
  virtual void observeAfterCheckTimeStep(
      const Integrator<Scalar>& integrator) override;

  /// Observe the end of the time step loop.
  virtual void observeEndTimeStep(
      const Integrator<Scalar>& integrator) override;

  /// Observe the end of the time integrator.
  virtual void observeEndIntegrator(
      const Integrator<Scalar>& integrator) override;
  //@}
};
}  // namespace Tempus
#endif  // Tempus_IntegratorObserverSubcycling_decl_hpp
