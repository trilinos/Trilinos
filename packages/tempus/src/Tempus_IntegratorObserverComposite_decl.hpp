//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorObserverComposite_decl_hpp
#define Tempus_IntegratorObserverComposite_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_IntegratorObserver.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This observer is a composite observer,
 *
 *  which takes other IntegratorObservers and sequentially calls each
 *  individual observer function.
 */
template <class Scalar>
class IntegratorObserverComposite
  : virtual public Tempus::IntegratorObserver<Scalar> {
 public:
  /// Default constructor
  IntegratorObserverComposite();

  /// Destructor
  virtual ~IntegratorObserverComposite();

  /// \name Override IntegratorObserver basic methods
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

  /// Observe after checking time step.
  virtual void observeAfterCheckTimeStep(
      const Integrator<Scalar>& integrator) override;

  /// Observe the end of the time step loop.
  virtual void observeEndTimeStep(
      const Integrator<Scalar>& integrator) override;

  /// Observe the end of the time integrator.
  virtual void observeEndIntegrator(
      const Integrator<Scalar>& integrator) override;

  // add observer to the composite observer list
  void addObserver(const Teuchos::RCP<IntegratorObserver<Scalar> >& observer);

  // clear all observer from the composite observer list
  void clearObservers();
  //@}

 private:
  std::vector<Teuchos::RCP<IntegratorObserver<Scalar> > > observers_;
};

}  // namespace Tempus
#endif  // Tempus_IntegratorObserverComposite_decl_hpp
