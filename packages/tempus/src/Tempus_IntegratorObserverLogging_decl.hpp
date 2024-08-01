//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorObserverLogging_decl_hpp
#define Tempus_IntegratorObserverLogging_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_IntegratorObserver.hpp"
#include <list>

namespace Tempus {

/** \brief This observer logs calls to observer functions.
 *  This observer simply logs and counts the calls to each of the
 *  observer functions.  This is useful in monirtoring and debugging
 *  the time integration.
 */
template <class Scalar>
class IntegratorObserverLogging
  : virtual public Tempus::IntegratorObserver<Scalar> {
 public:
  /// Constructor
  IntegratorObserverLogging();

  /// Destructor
  virtual ~IntegratorObserverLogging();

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
  //@}

  void resetLogCounters();

  Teuchos::RCP<const std::map<std::string, int> > getCounters();

  Teuchos::RCP<const std::list<std::string> > getOrder();

  /** \name String names logged in map
      Use these strings to validate a call stack with this observer
  */
  //@{
  const std::string nameObserveStartIntegrator_;
  const std::string nameObserveStartTimeStep_;
  const std::string nameObserveNextTimeStep_;
  const std::string nameObserveBeforeTakeStep_;
  const std::string nameObserveAfterTakeStep_;
  const std::string nameObserveAfterCheckTimeStep_;
  const std::string nameObserveEndTimeStep_;
  const std::string nameObserveEndIntegrator_;
  //@}

 private:
  /** \brief Asserts next call on the stack is correct and removes from stack

      This is a const method so that it can be called from the
      derived IntegratorObserver methods that are const.
  */
  void logCall(const std::string call) const;

  Teuchos::RCP<std::map<std::string, int> > counters_;
  Teuchos::RCP<std::list<std::string> > order_;
};

}  // namespace Tempus
#endif  // Tempus_IntegratorObserverLogging_decl_hpp
