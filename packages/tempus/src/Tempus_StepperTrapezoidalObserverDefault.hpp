//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperTrapezoidalObserverDefault_hpp
#define Tempus_StepperTrapezoidalObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperTrapezoidalObserverBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperTrapezoidal.hpp"

namespace Tempus {

/** \brief Default observer for StepperTrapezoidal.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperTrapezoidalObserverBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperTrapezoidalObserverDefault
  : virtual public Tempus::StepperTrapezoidalObserverBase<Scalar> {
 public:
  /// Constructor
  StepperTrapezoidalObserverDefault() {}

  /// Destructor
  virtual ~StepperTrapezoidalObserverDefault() {}

  /// Observe Subcycling Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperTrapezoidal<Scalar> > /* stepper */,
      const typename StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperTrapezoidalAppAction<Scalar>::BEGIN_STEP:
      case StepperTrapezoidalAppAction<Scalar>::BEFORE_SOLVE:
      case StepperTrapezoidalAppAction<Scalar>::AFTER_SOLVE:
      case StepperTrapezoidalAppAction<Scalar>::END_STEP: {
        // No-op.
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown action location.\n");
    }
  }
};

}  // namespace Tempus

#endif  // Tempus_StepperTrapezoidalObserverDefault_hpp
