//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperRKObserverDefault_hpp
#define Tempus_StepperRKObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperRKObserverBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperRKBase.hpp"

namespace Tempus {

/** \brief Default observer for StepperRK.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperRKObserverBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperRKObserverDefault
  : virtual public Tempus::StepperRKObserverBase<Scalar> {
 public:
  /// Constructor
  StepperRKObserverDefault() {}

  /// Destructor
  virtual ~StepperRKObserverDefault() {}

  /// Observe RK Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperRKBase<Scalar> > /* stepper */,
      const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperRKAppAction<Scalar>::BEGIN_STEP:
      case StepperRKAppAction<Scalar>::BEFORE_SOLVE:
      case StepperRKAppAction<Scalar>::AFTER_SOLVE:
      case StepperRKAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperRKObserverDefault_hpp
