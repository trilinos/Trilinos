//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperForwardEulerObserverDefault_hpp
#define Tempus_StepperForwardEulerObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperForwardEulerObserverBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperForwardEuler.hpp"

namespace Tempus {

/** \brief Default observer for StepperForwardEuler.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperForwardEulerObserverBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperForwardEulerObserverDefault
  : virtual public Tempus::StepperForwardEulerObserverBase<Scalar> {
 public:
  /// Constructor
  StepperForwardEulerObserverDefault() {}

  /// Destructor
  virtual ~StepperForwardEulerObserverDefault() {}

  /// Observe ForwardEuler Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperForwardEuler<Scalar> > /* stepper */,
      const typename StepperForwardEulerAppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperForwardEulerAppAction<Scalar>::BEGIN_STEP:
      case StepperForwardEulerAppAction<Scalar>::BEFORE_EXPLICIT_EVAL:
      case StepperForwardEulerAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperForwardEulerObserverDefault_hpp
