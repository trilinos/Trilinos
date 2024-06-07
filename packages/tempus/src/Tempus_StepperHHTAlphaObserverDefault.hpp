//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperHHTAlphaObserverDefault_hpp
#define Tempus_StepperHHTAlphaObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperHHTAlphaObserverBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperHHTAlpha.hpp"

namespace Tempus {

/** \brief Default observer for StepperHHTAlpha.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperHHTAlphaObserverBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperHHTAlphaObserverDefault
  : virtual public Tempus::StepperHHTAlphaObserverBase<Scalar> {
 public:
  /// Constructor
  StepperHHTAlphaObserverDefault() {}

  /// Destructor
  virtual ~StepperHHTAlphaObserverDefault() {}

  /// Observe HHTAlpha Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperHHTAlpha<Scalar> > /* stepper */,
      const typename StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperHHTAlphaAppAction<Scalar>::BEGIN_STEP:
      case StepperHHTAlphaAppAction<Scalar>::BEFORE_SOLVE:
      case StepperHHTAlphaAppAction<Scalar>::AFTER_SOLVE:
      case StepperHHTAlphaAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperHHTAlphaObserverDefault_hpp
