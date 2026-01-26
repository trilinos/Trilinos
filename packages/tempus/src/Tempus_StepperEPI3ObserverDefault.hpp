//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPI3ObserverDefault_hpp
#define Tempus_StepperEPI3ObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperEPI3ObserverBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperEPI3.hpp"

namespace Tempus {

/** \brief Default observer for StepperEPI3.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperEPI3ObserverBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperEPI3ObserverDefault
  : virtual public Tempus::StepperEPI3ObserverBase<Scalar> {
 public:
  /// Constructor
  StepperEPI3ObserverDefault() {}

  /// Destructor
  virtual ~StepperEPI3ObserverDefault() {}

  /// Observe EPI3 Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperEPI3<Scalar> > /* stepper */,
      const typename StepperEPI3AppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperEPI3AppAction<Scalar>::BEGIN_STEP:
      case StepperEPI3AppAction<Scalar>::BEFORE_EXPLICIT_EVAL:
      case StepperEPI3AppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperEPI3ObserverDefault_hpp
