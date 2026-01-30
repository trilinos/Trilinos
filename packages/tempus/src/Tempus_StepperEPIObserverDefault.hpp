//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPIObserverDefault_hpp
#define Tempus_StepperEPIObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperEPIObserverBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperEPI.hpp"

namespace Tempus {

/** \brief Default observer for StepperEPI.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperEPIObserverBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperEPIObserverDefault
  : virtual public Tempus::StepperEPIObserverBase<Scalar> {
 public:
  /// Constructor
  StepperEPIObserverDefault() {}

  /// Destructor
  virtual ~StepperEPIObserverDefault() {}

  /// Observe EPI Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperEPI<Scalar> > /* stepper */,
      const typename StepperEPIAppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperEPIAppAction<Scalar>::BEGIN_STEP:
      case StepperEPIAppAction<Scalar>::BEFORE_EXPLICIT_EVAL:
      case StepperEPIAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperEPIObserverDefault_hpp
