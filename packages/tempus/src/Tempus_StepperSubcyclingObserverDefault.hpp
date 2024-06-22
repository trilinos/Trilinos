//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperSubcyclingObserverDefault_hpp
#define Tempus_StepperSubcyclingObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperSubcyclingObserverBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperSubcycling.hpp"

namespace Tempus {

/** \brief Default observer for StepperSubcycling.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperSubcyclingObserverBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperSubcyclingObserverDefault
  : virtual public Tempus::StepperSubcyclingObserverBase<Scalar> {
 public:
  /// Constructor
  StepperSubcyclingObserverDefault() {}

  /// Destructor
  virtual ~StepperSubcyclingObserverDefault() {}

  /// Observe Subcycling Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperSubcycling<Scalar> > /* stepper */,
      const typename StepperSubcyclingAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperSubcyclingAppAction<Scalar>::BEGIN_STEP:
      case StepperSubcyclingAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperSubcyclingObserverDefault_hpp
