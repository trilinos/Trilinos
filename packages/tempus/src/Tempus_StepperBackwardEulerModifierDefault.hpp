//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBackwardEulerModifierDefault_hpp
#define Tempus_StepperBackwardEulerModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperBackwardEulerModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperBackwardEuler.hpp"

namespace Tempus {

/** \brief Default modifier for StepperBackwardEuler.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperBackwardEulerModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperBackwardEulerModifierDefault
  : virtual public Tempus::StepperBackwardEulerModifierBase<Scalar> {
 public:
  /// Constructor
  StepperBackwardEulerModifierDefault() {}

  /// Destructor
  virtual ~StepperBackwardEulerModifierDefault() {}

  /// Modify BackwardEuler Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperBackwardEuler<Scalar> > /* stepper */,
      const typename StepperBackwardEulerAppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperBackwardEulerAppAction<Scalar>::BEGIN_STEP:
      case StepperBackwardEulerAppAction<Scalar>::BEFORE_SOLVE:
      case StepperBackwardEulerAppAction<Scalar>::AFTER_SOLVE:
      case StepperBackwardEulerAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperBackwardEulerModifierDefault_hpp
