//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperTrapezoidalModifierDefault_hpp
#define Tempus_StepperTrapezoidalModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperTrapezoidalModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperTrapezoidal.hpp"

namespace Tempus {

/** \brief Default modifier for StepperTrapezoidal.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperTrapezoidalModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperTrapezoidalModifierDefault
  : virtual public Tempus::StepperTrapezoidalModifierBase<Scalar> {
 public:
  /// Constructor
  StepperTrapezoidalModifierDefault() {}

  /// Destructor
  virtual ~StepperTrapezoidalModifierDefault() {}

  /// Modify Trapezoidal Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperTrapezoidal<Scalar> > /* stepper */,
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

#endif  // Tempus_StepperTrapezoidalModifierDefault_hpp
