//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperRKModifierDefault_hpp
#define Tempus_StepperRKModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperRKModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperRKBase.hpp"

namespace Tempus {

/** \brief Default modifier for StepperRK.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperRKModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperRKModifierDefault
  : virtual public Tempus::StepperRKModifierBase<Scalar> {
 public:
  /// Constructor
  StepperRKModifierDefault() {}

  /// Destructor
  virtual ~StepperRKModifierDefault() {}

  /// Modify RK Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperRKBase<Scalar> > /* stepper */,
      const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperRKAppAction<Scalar>::BEGIN_STEP:
      case StepperRKAppAction<Scalar>::BEGIN_STAGE:
      case StepperRKAppAction<Scalar>::BEFORE_SOLVE:
      case StepperRKAppAction<Scalar>::AFTER_SOLVE:
      case StepperRKAppAction<Scalar>::BEFORE_EXPLICIT_EVAL:
      case StepperRKAppAction<Scalar>::END_STAGE:
      case StepperRKAppAction<Scalar>::END_STEP: {
        // No-op.
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error,
            "Error - unknown action location = " + std::to_string(actLoc) +
                "\n"
                "  Valid actions are\n"
                "    StepperRKAppAction<Scalar>::BEGIN_STEP           = " +
                std::to_string(StepperRKAppAction<Scalar>::BEGIN_STEP) +
                "\n"
                "    StepperRKAppAction<Scalar>::BEGIN_STAGE          = " +
                std::to_string(StepperRKAppAction<Scalar>::BEGIN_STAGE) +
                "\n"
                "    StepperRKAppAction<Scalar>::BEFORE_SOLVE         = " +
                std::to_string(StepperRKAppAction<Scalar>::BEFORE_SOLVE) +
                "\n"
                "    StepperRKAppAction<Scalar>::AFTER_SOLVE          = " +
                std::to_string(StepperRKAppAction<Scalar>::AFTER_SOLVE) +
                "\n"
                "    StepperRKAppAction<Scalar>::BEFORE_EXPLICIT_EVAL = " +
                std::to_string(
                    StepperRKAppAction<Scalar>::BEFORE_EXPLICIT_EVAL) +
                "\n"
                "    StepperRKAppAction<Scalar>::END_STAGE            = " +
                std::to_string(StepperRKAppAction<Scalar>::END_STAGE) +
                "\n"
                "    StepperRKAppAction<Scalar>::END_STEP             = " +
                std::to_string(StepperRKAppAction<Scalar>::END_STEP) + "\n");
    }
  }
};

}  // namespace Tempus

#endif  // Tempus_StepperRKModifierDefault_hpp
