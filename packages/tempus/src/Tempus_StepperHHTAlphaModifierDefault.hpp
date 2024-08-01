//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperHHTAlphaModifierDefault_hpp
#define Tempus_StepperHHTAlphaModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperHHTAlphaModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperHHTAlpha.hpp"

namespace Tempus {

/** \brief Default modifier for StepperHHTAlpha.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperHHTAlphaModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperHHTAlphaModifierDefault
  : virtual public Tempus::StepperHHTAlphaModifierBase<Scalar> {
 public:
  /// Constructor
  StepperHHTAlphaModifierDefault() {}

  /// Destructor
  virtual ~StepperHHTAlphaModifierDefault() {}

  /// Modify HHTAlpha Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperHHTAlpha<Scalar> > /* stepper */,
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

#endif  // Tempus_StepperHHTAlphaModifierDefault_hpp
