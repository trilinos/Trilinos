//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperOperatorSplitModifierDefault_hpp
#define Tempus_StepperOperatorSplitModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperOperatorSplitModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperOperatorSplit.hpp"

namespace Tempus {

/** \brief Default modifier for StepperOperatorSplit.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperOperatorSplitModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperOperatorSplitModifierDefault
  : virtual public Tempus::StepperOperatorSplitModifierBase<Scalar> {
 public:
  /// Constructor
  StepperOperatorSplitModifierDefault() {}

  /// Destructor
  virtual ~StepperOperatorSplitModifierDefault() {}

  /// Modify OperatorSplit Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperOperatorSplit<Scalar> > /* stepper */,
      const typename StepperOperatorSplitAppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    switch (actLoc) {
      case StepperOperatorSplitAppAction<Scalar>::BEGIN_STEP:
      case StepperOperatorSplitAppAction<Scalar>::BEFORE_STEPPER:
      case StepperOperatorSplitAppAction<Scalar>::AFTER_STEPPER:
      case StepperOperatorSplitAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperOperatorSplitModifierDefault_hpp
