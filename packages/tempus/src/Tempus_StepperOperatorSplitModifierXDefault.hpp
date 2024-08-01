//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperOperatorSplitModifierXDefault_hpp
#define Tempus_StepperOperatorSplitModifierXDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperOperatorSplitModifierXBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperOperatorSplit.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperOperatorSplit.
 *
 *  The default ModifierX provides no-op functionality for ModifierX.
 *  See StepperOperatorSplitModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperOperatorSplitModifierXDefault
  : virtual public Tempus::StepperOperatorSplitModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperOperatorSplitModifierXDefault() {}

  /// Destructor
  virtual ~StepperOperatorSplitModifierXDefault() {}

  /// Modify OperatorSplit Stepper.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperOperatorSplitModifierXBase<Scalar>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperOperatorSplitModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperOperatorSplitModifierXBase<Scalar>::X_BEFORE_STEPPER:
      case StepperOperatorSplitModifierXBase<Scalar>::X_AFTER_STEPPER:
      case StepperOperatorSplitModifierXBase<Scalar>::XDOT_END_STEP: {
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

#endif  // Tempus_StepperOperatorSplitModifierXDefault_hpp
