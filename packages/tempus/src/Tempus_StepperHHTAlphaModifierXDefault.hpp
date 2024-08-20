//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperHHTAlphaModifierX_hpp
#define Tempus_StepperHHTAlphaModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperHHTAlphaModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperHHTAlpha.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperHHTAlphaModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperHHTAlphaModifierXDefault
  : virtual public Tempus::StepperHHTAlphaModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperHHTAlphaModifierXDefault() {}

  /// Destructor
  virtual ~StepperHHTAlphaModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperHHTAlphaModifierXBase<Scalar>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperHHTAlphaModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperHHTAlphaModifierXBase<Scalar>::X_BEFORE_SOLVE:
      case StepperHHTAlphaModifierXBase<Scalar>::X_AFTER_SOLVE:
      case StepperHHTAlphaModifierXBase<Scalar>::X_END_STEP: {
        // No-op.
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown modifier type.\n");
    }
  }
};

}  // namespace Tempus

#endif  // Tempus_StepperHHTAlphaModifierX_hpp
