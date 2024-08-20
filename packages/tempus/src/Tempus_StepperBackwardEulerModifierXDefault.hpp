//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBackwardEulerModifierX_hpp
#define Tempus_StepperBackwardEulerModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperBackwardEulerModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperBackwardEuler.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperBackwardEulerModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperBackwardEulerModifierXDefault
  : virtual public Tempus::StepperBackwardEulerModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperBackwardEulerModifierXDefault() {}

  /// Destructor
  virtual ~StepperBackwardEulerModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperBackwardEulerModifierXBase<Scalar>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperBackwardEulerModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperBackwardEulerModifierXBase<Scalar>::X_BEFORE_SOLVE:
      case StepperBackwardEulerModifierXBase<Scalar>::X_AFTER_SOLVE:
      case StepperBackwardEulerModifierXBase<Scalar>::XDOT_END_STEP: {
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

#endif  // Tempus_StepperBackwardEulerModifierX_hpp
