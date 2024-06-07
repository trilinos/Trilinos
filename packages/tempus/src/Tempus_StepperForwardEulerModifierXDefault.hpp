//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperForwardEulerModifierX_hpp
#define Tempus_StepperForwardEulerModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperForwardEulerModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperForwardEuler.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperForwardEulerModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperForwardEulerModifierXDefault
  : virtual public Tempus::StepperForwardEulerModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperForwardEulerModifierXDefault() {}

  /// Destructor
  virtual ~StepperForwardEulerModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperForwardEulerModifierXBase<Scalar>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperForwardEulerModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperForwardEulerModifierXBase<Scalar>::X_BEFORE_EXPLICIT_EVAL:
      case StepperForwardEulerModifierXBase<Scalar>::XDOT_END_STEP: {
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

#endif  // Tempus_StepperForwardEulerModifierX_hpp
