//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperTrapezoidalModifierX_hpp
#define Tempus_StepperTrapezoidalModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperTrapezoidalModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperTrapezoidal.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperTrapezoidalModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperTrapezoidalModifierXDefault
  : virtual public Tempus::StepperTrapezoidalModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperTrapezoidalModifierXDefault() {}

  /// Destructor
  virtual ~StepperTrapezoidalModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperTrapezoidalModifierXBase<Scalar>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperTrapezoidalModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperTrapezoidalModifierXBase<Scalar>::X_BEFORE_SOLVE:
      case StepperTrapezoidalModifierXBase<Scalar>::X_AFTER_SOLVE:
      case StepperTrapezoidalModifierXBase<Scalar>::XDOT_END_STEP: {
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

#endif  // Tempus_StepperTrapezoidalModifierX_hpp
