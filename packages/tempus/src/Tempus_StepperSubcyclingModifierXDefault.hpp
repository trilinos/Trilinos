//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperSubcyclingModifierX_hpp
#define Tempus_StepperSubcyclingModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperSubcyclingModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperSubcycling.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperSubcyclingModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperSubcyclingModifierXDefault
  : virtual public Tempus::StepperSubcyclingModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperSubcyclingModifierXDefault() {}

  /// Destructor
  virtual ~StepperSubcyclingModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperSubcyclingModifierXBase<Scalar>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperSubcyclingModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperSubcyclingModifierXBase<Scalar>::XDOT_END_STEP: {
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

#endif  // Tempus_StepperSubcyclingModifierX_hpp
