//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPIModifierX_hpp
#define Tempus_StepperEPIModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperEPIModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperEPI.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperEPIModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperEPIModifierXDefault
  : virtual public Tempus::StepperEPIModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperEPIModifierXDefault() {}

  /// Destructor
  virtual ~StepperEPIModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperEPIModifierXBase<Scalar>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperEPIModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperEPIModifierXBase<Scalar>::X_BEFORE_EXPLICIT_EVAL:
      case StepperEPIModifierXBase<Scalar>::XDOT_END_STEP: {
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

#endif  // Tempus_StepperEPIModifierX_hpp
