//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPI3ModifierX_hpp
#define Tempus_StepperEPI3ModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperEPI3ModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperEPI3.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperEPI3ModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperEPI3ModifierXDefault
  : virtual public Tempus::StepperEPI3ModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperEPI3ModifierXDefault() {}

  /// Destructor
  virtual ~StepperEPI3ModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
      Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */, const Scalar /* time */,
      const Scalar /* dt */,
      const typename StepperEPI3ModifierXBase<Scalar>::MODIFIER_TYPE
          modType)
  {
    switch (modType) {
      case StepperEPI3ModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperEPI3ModifierXBase<Scalar>::X_BEFORE_EXPLICIT_EVAL:
      case StepperEPI3ModifierXBase<Scalar>::XDOT_END_STEP: {
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

#endif  // Tempus_StepperEPI3ModifierX_hpp
