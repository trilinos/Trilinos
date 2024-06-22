//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitAFormModifierX_hpp
#define Tempus_StepperNewmarkImplicitAFormModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperNewmarkImplicitAFormModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperNewmarkImplicitAForm.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperNewmarkImplicitAFormModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperNewmarkImplicitAFormModifierXDefault
  : virtual public Tempus::StepperNewmarkImplicitAFormModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperNewmarkImplicitAFormModifierXDefault() {}

  /// Destructor
  virtual ~StepperNewmarkImplicitAFormModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
                      const Scalar /* time */, const Scalar /* dt */,
                      const typename StepperNewmarkImplicitAFormModifierXBase<
                          Scalar>::MODIFIER_TYPE modType)
  {
    switch (modType) {
      case StepperNewmarkImplicitAFormModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperNewmarkImplicitAFormModifierXBase<Scalar>::X_BEFORE_SOLVE:
      case StepperNewmarkImplicitAFormModifierXBase<Scalar>::X_AFTER_SOLVE:
      case StepperNewmarkImplicitAFormModifierXBase<Scalar>::X_END_STEP: {
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

#endif  // Tempus_StepperNewmarkImplicitAFormModifierX_hpp
