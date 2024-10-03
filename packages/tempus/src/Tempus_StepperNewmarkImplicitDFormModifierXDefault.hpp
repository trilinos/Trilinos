//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitDFormModifierX_hpp
#define Tempus_StepperNewmarkImplicitDFormModifierX_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperNewmarkImplicitDFormModifierXBase.hpp"

namespace Tempus {

/** \brief Default ModifierX for StepperNewmarkImplicitDForm.
 *
 *  The default provides no-op functionality for ModifierX.
 *  See StepperNewmarkImplicitDFormModifierXBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperNewmarkImplicitDFormModifierXDefault
  : virtual public Tempus::StepperNewmarkImplicitDFormModifierXBase<Scalar> {
 public:
  /// Constructor
  StepperNewmarkImplicitDFormModifierXDefault() {}

  /// Destructor
  virtual ~StepperNewmarkImplicitDFormModifierXDefault() {}

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
                      const Scalar /* time */, const Scalar /* dt */,
                      const typename StepperNewmarkImplicitDFormModifierXBase<
                          Scalar>::MODIFIER_TYPE modType)
  {
    switch (modType) {
      case StepperNewmarkImplicitDFormModifierXBase<Scalar>::X_BEGIN_STEP:
      case StepperNewmarkImplicitDFormModifierXBase<Scalar>::X_BEFORE_SOLVE:
      case StepperNewmarkImplicitDFormModifierXBase<Scalar>::X_AFTER_SOLVE:
      case StepperNewmarkImplicitDFormModifierXBase<Scalar>::X_END_STEP: {
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

#endif  // Tempus_StepperNewmarkImplicitDFormModifierX_hpp
