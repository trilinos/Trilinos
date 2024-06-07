//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitDFormModifierDefault_hpp
#define Tempus_StepperNewmarkImplicitDFormModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperNewmarkImplicitDFormModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperNewmarkImplicitDForm.hpp"

namespace Tempus {

/** \brief Default modifier for StepperNewmarkImplicitDForm.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperNewmarkImplicitDFormModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperNewmarkImplicitDFormModifierDefault
  : virtual public Tempus::StepperNewmarkImplicitDFormModifierBase<Scalar> {
 public:
  /// Constructor
  StepperNewmarkImplicitDFormModifierDefault() {}

  /// Destructor
  virtual ~StepperNewmarkImplicitDFormModifierDefault() {}

  /// Modify NewmarkImplicitDForm Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperNewmarkImplicitDForm<Scalar> > /* stepper */,
      const typename StepperNewmarkImplicitDFormAppAction<
          Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperNewmarkImplicitDFormAppAction<Scalar>::BEGIN_STEP:
      case StepperNewmarkImplicitDFormAppAction<Scalar>::BEFORE_SOLVE:
      case StepperNewmarkImplicitDFormAppAction<Scalar>::AFTER_SOLVE:
      case StepperNewmarkImplicitDFormAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperNewmarkImplicitDFormModifierDefault_hpp
