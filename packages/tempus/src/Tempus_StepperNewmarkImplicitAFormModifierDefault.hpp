//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitAFormModifierDefault_hpp
#define Tempus_StepperNewmarkImplicitAFormModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperNewmarkImplicitAFormModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperNewmarkImplicitAForm.hpp"

namespace Tempus {

/** \brief Default modifier for StepperNewmarkImplicitAForm.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperNewmarkImplicitAFormModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperNewmarkImplicitAFormModifierDefault
  : virtual public Tempus::StepperNewmarkImplicitAFormModifierBase<Scalar> {
 public:
  /// Constructor
  StepperNewmarkImplicitAFormModifierDefault() {}

  /// Destructor
  virtual ~StepperNewmarkImplicitAFormModifierDefault() {}

  /// Modify NewmarkImplicitAForm Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperNewmarkImplicitAForm<Scalar> > /* stepper */,
      const typename StepperNewmarkImplicitAFormAppAction<
          Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperNewmarkImplicitAFormAppAction<Scalar>::BEGIN_STEP:
      case StepperNewmarkImplicitAFormAppAction<Scalar>::BEFORE_SOLVE:
      case StepperNewmarkImplicitAFormAppAction<Scalar>::AFTER_SOLVE:
      case StepperNewmarkImplicitAFormAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperNewmarkImplicitAFormModifierDefault_hpp
