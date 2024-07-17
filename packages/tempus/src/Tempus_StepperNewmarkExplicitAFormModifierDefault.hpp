//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkExplicitAFormModifierDefault_hpp
#define Tempus_StepperNewmarkExplicitAFormModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperNewmarkExplicitAFormModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperNewmarkExplicitAForm.hpp"

namespace Tempus {

/** \brief Default modifier for StepperNewmarkExplicitAForm.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperNewmarkExplicitAFormModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperNewmarkExplicitAFormModifierDefault
  : virtual public Tempus::StepperNewmarkExplicitAFormModifierBase<Scalar> {
 public:
  /// Constructor
  StepperNewmarkExplicitAFormModifierDefault() {}

  /// Destructor
  virtual ~StepperNewmarkExplicitAFormModifierDefault() {}

  /// Modify NewmarkExplicitAForm Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperNewmarkExplicitAForm<Scalar> > /* stepper */,
      const typename StepperNewmarkExplicitAFormAppAction<
          Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperNewmarkExplicitAFormAppAction<Scalar>::BEGIN_STEP:
      case StepperNewmarkExplicitAFormAppAction<Scalar>::BEFORE_EXPLICIT_EVAL:
      case StepperNewmarkExplicitAFormAppAction<Scalar>::AFTER_EXPLICIT_EVAL:
      case StepperNewmarkExplicitAFormAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperNewmarkExplicitAFormModifierDefault_hpp
