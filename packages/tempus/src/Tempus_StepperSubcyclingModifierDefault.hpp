//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperSubcyclingModifierDefault_hpp
#define Tempus_StepperSubcyclingModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperSubcyclingModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperSubcycling.hpp"

namespace Tempus {

/** \brief Default modifier for StepperSubcycling.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperSubcyclingModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template <class Scalar>
class StepperSubcyclingModifierDefault
  : virtual public Tempus::StepperSubcyclingModifierBase<Scalar> {
 public:
  /// Constructor
  StepperSubcyclingModifierDefault() {}

  /// Destructor
  virtual ~StepperSubcyclingModifierDefault() {}

  /// Modify Subcycling Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperSubcycling<Scalar> > /* stepper */,
      const typename StepperSubcyclingAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperSubcyclingAppAction<Scalar>::BEGIN_STEP:
      case StepperSubcyclingAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperSubcyclingModifierDefault_hpp
