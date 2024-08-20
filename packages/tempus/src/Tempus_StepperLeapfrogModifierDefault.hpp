//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperLeapfrogModifierDefault_hpp
#define Tempus_StepperLeapfrogModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperLeapfrogModifierBase.hpp"

namespace Tempus {

/** \brief Default modifier for StepperLeapfrog.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperLeapfrogModifierBase for details on the algorithm.
 */
template <class Scalar>
class StepperLeapfrogModifierDefault
  : virtual public Tempus::StepperLeapfrogModifierBase<Scalar> {
 public:
  /// Constructor
  StepperLeapfrogModifierDefault() {}

  /// Destructor
  virtual ~StepperLeapfrogModifierDefault() {}

  /// Modify Leapfrog Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperLeapfrog<Scalar> > /* stepper */,
      const typename StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperLeapfrogAppAction<Scalar>::BEGIN_STEP:
      case StepperLeapfrogAppAction<Scalar>::BEFORE_X_UPDATE:
      case StepperLeapfrogAppAction<Scalar>::BEFORE_EXPLICIT_EVAL:
      case StepperLeapfrogAppAction<Scalar>::BEFORE_XDOT_UPDATE:
      case StepperLeapfrogAppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperLeapfrogModifierDefault_hpp
