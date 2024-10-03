//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBDF2ModifierDefault_hpp
#define Tempus_StepperBDF2ModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBDF2ModifierBase.hpp"

namespace Tempus {

/** \brief Default modifier for StepperBDF2.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperBDF2ModifierBase for details on the algorithm.
 */
template <class Scalar>
class StepperBDF2ModifierDefault
  : virtual public Tempus::StepperBDF2ModifierBase<Scalar> {
 public:
  /// Constructor
  StepperBDF2ModifierDefault() {}

  /// Destructor
  virtual ~StepperBDF2ModifierDefault() {}

  /// Modify BDF2 Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperBDF2<Scalar> > /* stepper */,
      const typename StepperBDF2AppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch (actLoc) {
      case StepperBDF2AppAction<Scalar>::BEGIN_STEP:
      case StepperBDF2AppAction<Scalar>::BEFORE_SOLVE:
      case StepperBDF2AppAction<Scalar>::AFTER_SOLVE:
      case StepperBDF2AppAction<Scalar>::END_STEP: {
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

#endif  // Tempus_StepperBDF2ModifierDefault_hpp
