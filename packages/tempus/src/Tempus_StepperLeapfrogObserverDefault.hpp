//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperLeapfrogObserverDefault_hpp
#define Tempus_StepperLeapfrogObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperLeapfrogObserverBase.hpp"

namespace Tempus {

/** \brief Default observer for StepperLeapfrog.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperLeapfrogObserverBase for details on the algorithm.
 */
template <class Scalar>
class StepperLeapfrogObserverDefault
  : virtual public Tempus::StepperLeapfrogObserverBase<Scalar> {
 public:
  /// Constructor
  StepperLeapfrogObserverDefault() {}

  /// Destructor
  virtual ~StepperLeapfrogObserverDefault() {}

  /// Observe Leapfrog Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperLeapfrog<Scalar> > /* stepper */,
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

#endif  // Tempus_StepperLeapfrogObserverDefault_hpp
