//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBDF2ObserverDefault_hpp
#define Tempus_StepperBDF2ObserverDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBDF2ObserverBase.hpp"

namespace Tempus {

/** \brief Default observer for StepperBDF2.
 *
 *  The default observer provides no-op functionality for the observer.
 *  See StepperBDF2ObserverBase for details on the algorithm.
 */
template <class Scalar>
class StepperBDF2ObserverDefault
  : virtual public Tempus::StepperBDF2ObserverBase<Scalar> {
 public:
  /// Constructor
  StepperBDF2ObserverDefault() {}

  /// Destructor
  virtual ~StepperBDF2ObserverDefault() {}

  /// Observe BDF2 Stepper at end of takeStep.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperBDF2<Scalar> > /* stepper */,
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

#endif  // Tempus_StepperBDF2ObserverDefault_hpp
