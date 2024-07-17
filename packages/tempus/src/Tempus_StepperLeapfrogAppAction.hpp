//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperLeapfrogAppAction_hpp
#define Tempus_StepperLeapfrogAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperLeapfrog;

/** \brief Application Action for StepperLeapfrog.
 *
 *  This class provides a means to apply various actions with the Leapfrog time
 * step. The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperLeapfrogAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperLeapfrog.
 */
template <class Scalar>
class StepperLeapfrogAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,            ///< At the beginning of the step.
    BEFORE_X_UPDATE,       ///< Before updating x
    BEFORE_EXPLICIT_EVAL,  ///< Before the explicit ME evaluation.
    BEFORE_XDOT_UPDATE,    ///< Before updating xDot
    END_STEP               ///< At the end of the step.
  };

  /// Constructor
  StepperLeapfrogAppAction() {}

  /// Destructor
  virtual ~StepperLeapfrogAppAction() {}

  /// Execute application action for Leapfrog Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperLeapfrog<Scalar> > stepper,
      const typename StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperLeapfrogAppAction_hpp
