//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBackwardEulerAppAction_hpp
#define Tempus_StepperBackwardEulerAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperBackwardEuler;

/** \brief Application Action for StepperBackwardEuler.
 *
 *  This class provides a means to apply various actions with the BackwardEuler
 * time step. The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperBackwardEulerAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperBackwardEuler.
 */
template <class Scalar>
class StepperBackwardEulerAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,    ///< At the beginning of the step.
    BEFORE_SOLVE,  ///< Before the implicit solve.
    AFTER_SOLVE,   ///< After the implicit solve.
    END_STEP       ///< At the end of the step.
  };

  /// Constructor
  StepperBackwardEulerAppAction() {}

  /// Destructor
  virtual ~StepperBackwardEulerAppAction() {}

  /// Execute application action for BackwardEuler Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperBackwardEuler<Scalar> > stepper,
      const typename StepperBackwardEulerAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperBackwardEulerAppAction_hpp
