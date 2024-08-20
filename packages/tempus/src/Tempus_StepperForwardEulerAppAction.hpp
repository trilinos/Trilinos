//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperForwardEulerAppAction_hpp
#define Tempus_StepperForwardEulerAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperForwardEuler;

/** \brief Application Action for StepperForwardEuler.
 *
 *  This class provides a means to apply various actions with the ForwardEuler
 * time step. The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperForwardEulerAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperForwardEuler.
 */
template <class Scalar>
class StepperForwardEulerAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,            ///< At the beginning of the step.
    BEFORE_EXPLICIT_EVAL,  ///< Before the explicit evaluation.
    END_STEP               ///< At the end of the step.
  };

  /// Constructor
  StepperForwardEulerAppAction() {}

  /// Destructor
  virtual ~StepperForwardEulerAppAction() {}

  /// Execute application action for ForwardEuler Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperForwardEuler<Scalar> > stepper,
      const typename StepperForwardEulerAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperForwardEulerAppAction_hpp
