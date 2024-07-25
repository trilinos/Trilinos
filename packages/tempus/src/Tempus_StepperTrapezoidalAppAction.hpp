//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperTrapezoidalAppAction_hpp
#define Tempus_StepperTrapezoidalAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperTrapezoidal;

/** \brief Application Action for StepperTrapezoidal.
 *
 *  This class provides a means to apply various actions with the Trapezoidal
 * time step. The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperTrapezoidalAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperTrapezoidal.
 */
template <class Scalar>
class StepperTrapezoidalAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,    ///< At the beginning of the step.
    BEFORE_SOLVE,  ///< Before the solve
    AFTER_SOLVE,   ///< After then solve
    END_STEP       ///< At the end of the step.
  };

  /// Constructor
  StepperTrapezoidalAppAction() {}

  /// Destructor
  virtual ~StepperTrapezoidalAppAction() {}

  /// Execute application action for Trapezoidal Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperTrapezoidal<Scalar> > stepper,
      const typename StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperTrapezoidalAppAction_hpp
