//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperRKAppAction_hpp
#define Tempus_StepperRKAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperRKBase;

/** \brief Application Action for StepperRKBase.
 *
 *  This class provides a means to apply various actions with the
 *  RK time step.  The data available to this class is solution
 *  variables (through SolutionHistory), and stepper data (through
 *  the Stepper).  It allows the application to just observe this
 *  data, i.e., use but not change any of it (USER BEWARE!).
 *
 *  The locations of the RK AppActions (StepperRKAppAction::ACTION_LOCATION)
 *  in takeStep are documented in each of the RK Algorithm sections:
 *  StepperExplicitRK, StepperDIRK and StepperIMEX_RK.
 */
template <class Scalar>
class StepperRKAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,            ///< At the beginning of the step.
    BEGIN_STAGE,           ///< At the beginning of the stage.
    BEFORE_SOLVE,          ///< Before the implicit solve.
    AFTER_SOLVE,           ///< After the implicit solve.
    BEFORE_EXPLICIT_EVAL,  ///< Before the explicit evaluation.
    END_STAGE,             ///< At the end of the stage.
    END_STEP               ///< At the end of the step.
  };

  /// Constructor
  StepperRKAppAction() {}

  /// Destructor
  virtual ~StepperRKAppAction() {}

  /// Execute application action for RK Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperRKBase<Scalar> > stepper,
      const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperRKAppAction_hpp
