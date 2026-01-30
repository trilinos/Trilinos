//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPIAppAction_hpp
#define Tempus_StepperEPIAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperEPI;

/** \brief Application Action for StepperEPI.
 *
 *  This class provides a means to apply various actions with the EPI
 * time step. The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperEPIAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperEPI.
 */
template <class Scalar>
class StepperEPIAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,            ///< At the beginning of the step.
    BEFORE_EXPLICIT_EVAL,  ///< Before the explicit evaluation.
    END_STEP               ///< At the end of the step.
  };

  /// Constructor
  StepperEPIAppAction() {}

  /// Destructor
  virtual ~StepperEPIAppAction() {}

  /// Execute application action for EPI Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperEPI<Scalar> > stepper,
      const typename StepperEPIAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperEPIAppAction_hpp
