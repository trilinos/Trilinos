//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperSubcyclingAppAction_hpp
#define Tempus_StepperSubcyclingAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperSubcycling;

/** \brief Application Action for StepperSubcycling.
 *
 *  This class provides a means to apply various actions with the Subcycling
 *  time step.  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperSubcyclingAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperSubcycling.
 */
template <class Scalar>
class StepperSubcyclingAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,  ///< At the beginning of the step.
    END_STEP     ///< At the end of the step.
  };

  /// Constructor
  StepperSubcyclingAppAction() {}

  /// Destructor
  virtual ~StepperSubcyclingAppAction() {}

  /// Execute application action for Subcycling Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperSubcycling<Scalar> > stepper,
      const typename StepperSubcyclingAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperSubcyclingAppAction_hpp
