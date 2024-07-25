//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBDF2AppAction_hpp
#define Tempus_StepperBDF2AppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperBDF2;

/** \brief Application Action for StepperBDF2.
 *
 *  This class provides a means to apply various actions with the BDF2 time
 * step. The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperBDF2AppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperBDF2.
 */
template <class Scalar>
class StepperBDF2AppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,    ///< At the beginning of the step.
    BEFORE_SOLVE,  ///< Before the solve.
    AFTER_SOLVE,   ///<  After the solve.
    END_STEP,      ///< At the end of the step.
  };

  /// Constructor
  StepperBDF2AppAction() {}

  /// Destructor
  virtual ~StepperBDF2AppAction() {}

  /// Execute application action for BDF2 Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperBDF2<Scalar> > stepper,
      const typename StepperBDF2AppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperBDF2AppAction_hpp
