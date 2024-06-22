//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperOperatorSplitAppAction_hpp
#define Tempus_StepperOperatorSplitAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperOperatorSplit;

/** \brief StepperOperatorSplitAppAction class for StepperOperatorSplit.
 *
 *  This class provides a means to apply various actions with the OperatorSplit
 * time step. The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperOperatorSplitAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperOperatorSplit.
 *
 * <b>Design Considerations</b>
 *   - StepperOperatorSplitAppAction is not stateless!  Developers may touch the
 *     solution state!  Developers need to be careful not to break the
 *     restart (checkpoint) capability.
 */
template <class Scalar>
class StepperOperatorSplitAppAction {
 public:
  enum ACTION_LOCATION {
    BEGIN_STEP,      ///< At the beginning of the step.
    BEFORE_STEPPER,  ///< Before a stepper evaluation.
    AFTER_STEPPER,   ///< After a stepper evaluation.
    END_STEP         ///< At the end of the step.
  };

  /// Constructor
  StepperOperatorSplitAppAction() {}

  /// Destructor
  virtual ~StepperOperatorSplitAppAction() {}

  /// Execute application action for OperatorSplit Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperOperatorSplit<Scalar> > stepper,
      const typename StepperOperatorSplitAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};
}  // namespace Tempus
#endif  // Tempus_StepperOperatorSplitAppAction_hpp
