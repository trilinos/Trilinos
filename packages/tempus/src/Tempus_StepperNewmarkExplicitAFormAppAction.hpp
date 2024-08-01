//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkExplicitAFormAppAction_hpp
#define Tempus_StepperNewmarkExplicitAFormAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperNewmarkExplicitAForm;

/** \brief Application Action for StepperNewmarkExplicitAForm.
 *
 *  This class provides a means to apply various actions with the
 * NewmarkExplicitAForm time step. The data available to this class is solution
 * variables (through SolutionHistory), and stepper data (through the Stepper).
 * It allows the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperNewmarkExplicitAFormAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperNewmarkExplicitAForm.
 */
template <class Scalar>
class StepperNewmarkExplicitAFormAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,            ///< At the beginning of the step.
    BEFORE_EXPLICIT_EVAL,  ///< Before the explicit evaluation.
    AFTER_EXPLICIT_EVAL,   ///< After the explicit evaluation.
    END_STEP               ///< At the end of the step.
  };

  /// Constructor
  StepperNewmarkExplicitAFormAppAction() {}

  /// Destructor
  virtual ~StepperNewmarkExplicitAFormAppAction() {}

  /// Execute application action for NewmarkExplicitAForm Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperNewmarkExplicitAForm<Scalar> > stepper,
      const typename StepperNewmarkExplicitAFormAppAction<
          Scalar>::ACTION_LOCATION actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperNewmarkExplicitAFormAppAction_hpp
