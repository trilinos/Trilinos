//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitDFormAppAction_hpp
#define Tempus_StepperNewmarkImplicitDFormAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperNewmarkImplicitDForm;

/** \brief Application Action for StepperNewmarkImplicitDForm.
 *
 *  This class provides a means to apply various actions with the
 * NewmarkImplicitDForm time step. The data available to this class is solution
 * variables (through SolutionHistory), and stepper data (through the Stepper).
 * It allows the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperNewmarkImplicitDFormAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperNewmarkImplicitDForm.
 */
template <class Scalar>
class StepperNewmarkImplicitDFormAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,    ///< At the beginning of the step.
    BEFORE_SOLVE,  ///< Before the implicit solve.
    AFTER_SOLVE,   ///< After the implicit solve.
    END_STEP       ///< At the end of the step.
  };

  /// Constructor
  StepperNewmarkImplicitDFormAppAction() {}

  /// Destructor
  virtual ~StepperNewmarkImplicitDFormAppAction() {}

  /// Execute application action for NewmarkImplicitDForm Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperNewmarkImplicitDForm<Scalar> > stepper,
      const typename StepperNewmarkImplicitDFormAppAction<
          Scalar>::ACTION_LOCATION actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperNewmarkImplicitDFormAppAction_hpp
