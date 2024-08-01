//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitAFormAppAction_hpp
#define Tempus_StepperNewmarkImplicitAFormAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperNewmarkImplicitAForm;

/** \brief Application Action for StepperNewmarkImplicitAForm.
 *
 *  This class provides a means to apply various actions with the
 * NewmarkImplicitAForm time step. The data available to this class is solution
 * variables (through SolutionHistory), and stepper data (through the Stepper).
 * It allows the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperNewmarkImplicitAFormAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperNewmarkImplicitAForm.
 */
template <class Scalar>
class StepperNewmarkImplicitAFormAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,    ///< At the beginning of the step.
    BEFORE_SOLVE,  ///< Before the implicit solve.
    AFTER_SOLVE,   ///< After the implicit solve.
    END_STEP       ///< At the end of the step.
  };

  /// Constructor
  StepperNewmarkImplicitAFormAppAction() {}

  /// Destructor
  virtual ~StepperNewmarkImplicitAFormAppAction() {}

  /// Execute application action for NewmarkImplicitAForm Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperNewmarkImplicitAForm<Scalar> > stepper,
      const typename StepperNewmarkImplicitAFormAppAction<
          Scalar>::ACTION_LOCATION actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperNewmarkImplicitAFormAppAction_hpp
