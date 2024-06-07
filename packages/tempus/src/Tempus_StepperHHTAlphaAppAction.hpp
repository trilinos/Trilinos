//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperHHTAlphaAppAction_hpp
#define Tempus_StepperHHTAlphaAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

// Forward Declaration
template <class Scalar>
class StepperHHTAlpha;

/** \brief Application Action for HHT Alpha.
 *
 *  This class provides a means to apply various actions with the HHT Alpha time
 * step. The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 */
template <class Scalar>
class StepperHHTAlphaAppAction {
 public:
  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,    ///< At the beginning of the step.
    BEFORE_SOLVE,  ///< Before the implicit solve.
    AFTER_SOLVE,   ///< After the implicit solve.
    END_STEP       ///< At the end of the step.
  };

  /// Constructor
  StepperHHTAlphaAppAction() {}

  /// Destructor
  virtual ~StepperHHTAlphaAppAction() {}

  /// Execute application action for HHTAlpha Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperHHTAlpha<Scalar> > stepper,
      const typename StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperHHTAlphaAppAction_hpp
