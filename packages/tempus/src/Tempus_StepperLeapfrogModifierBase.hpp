//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperLeapfrogModifierBase_hpp
#define Tempus_StepperLeapfrogModifierBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperLeapfrogAppAction.hpp"

namespace Tempus {

/** \brief Base modifier for StepperLeapfrog.
 *
 *  This class provides a means to modify values (e.g., solution variables
 *  through SolutionHistory, and stepper member data through the Stepper),
 *  and can be very powerful and easy to make changes to the stepper and
 *  the solution.
 *
 *  Users deriving from this class can access a lot of data, and it is
 *  expected that those users know what changes are allowable without
 *  affecting the Stepper correctness, performance, accuracy and stability.
 *  Thus the user should be careful when accessing data through classes
 *  derived from the default modifier (i.e., USER BEWARE!!).
 *
 *  The locations for these AppAction calls
 *  (StepperLeapfrogAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperLeapfrog.
 */

template <class Scalar>
class StepperLeapfrogModifierBase
  : virtual public Tempus::StepperLeapfrogAppAction<Scalar> {
 private:
  /* \brief Adaptor execute function
   *
   *  This is an adaptor function to bridge between the AppAction
   *  interface and the Modifier interface.  It is meant to be private
   *  and non-virtual as deriving from this class should only need to
   *  implement the modify function.
   *
   *  For the Modifier interface, this adaptor is a "simple pass through".
   */
  void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperLeapfrog<Scalar> > stepper,
      const typename StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    this->modify(sh, stepper, actLoc);
  }

 public:
  /// Modify Leapfrog Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperLeapfrog<Scalar> > /* stepper */,
      const typename StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperLeapfrogModifierBase_hpp
