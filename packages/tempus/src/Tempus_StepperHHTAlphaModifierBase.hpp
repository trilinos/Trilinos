//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperHHTAlphaModifierBase_hpp
#define Tempus_StepperHHTAlphaModifierBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperHHTAlphaAppAction.hpp"

namespace Tempus {

/** \brief Base modifier for StepperHHTAlpha.
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
 */
template <class Scalar>
class StepperHHTAlphaModifierBase
  : virtual public Tempus::StepperHHTAlphaAppAction<Scalar> {
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
      Teuchos::RCP<StepperHHTAlpha<Scalar> > stepper,
      const typename StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    this->modify(sh, stepper, actLoc);
  }

 public:
  /// Modify HHT Alpha Stepper.
  virtual void modify(
      Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<StepperHHTAlpha<Scalar> > /* stepper */,
      const typename StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperHHTAlphaModifierBase_hpp
