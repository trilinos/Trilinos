//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperHHTAlphaObserverBase_hpp
#define Tempus_StepperHHTAlphaObserverBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperHHTAlphaAppAction.hpp"

namespace Tempus {

/** \brief Base observer for StepperHHTAlpha.
 *
 *  This class provides a means to observe values (e.g., solution variables
 *  through SolutionHistory, and stepper member data through the Stepper),
 *  and cannot modify them.
 *
 *  Users deriving from this class can observer a lot of data, and it is
 *  expected that users will NOT modify any of that data.  If the user
 *  wishes to modify the solution and/or stepper data during the
 *  Stepper::takeStep, they should use the Modifier class (with care!).
 */
template <class Scalar>
class StepperHHTAlphaObserverBase
  : virtual public Tempus::StepperHHTAlphaAppAction<Scalar> {
 private:
  /* \brief Adaptor execute function
   *
   *  This is an adaptor function to bridge between the AppAction
   *  interface and this derived interface.  It is meant to be private
   *  and non-virtual as deriving from this class should only need to
   *  implement the observe function.
   *
   *  For the Observer interface, this adaptor simply "applies" constantness
   *  to the arguments.
   */
  void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperHHTAlpha<Scalar> > stepper,
      const typename StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    this->observe(sh, stepper, actLoc);
  }

 public:
  /// Observe HHTAlpha Stepper.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperHHTAlpha<Scalar> > /* stepper */,
      const typename StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperHHTAlphaObserverBase_hpp
