//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperSubcyclingObserverBase_hpp
#define Tempus_StepperSubcyclingObserverBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperSubcyclingAppAction.hpp"

namespace Tempus {

/** \brief Base observer for StepperSubcycling.
 *
 *  This class provides a means to observe values (e.g., solution variables
 *  through SolutionHistory, and stepper member data through the Stepper),
 *  and cannot modify them.
 *
 *  Users deriving from this class can observer a lot of data, and it is
 *  expected that users will NOT modify any of that data.  If the user
 *  wishes to modify the solution and/or stepper data during the
 *  Stepper::takeStep, they should use the Modifier class (with care!).
 *
 *  The locations for these AppAction calls
 *  (StepperSubcyclingAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperSubcycling.
 */
template <class Scalar>
class StepperSubcyclingObserverBase
  : virtual public Tempus::StepperSubcyclingAppAction<Scalar> {
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
      Teuchos::RCP<StepperSubcycling<Scalar> > stepper,
      const typename StepperSubcyclingAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    this->observe(sh, stepper, actLoc);
  }

 public:
  /// Observe Subcycling Stepper.
  virtual void observe(
      Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
      Teuchos::RCP<const StepperSubcycling<Scalar> > /* stepper */,
      const typename StepperSubcyclingAppAction<Scalar>::ACTION_LOCATION
          actLoc) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperSubcyclingObserverBase_hpp
