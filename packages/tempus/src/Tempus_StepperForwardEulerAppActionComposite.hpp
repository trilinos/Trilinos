//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperForwardEulerAppActionComposite_hpp
#define Tempus_StepperForwardEulerAppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperForwardEulerAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperForwardEulerAppActionComposite
  : virtual public Tempus::StepperForwardEulerAppAction<Scalar> {
 public:
  /// Default constructor
  StepperForwardEulerAppActionComposite();

  /// Destructor
  virtual ~StepperForwardEulerAppActionComposite();

  /// Execute application action for ForwardEuler Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperForwardEuler<Scalar> > stepper,
      const typename StepperForwardEulerAppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addForwardEulerAppAction(
      Teuchos::RCP<StepperForwardEulerAppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearForwardEulerAppActions();
  {
    appActions_.clear();
  }

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

 private:
  std::vector<Teuchos::RCP<StepperForwardEulerAppAction<Scalar> > > appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperForwardEulerAppActionComposite_hpp
