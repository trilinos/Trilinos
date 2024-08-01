//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperSubcyclingAppActionComposite_hpp
#define Tempus_StepperSubcyclingAppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperSubcyclingAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperSubcyclingAppActionComposite
  : virtual public Tempus::StepperSubcyclingAppAction<Scalar> {
 public:
  /// Default constructor
  StepperSubcyclingAppActionComposite();

  /// Destructor
  virtual ~StepperSubcyclingAppActionComposite();

  /// Execute application action for Subcycling Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperSubcycling<Scalar> > stepper,
      const typename StepperSubcyclingAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addSubcyclingAppAction(
      Teuchos::RCP<StepperSubcyclingAppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearSubcyclingAppActions();
  {
    appActions_.clear();
  }

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

 private:
  std::vector<Teuchos::RCP<StepperSubcyclingAppAction<Scalar> > > appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperSubcyclingAppActionComposite_hpp
