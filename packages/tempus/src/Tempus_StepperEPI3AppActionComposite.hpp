//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPI3AppActionComposite_hpp
#define Tempus_StepperEPI3AppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperEPI3AppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperEPI3AppActionComposite
  : virtual public Tempus::StepperEPI3AppAction<Scalar> {
 public:
  /// Default constructor
  StepperEPI3AppActionComposite();

  /// Destructor
  virtual ~StepperEPI3AppActionComposite();

  /// Execute application action for EPI3 Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperEPI3<Scalar> > stepper,
      const typename StepperEPI3AppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addEPI3AppAction(
      Teuchos::RCP<StepperEPI3AppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearEPI3AppActions();
  {
    appActions_.clear();
  }

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

 private:
  std::vector<Teuchos::RCP<StepperEPI3AppAction<Scalar> > > appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperEPI3AppActionComposite_hpp
