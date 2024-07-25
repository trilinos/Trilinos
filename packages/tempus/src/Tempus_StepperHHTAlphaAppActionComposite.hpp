//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperHHTAlphaAppActionComposite_hpp
#define Tempus_StepperHHTAlphaAppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperHHTAlphaAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperHHTAlphaAppActionComposite
  : virtual public Tempus::StepperHHTAlphaAppAction<Scalar> {
 public:
  /// Default constructor
  StepperHHTAlphaAppActionComposite();

  /// Destructor
  virtual ~StepperHHTAlphaAppActionComposite();

  /// Execute application action for HHTAlpha Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperHHTAlpha<Scalar> > stepper,
      const typename StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addHHTAlphaAppAction(
      Teuchos::RCP<StepperHHTAlphaAppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearHHTAlphaAppActions();
  {
    appActions_.clear();
  }

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

 private:
  std::vector<Teuchos::RCP<StepperHHTAlphaAppAction<Scalar> > > appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperHHTAlphaAppActionComposite_hpp
