//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkExplicitAFormAppActionComposite_hpp
#define Tempus_StepperNewmarkExplicitAFormAppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperNewmarkExplicitAFormAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperNewmarkExplicitAFormAppActionComposite
  : virtual public Tempus::StepperNewmarkExplicitAFormAppAction<Scalar> {
 public:
  /// Default constructor
  StepperNewmarkExplicitAFormAppActionComposite();

  /// Destructor
  virtual ~StepperNewmarkExplicitAFormAppActionComposite();

  /// Execute application action for NewmarkExplicitAForm Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperNewmarkExplicitAForm<Scalar> > stepper,
      const typename StepperNewmarkExplicitAFormAppAction<
          Scalar>::ACTION_LOCATION actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addNewmarkExplicitAFormAppAction(
      Teuchos::RCP<StepperNewmarkExplicitAFormAppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearAppActions();
  {
    appActions_.clear();
  }

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

 private:
  std::vector<Teuchos::RCP<StepperNewmarkExplicitAFormAppAction<Scalar> > >
      appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperNewmarkExplicitAFormAppActionComposite_hpp
