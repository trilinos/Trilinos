//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitDFormAppActionComposite_hpp
#define Tempus_StepperNewmarkImplicitDFormAppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperNewmarkImplicitDFormAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperNewmarkImplicitDFormAppActionComposite
  : virtual public Tempus::StepperNewmarkImplicitDFormAppAction<Scalar> {
 public:
  /// Default constructor
  StepperNewmarkImplicitDFormAppActionComposite();

  /// Destructor
  virtual ~StepperNewmarkImplicitDFormAppActionComposite();

  /// Execute application action for NewmarkImplicitDForm Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperNewmarkImplicitDForm<Scalar> > stepper,
      const typename StepperNewmarkImplicitDFormAppAction<
          Scalar>::ACTION_LOCATION actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addNewmarkImplicitDFormAppAction(
      Teuchos::RCP<StepperNewmarkImplicitDFormAppAction<Scalar> > appAction);
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
  std::vector<Teuchos::RCP<StepperNewmarkImplicitDFormAppAction<Scalar> > >
      appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperNewmarkImplicitDFormAppActionComposite_hpp
