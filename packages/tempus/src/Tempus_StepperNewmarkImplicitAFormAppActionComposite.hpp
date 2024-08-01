//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitAFormAppActionComposite_hpp
#define Tempus_StepperNewmarkImplicitAFormAppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperNewmarkImplicitAFormAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperNewmarkImplicitAFormAppActionComposite
  : virtual public Tempus::StepperNewmarkImplicitAFormAppAction<Scalar> {
 public:
  /// Default constructor
  StepperNewmarkImplicitAFormAppActionComposite();

  /// Destructor
  virtual ~StepperNewmarkImplicitAFormAppActionComposite();

  /// Execute application action for NewmarkImplicitAForm Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperNewmarkImplicitAForm<Scalar> > stepper,
      const typename StepperNewmarkImplicitAFormAppAction<
          Scalar>::ACTION_LOCATION actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addNewmarkImplicitAFormAppAction(
      Teuchos::RCP<StepperNewmarkImplicitAFormAppAction<Scalar> > appAction);
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
  std::vector<Teuchos::RCP<StepperNewmarkImplicitAFormAppAction<Scalar> > >
      appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperNewmarkImplicitAFormAppActionComposite_hpp
