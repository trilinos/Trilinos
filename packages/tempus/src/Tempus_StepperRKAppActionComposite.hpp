//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperRKAppActionComposite_hpp
#define Tempus_StepperRKAppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperRKAppAction.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Individual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperRKAppActionComposite
  : virtual public Tempus::StepperRKAppAction<Scalar> {
 public:
  /// Default constructor
  StepperRKAppActionComposite() {}

  /// Destructor
  virtual ~StepperRKAppActionComposite() {}

  /// Execute application action for RK Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperRKBase<Scalar> > stepper,
      const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addRKAppAction(Teuchos::RCP<StepperRKAppAction<Scalar> > appAction)
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearRKAppActions() { appActions_.clear(); }

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

 private:
  std::vector<Teuchos::RCP<StepperRKAppAction<Scalar> > > appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperRKAppActionComposite_hpp
