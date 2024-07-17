//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperOperatorSplitAppActionComposite_hpp
#define Tempus_StepperOperatorSplitAppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperOperatorSplitAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperOperatorSplitAppActionComposite
  : virtual public Tempus::StepperOperatorSplitAppAction<Scalar> {
 public:
  /// Default constructor
  StepperOperatorSplitAppActionComposite();

  /// Destructor
  virtual ~StepperOperatorSplitAppActionComposite();

  /// Execute application action for OperatorSplit Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperOperatorSplit<Scalar> > stepper,
      const typename StepperOperatorSplitAppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addOperatorSplitAppAction(
      Teuchos::RCP<StepperOperatorSplitAppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearOperatorSplitAppActions();
  {
    appActions_.clear();
  }

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

 private:
  std::vector<Teuchos::RCP<StepperOperatorSplitAppAction<Scalar> > >
      appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperOperatorSplitAppActionComposite_hpp
