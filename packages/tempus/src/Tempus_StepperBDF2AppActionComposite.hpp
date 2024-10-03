//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperBDF2AppActionComposite_hpp
#define Tempus_StepperBDF2AppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperBDF2AppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperBDF2AppActionComposite
  : virtual public Tempus::StepperBDF2AppAction<Scalar> {
 public:
  /// Default constructor
  StepperBDF2AppActionComposite();

  /// Destructor
  virtual ~StepperBDF2AppActionComposite();

  /// Execute application action for BDF2 Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperBDF2<Scalar> > stepper,
      const typename StepperBDF2AppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addBDF2AppAction(Teuchos::RCP<StepperBDF2AppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearBDF2AppActions();
  {
    appActions_.clear();
  }

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

 private:
  std::vector<Teuchos::RCP<StepperBDF2AppAction<Scalar> > > appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperBDF2AppActionComposite_hpp
