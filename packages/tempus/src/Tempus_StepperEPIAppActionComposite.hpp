//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPIAppActionComposite_hpp
#define Tempus_StepperEPIAppActionComposite_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperEPIAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class StepperEPIAppActionComposite
  : virtual public Tempus::StepperEPIAppAction<Scalar> {
 public:
  /// Default constructor
  StepperEPIAppActionComposite();

  /// Destructor
  virtual ~StepperEPIAppActionComposite();

  /// Execute application action for EPI Stepper.
  virtual void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperEPI<Scalar> > stepper,
      const typename StepperEPIAppAction<Scalar>::ACTION_LOCATION
          actLoc)
  {
    for (auto& a : appActions_) a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addEPIAppAction(
      Teuchos::RCP<StepperEPIAppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearEPIAppActions();
  {
    appActions_.clear();
  }

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

 private:
  std::vector<Teuchos::RCP<StepperEPIAppAction<Scalar> > > appActions_;
};

}  // namespace Tempus
#endif  // Tempus_StepperEPIAppActionComposite_hpp
