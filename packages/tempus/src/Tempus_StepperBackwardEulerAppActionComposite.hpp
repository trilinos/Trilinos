// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBackwardEulerAppActionComposite_hpp
#define Tempus_StepperBackwardEulerAppActionComposite_hpp

#include "Tempus_StepperBackwardEulerAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template<class Scalar>
class StepperBackwardEulerAppActionComposite
  : virtual public Tempus::StepperBackwardEulerAppAction<Scalar>
{
public:

  /// Default constructor
  StepperBackwardEulerAppActionComposite();

  /// Destructor
  virtual ~StepperBackwardEulerAppActionComposite();

  /// Execute application action for BackwardEuler Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperBackwardEuler<Scalar> > stepper,
    const typename StepperBackwardEulerAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    for(auto& a : appActions_)
      a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addBackwardEulerAppAction(Teuchos::RCP<StepperBackwardEulerAppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearBackwardEulerAppActions();
  { appActions_.clear();}

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

private:

  std::vector<Teuchos::RCP<StepperBackwardEulerAppAction<Scalar > > > appActions_;

};

} // namespace Tempus
#endif // Tempus_StepperBackwardEulerAppActionComposite_hpp
