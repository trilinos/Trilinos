// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperLeapfrogAppActionComposite_hpp
#define Tempus_StepperLeapfrogAppActionComposite_hpp

#include "Tempus_StepperLeapfrogAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template<class Scalar>
class StepperLeapfrogAppActionComposite
  : virtual public Tempus::StepperLeapfrogAppAction<Scalar>
{
public:

  /// Default constructor
  StepperLeapfrogAppActionComposite();

  /// Destructor
  virtual ~StepperLeapfrogAppActionComposite();

  /// Execute application action for Leapfrog Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperLeapfrog<Scalar> > stepper,
    const typename StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    for(auto& a : appActions_)
      a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addLeapfrogAppAction(Teuchos::RCP<StepperLeapfrogAppAction<Scalar> > appAction)
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearLeapfrogAppActions()
  { appActions_.clear();}

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

private:

  std::vector<Teuchos::RCP<StepperLeapfrogAppAction<Scalar > > > appActions_;

};

} // namespace Tempus
#endif // Tempus_StepperLeapfrogAppActionComposite_hpp
