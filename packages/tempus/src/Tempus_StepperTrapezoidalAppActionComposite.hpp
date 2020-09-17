// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperTrapezoidalAppActionComposite_hpp
#define Tempus_StepperTrapezoidalAppActionComposite_hpp

#include "Tempus_StepperTrapezoidalAppAction.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This composite AppAction loops over added AppActions.
 *
 *  Inidividual AppActions are executed in the order in which they
 *  were added.
 */
template<class Scalar>
class StepperTrapezoidalAppActionComposite
  : virtual public Tempus::StepperTrapezoidalAppAction<Scalar>
{
public:

  /// Default constructor
  StepperTrapezoidalAppActionComposite();

  /// Destructor
  virtual ~StepperTrapezoidalAppActionComposite();

  /// Execute application action for Subcycling Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperTrapezoidal<Scalar> > stepper,
    const typename StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    for(auto& a : appActions_)
      a->execute(sh, stepper, actLoc);
  }

  // Add AppAction to the AppAction vector.
  void addTrapezoidalAppAction(Teuchos::RCP<StepperTrapezoidalAppAction<Scalar> > appAction);
  {
    appActions_.push_back(appAction);
  }

  // Clear the AppAction vector.
  void clearTrapezoidalAppActions();
  { appActions_.clear();}

  // Return the size of the AppAction vector.
  std::size_t getSize() const { return appActions_.size(); }

private:

  std::vector<Teuchos::RCP<StepperTrapezoidalAppAction<Scalar > > > appActions_;

};

} // namespace Tempus
#endif // Tempus_StepperTrapezoidalAppActionComposite_hpp
