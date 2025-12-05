// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperEPI3ModifierBase_hpp
#define Tempus_StepperEPI3ModifierBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperEPI3AppAction.hpp"


namespace Tempus {

/** \brief Base modifier for StepperEPI3.
 *
 *  This class provides a means to modify values (e.g., solution variables
 *  through SolutionHistory, and stepper member data through the Stepper),
 *  and can be very powerful and easy to make changes to the stepper and
 *  the solution.
 *
 *  Users deriving from this class can access a lot of data, and it is
 *  expected that those users know what changes are allowable without
 *  affecting the Stepper correctness, performance, accuracy and stability.
 *  Thus the user should be careful when accessing data through classes
 *  derived from the default modifier (i.e., USER BEWARE!!).
 *
 *  The locations for these AppAction calls
 *  (StepperEPI3AppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperEPI3.
 */
template<class Scalar>
class StepperEPI3ModifierBase
  : virtual public Tempus::StepperEPI3AppAction<Scalar>
{
private:

  /* \brief Adaptor execute function
   *
   *  This is an adaptor function to bridge between the AppAction
   *  interface and the Modifier interface.  It is meant to be private
   *  and non-virtual as deriving from this class should only need to
   *  implement the modify function.
   *
   *  For the Modifier interface, this adaptor is a "simple pass through".
   */
  void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperEPI3<Scalar> > stepper,
    const typename StepperEPI3AppAction<Scalar>::ACTION_LOCATION actLoc)
  { this->modify(sh, stepper, actLoc); }

public:

  /// Modify EPI3 Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperEPI3<Scalar> > /* stepper */,
    const typename StepperEPI3AppAction<Scalar>::ACTION_LOCATION actLoc) = 0;

};

} // namespace Tempus

#endif // Tempus_StepperEPI3ModifierBase_hpp
