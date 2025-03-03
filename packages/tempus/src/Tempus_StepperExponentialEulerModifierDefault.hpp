// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExponentialEulerModifierDefault_hpp
#define Tempus_StepperExponentialEulerModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExponentialEulerModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperExponentialEuler.hpp"


namespace Tempus {

/** \brief Default modifier for StepperExponentialEuler.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperExponentialEulerModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template<class Scalar>
class StepperExponentialEulerModifierDefault
  : virtual public Tempus::StepperExponentialEulerModifierBase<Scalar>
{
public:

  /// Constructor
  StepperExponentialEulerModifierDefault(){}

  /// Destructor
  virtual ~StepperExponentialEulerModifierDefault(){}

  /// Modify ExponentialEuler Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperExponentialEuler<Scalar> > /* stepper */,
    const typename StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperExponentialEulerAppAction<Scalar>::BEGIN_STEP:
      case StepperExponentialEulerAppAction<Scalar>::BEFORE_EXP:
      case StepperExponentialEulerAppAction<Scalar>::AFTER_EXP:
      case StepperExponentialEulerAppAction<Scalar>::END_STEP:
      {
        // No-op.
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }
  }

};

} // namespace Tempus

#endif // Tempus_StepperExponentialEulerModifierDefault_hpp
