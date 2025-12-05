// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperEPI3ModifierDefault_hpp
#define Tempus_StepperEPI3ModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperEPI3ModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperEPI3.hpp"


namespace Tempus {

/** \brief Default modifier for StepperEPI3.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperEPI3ModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template<class Scalar>
class StepperEPI3ModifierDefault
  : virtual public Tempus::StepperEPI3ModifierBase<Scalar>
{
public:

  /// Constructor
  StepperEPI3ModifierDefault(){}

  /// Destructor
  virtual ~StepperEPI3ModifierDefault(){}

  /// Modify EPI3 Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperEPI3<Scalar> > /* stepper */,
    const typename StepperEPI3AppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperEPI3AppAction<Scalar>::BEGIN_STEP:
      case StepperEPI3AppAction<Scalar>::BEFORE_EXP:
      case StepperEPI3AppAction<Scalar>::AFTER_EXP:
      case StepperEPI3AppAction<Scalar>::END_STEP:
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

#endif // Tempus_StepperEPI3ModifierDefault_hpp
