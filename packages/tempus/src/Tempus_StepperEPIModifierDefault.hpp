// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperEPIModifierDefault_hpp
#define Tempus_StepperEPIModifierDefault_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperEPIModifierBase.hpp"

// Applications can uncomment this include in their implementation,
// if they need access to the stepper methods.
//#include "Tempus_StepperEPI.hpp"


namespace Tempus {

/** \brief Default modifier for StepperEPI.
 *
 *  The default modifier provides no-op functionality for the modifier.
 *  See StepperEPIModifierBase for details on the algorithm.
 *
 *  Applications can copy this implementation, rename, implement their
 *  action, and set on the stepper to get app-specific functionality.
 */
template<class Scalar>
class StepperEPIModifierDefault
  : virtual public Tempus::StepperEPIModifierBase<Scalar>
{
public:

  /// Constructor
  StepperEPIModifierDefault(){}

  /// Destructor
  virtual ~StepperEPIModifierDefault(){}

  /// Modify EPI Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperEPI<Scalar> > /* stepper */,
    const typename StepperEPIAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    switch(actLoc) {
      case StepperEPIAppAction<Scalar>::BEGIN_STEP:
      case StepperEPIAppAction<Scalar>::BEFORE_EXP:
      case StepperEPIAppAction<Scalar>::AFTER_EXP:
      case StepperEPIAppAction<Scalar>::END_STEP:
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

#endif // Tempus_StepperEPIModifierDefault_hpp
