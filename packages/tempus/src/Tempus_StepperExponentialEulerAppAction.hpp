// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExponentialEulerAppAction_hpp
#define Tempus_StepperExponentialEulerAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration
template<class Scalar> class StepperExponentialEuler;

/** \brief Application Action for StepperExponentialEuler.
 *
 *  This class provides a means to apply various actions with the ExponentialEuler time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperExponentialEulerAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperExponentialEuler.
 */
template<class Scalar>
class StepperExponentialEulerAppAction
{
public:

  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,     ///< At the beginning of the step.
    BEFORE_EXP,     ///< Before the exponential solve.
    AFTER_EXP,      ///< After the exponential solve.
    END_STEP        ///< At the end of the step.
  };

  /// Constructor
  StepperExponentialEulerAppAction() {}

  /// Destructor
  virtual ~StepperExponentialEulerAppAction() {}

  /// Execute application action for ExponentialEuler Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperExponentialEuler<Scalar> > stepper,
    const typename StepperExponentialEulerAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperExponentialEulerAppAction_hpp
