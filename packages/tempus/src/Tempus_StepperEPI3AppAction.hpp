// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperEPI3AppAction_hpp
#define Tempus_StepperEPI3AppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration
template<class Scalar> class StepperEPI3;

/** \brief Application Action for StepperEPI3.
 *
 *  This class provides a means to apply various actions with the EPI3 time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperEPI3AppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperEPI3.
 */
template<class Scalar>
class StepperEPI3AppAction
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
  StepperEPI3AppAction() {}

  /// Destructor
  virtual ~StepperEPI3AppAction() {}

  /// Execute application action for EPI3 Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperEPI3<Scalar> > stepper,
    const typename StepperEPI3AppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperEPI3AppAction_hpp
