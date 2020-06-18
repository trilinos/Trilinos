// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperOperatorSplitAppAction_hpp
#define Tempus_StepperOperatorSplitAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this Observer <--> Stepper)
template<class Scalar> class StepperOperatorSplit;

/** \brief StepperOperatorSplitAppAction class for StepperOperatorSplit.
 *
 * This is a means for application developers to perform tasks
 * during the time steps, e.g.,
 *   - Compute specific quantities
 *   - Output information
 *   - "Massage" the working solution state
 *   - ...
 *
 * <b>Design Considerations</b>
 *   - StepperOperatorSplitAppAction is not stateless!  Developers may touch the
 *     solution state!  Developers need to be careful not to break the
 *     restart (checkpoint) capability.
 */
template<class Scalar>
class StepperOperatorSplitAppAction
{
public:

  enum ACTION_LOCATION {
    BEGIN_STEP,     ///< At the beginning of the step.
    BEFORE_STEPPER, ///< Before a stepper evaluation.
    AFTER_STEPPER,  ///< After a stepper evaluation.
    END_STEP        ///< At the end of the step.
  };

  /// Constructor
  StepperOperatorSplitAppAction(){}

  /// Destructor
  virtual ~StepperOperatorSplitAppAction(){}

  /// Execute application action for OperatorSplit Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperOperatorSplit<Scalar> > stepper,
    const typename StepperOperatorSplitAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};
} // namespace Tempus
#endif // Tempus_StepperOperatorSplitAppAction_hpp




