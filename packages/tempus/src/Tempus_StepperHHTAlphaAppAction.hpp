// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperHHTAlphaAppAction_hpp
#define Tempus_StepperHHTAlphaAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this AppAction <--> Stepper)
template<class Scalar> class StepperHHTAlpha;

/** \brief Application Action for HHT Alpha.
 *
 *  This class provides a means to apply various actions with the HHT Alpha time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data (i.e., use but not change the
 *  data) to change any of it (USER BEWARE!).
 *
 *  Below is the HHT Alpha algorithm and includes the locations where the
 *  application can take actions (in italicized).
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{HHT Alpha with application-action locations indicated.}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State Compute the predictor (e.g., apply stepper to $x_n$).
 *    \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}
 *    \State Solve for $x_{n+1}$ using the HHT one-step update.                                                                                        *    \State {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *    \State Update $\dot x_{n+1}$.           
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperHHTAlphaAppAction
{
public:

  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,     ///< At the beginning of the step.
    BEFORE_SOLVE,   ///< Before the implicit solve.
    AFTER_SOLVE,    ///< After the implicit solve.
    END_STEP        ///< At the end of the step.
  };

  /// Constructor
  StepperHHTAlphaAppAction(){}

  /// Destructor
  virtual ~StepperHHTAlphaAppAction(){}

  /// Execute application action for HHTAlpha Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperHHTAlpha<Scalar> > stepper,
    const typename StepperHHTAlphaAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperHHTAlphaAppAction_hpp
