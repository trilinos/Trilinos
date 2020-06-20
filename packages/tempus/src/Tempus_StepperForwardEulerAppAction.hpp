// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperForwardEulerAppAction_hpp
#define Tempus_StepperForwardEulerAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this AppAction <--> Stepper)
template<class Scalar> class StepperForwardEuler;

/** \brief Application Action for StepperForwardEuler.
 *
 *  This class provides a means to apply various actions with the ForwardEuler time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data (i.e., use but not change the
 *  data) to change any of it (USER BEWARE!).
 *
 *  Below is the ForwardEuler algorithm and includes the locations where the
 *  application can take actions (in italicized).
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{Forward Euler with the locations of the application actions indicated.}
 *  \begin{algorithmic}[1]
 *    \State Start with $x_n$, $\Delta t_n$
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State Form $f(x_{n},t_{n})$
 *    \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_EXPLICIT\_EVAL)}
 *    \State Form $x_n \leftarrow x_{n} + \Delta t_n f(x_{n},t_n)$
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperForwardEulerAppAction
{
public:

  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,     ///< At the beginning of the step.
    BEFORE_EXPLICIT_EVAL,   ///< Before the explicit evaluation.
    END_STEP        ///< At the end of the step.
  };

  /// Constructor
  StepperForwardEulerAppAction(){}

  /// Destructor
  virtual ~StepperForwardEulerAppAction(){}

  /// Execute application action for ForwardEuler Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperForwardEuler<Scalar> > stepper,
    const typename StepperForwardEulerAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperForwardEulerAppAction_hpp
