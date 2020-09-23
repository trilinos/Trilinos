// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkExplicitAFormAppAction_hpp
#define Tempus_StepperNewmarkExplicitAFormAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this AppAction <--> Stepper)
template<class Scalar> class StepperNewmarkExplicitAForm;

/** \brief Application Action for StepperNewmarkExplicitAForm.
 *
 *  This class provides a means to apply various actions with the NewmarkExplicitAForm time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data (i.e., use but not change the
 *  data) to change any of it (USER BEWARE!).
 *
 *  Below is the NewmarkExplicitAForm algorithm and includes the locations where the
 *  application can take actions (in italicized).
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{Newmark Explicit-A with application-action locations indicated.}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State Compute displacement and velocity predictors 
 *    \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_EXPLICIT\_EVAL)}
 *    \State Form $\ddot{x}_n \leftarrow x_{n} + \Delta t_n f(x_{n},t_n)$
 *    \State {\it appAction.execute(solutionHistory, stepper, AFTER\_EXPLICIT\_EVAL)}
 *    \State Correct velocity 
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperNewmarkExplicitAFormAppAction
{
public:

  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,     ///< At the beginning of the step.
    BEFORE_EXPLICIT_EVAL,   ///< Before the explicit evaluation.
    AFTER_EXPLICIT_EVAL,    ///< After the explicit evaluation.
    END_STEP        ///< At the end of the step.
  };

  /// Constructor
  StepperNewmarkExplicitAFormAppAction(){}

  /// Destructor
  virtual ~StepperNewmarkExplicitAFormAppAction(){}

  /// Execute application action for NewmarkExplicitAForm Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperNewmarkExplicitAForm<Scalar> > stepper,
    const typename StepperNewmarkExplicitAFormAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperNewmarkExplicitAFormAppAction_hpp
