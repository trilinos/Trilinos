// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkImplicitAFormAppAction_hpp
#define Tempus_StepperNewmarkImplicitAFormAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this AppAction <--> Stepper)
template<class Scalar> class StepperNewmarkImplicitAForm;

/** \brief Application Action for StepperNewmarkImplicitAForm.
 *
 *  This class provides a means to apply various actions with the NewmarkImplicitAForm time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data (i.e., use but not change the
 *  data) to change any of it (USER BEWARE!).
 *
 *  Below is the NewmarkImplicitAForm algorithm and includes the locations where the
 *  application can take actions (in italicized).
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{Newmark Implicit-A with application-action locations indicated.}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State Compute the predictor (e.g., apply stepper to $x_n$).
 *    \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}
 *    \State Solve $\mathcal{F}_n(\dot{x}=(x_n-x_{n-1})/\Delta t_n, x_n, t_n)=0$ for $x_n$
 *    \State {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *    \State $\dot{x}_n \leftarrow (x_n-x_{n-1})/\Delta t_n$
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperNewmarkImplicitAFormAppAction
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
  StepperNewmarkImplicitAFormAppAction(){}

  /// Destructor
  virtual ~StepperNewmarkImplicitAFormAppAction(){}

  /// Execute application action for NewmarkImplicitAForm Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperNewmarkImplicitAForm<Scalar> > stepper,
    const typename StepperNewmarkImplicitAFormAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperNewmarkImplicitAFormAppAction_hpp
