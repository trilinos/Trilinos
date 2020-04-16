// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKAppAction_hpp
#define Tempus_StepperRKAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperRKBase.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this AppAction <--> Stepper)
template<class Scalar> class StepperRKBase;

/** \brief Application Action for StepperRKBase.
 *
 *  This class provides a means to apply various actions with the RK time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data (i.e., use but not change the
 *  data) to change any of it (USER BEWARE!).
 *
 *  Below is the IMEX RK algorithm and includes the locations where the
 *  application can take actions (in italicized).
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{IMEX RK with the application-action locations indicated.}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State $X \leftarrow x_{n-1}$ \Comment Set initial guess to last timestep.
 *    \For {$i = 0 \ldots s-1$}
 *      \If { $\hat{a}_{k,i} = 0 \;\forall k = (i+1,\ldots, s-1)$, $\hat{b}(i) = 0$, $\hat{b}^\ast(i) = 0$ and \newline \hspace*{0.37in}
 *            $a_{k,i} = 0 \;\forall k = (i+1,\ldots, s-1)$, $b(i) = 0$, $b^\ast(i) = 0$}
 *        \State $g_i \leftarrow 0$ \Comment{Not needed for later calculations.}
 *        \State $f_i \leftarrow 0$ \Comment{Not needed for later calculations.}
 *        \State {\bf continue}
 *      \EndIf
 *      \State $\tilde{X} \leftarrow x_{n-1} - \Delta t\,\sum_{j=1}^{i-1} \left(
 *            \hat{a}_{ij}\, f_j + a_{ij}\, g_j \right)$
 *      \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STAGE)}
 *      \State \Comment Implicit Tableau
 *      \If {$a_{ii} = 0$}
 *        \State $X \leftarrow \tilde{X}$
 *        \If {$i=0$ and ``Use FSAL''} \Comment{Save an evaluation.}
 *          \State $g_0 \leftarrow g_{s-1}$
 *                          \Comment{Use $g_{s-1}$ from $n-1$ time step.}
 *        \Else
 *          \State $g_i \leftarrow M(X, t_i)^{-1}\, G(X, t_i)$
 *        \EndIf
 *      \Else
 *        \If {``Zero initial guess.''}
 *          \State $X \leftarrow 0$
 *            \Comment{Else use previous stage value as initial guess.}
 *        \EndIf
 *        \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}
 *        \State Solve $\mathcal{G}\left(\tilde{\dot{X}}
 *            = \frac{X-\tilde{X}}{a_{ii} \Delta t},X,t_i\right) = 0$ for $X$
 *        \State {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *        \State $\tilde{\dot{X}} \leftarrow \frac{X - \tilde{X}}{a_{ii} \Delta t}$
 *        \State $g_i \leftarrow - \tilde{\dot{X}}$
 *      \EndIf
 *      \State \Comment Explicit Tableau
 *      \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_EXPLICIT\_EVAL)}
 *      \State $f_i \leftarrow M(X,\hat{t}_i)^{-1}\, F(X,\hat{t}_i)$
 *      \State $\dot{X} \leftarrow - g_i - f_i$ [Optionally]
 *      \State {\it appAction.execute(solutionHistory, stepper, END\_STAGE)}
 *    \EndFor
 *    \State $x_n \leftarrow x_{n-1} - \Delta t\,\sum_{i=1}^{s}\hat{b}_i\,f_i
 *                                   - \Delta t\,\sum_{i=1}^{s}     b_i \,g_i$
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperRKAppAction
{
public:

  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,           ///< At the beginning of the step.
    BEGIN_STAGE,          ///< At the beginning of the stage.
    BEFORE_SOLVE,         ///< Before the implicit solve.
    AFTER_SOLVE,          ///< After the implicit solve.
    BEFORE_EXPLICIT_EVAL, ///< Before the explicit evaluation.
    END_STAGE,            ///< At the end of the stage.
    END_STEP              ///< At the end of the step.
  };

  /// Constructor
  StepperRKAppAction(){}

  /// Destructor
  virtual ~StepperRKAppAction(){}

  /// Execute application action for RK Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperRKBase<Scalar> > stepper,
    const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperRKAppAction_hpp
