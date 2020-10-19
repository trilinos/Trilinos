// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBDF2AppAction_hpp
#define Tempus_StepperBDF2AppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBDF2.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this AppAction <--> Stepper)
template<class Scalar> class StepperBDF2;

/** \brief Application Action for StepperBDF2.
 *
 *  This class provides a means to apply various actions with the BDF2 time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data (i.e., use but not change the
 *  data) to change any of it (USER BEWARE!).
 *
 *  Below is the BDF2 algorithm and includes the locations where the
 *  application can take actions (in italicized).
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{BDF2 with the locations of the application actions indicated.  Note
 *  that the following algorithm in only applied after the first two steps (where
 *  the appAction for the start-up stepper is used)}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}                                          
 *    \State Set old values of $x_{n}$, $x_{n-1}$ to the new values of $x_{n-1}$, $x_{n-2}$ respectively. 
 *    \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}                                        
 *    \State Solve $F( (3 x_{n} - 4 x _{n-1} + x_{n-2})/(3\Delta t),x_{n},t_{n}) for $x_n$
 *    \State {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}                                         
 *    \State $\dot{x}_{n} \leftarrow (3 x_{n} - 4 x_{n-1} + x_{n-2})/(3\Delta t)$                                    
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}  
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperBDF2AppAction
{
public:

  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,    ///< At the beginning of the step.
    BEFORE_SOLVE,  ///< Before the solve.
    AFTER_SOLVE,   ///<  After the solve.
    END_STEP,      ///< At the end of the step.
  };

  /// Constructor
  StepperBDF2AppAction(){}

  /// Destructor
  virtual ~StepperBDF2AppAction(){}

  /// Execute application action for BDF2 Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperBDF2<Scalar> > stepper,
    const typename StepperBDF2AppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperBDF2AppAction_hpp
