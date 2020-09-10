// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperTrapezoidalAppAction_hpp
#define Tempus_StepperTrapezoidalAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this AppAction <--> Stepper)
template<class Scalar> class StepperTrapezoidal;

/** \brief Application Action for StepperTrapezoidal.
*
*  This class provides a means to apply various actions with the Trapezoidal time step.
*  The data available to this class is solution variables (through
*  SolutionHistory), and stepper data (through the Stepper).  It allows
*  the application to just observe this data (i.e., use but not change the
*  data) to change any of it (USER BEWARE!).
*
*  Below is the Trapezoidal algorithm and includes the locations where the
*  application can take actions (in italicized).
*
*  \f{algorithm}{                                                                      
*  \renewcommand{\thealgorithm}{}                                                      
*  \caption{Trapezoidal stepper with application-action locations indicated.}               
*  \begin{algorithmic}[1]                                                              
*    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}             
*    \State Compute $y \leftarrow x_{n-1} + \frac{\Delta t_{n-1}}{2}f(x_{n-1},t_{n-1})$     
*    \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}           
*    \State Solve $x - y - \frac{\Delta t_{n-1}}{2}f(x,t_{n}) = 0$ for $x$              
*    \State {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}            
*    \State $x_n \leftarrow x$ and $\dot x_n \leftarrow f(x,t_n)$                            
*    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}               
*  \end{algorithmic}                                                                   
*  \f} 
*/
template<class Scalar>
class StepperTrapezoidalAppAction
{
public:

  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,     ///< At the beginning of the step.
    BEFORE_SOLVE,   ///< Before the solve
    AFTER_SOLVE,    ///< After then solve
    END_STEP        ///< At the end of the step.
  };

  /// Constructor
  StepperTrapezoidalAppAction(){}

  /// Destructor
  virtual ~StepperTrapezoidalAppAction(){}

  /// Execute application action for Trapezoidal Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperTrapezoidal<Scalar> > stepper,
    const typename StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperTrapezoidalAppAction_hpp
