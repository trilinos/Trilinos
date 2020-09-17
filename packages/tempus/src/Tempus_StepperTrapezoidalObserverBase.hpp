// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperTrapezoidalObserverBase_hpp
#define Tempus_StepperTrapezoidalObserverBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperTrapezoidalAppAction.hpp"


namespace Tempus {

/** \brief Base observer for StepperTrapezoidal.
*
*  This class provides a means to observe values (e.g., solution variables
*  through SolutionHistory, and stepper member data through the Stepper),
*  and cannot modify them.
*
*  Users deriving from this class can observer a lot of data, and it is
*  expected that users will NOT modify any of that data.  If the user
*  wishes to modify the solution and/or stepper data during the
*  Stepper::takeStep, they should use the Modifier class (with care!).
*
*  Below is the Trapezoidal algorithm with the locations of the observe calls
*  italicized.
*
*  \f{algorithm}{                                                                                                   
*  \renewcommand{\thealgorithm}{}                                                                                   
*  \caption{Trapezoidal stepper with observe call locations indicated.}                                       
*  \begin{algorithmic}[1]                                                                                       
*    \State {\it observer.observe(solutionHistory, stepper, BEGIN\_STEP)}                                          
*    \State Compute $y \leftarrow x_{n-1} + \frac{\Delta t_{n-1}}{2}f(x_{n-1},t_{n-1})$                           
*    \State {\it observer.observe(solutionHistory, stepper, BEFORE\_SOLVE)}                                   
*    \State Solve $x - y - \frac{\Delta t_{n-1}}{2}f(x,t_{n}) = 0$ for $x$                    
*    \State {\it observer.observe(solutionHistory, stepper, AFTER\_SOLVE)}                                    
*    \State $x_n \leftarrow x$ and $\dot x_n \leftarrow f(x,t_n)$                                               
*    \State {\it observer.observe(solutionHistory, stepper, END\_STEP)}                                       
*  \end{algorithmic}                                                                                          
*  \f}       
*/
template<class Scalar>
class StepperTrapezoidalObserverBase
  : virtual public Tempus::StepperTrapezoidalAppAction<Scalar>
{
private:

  /* \brief Adaptor execute function
   *
   *  This is an adaptor function to bridge between the AppAction
   *  interface and this derived interface.  It is meant to be private
   *  and non-virtual as deriving from this class should only need to
   *  implement the observe function.
   *
   *  For the Observer interface, this adaptor simply "applies" constantness
   *  to the arguments.
   */
  void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperTrapezoidal<Scalar> > stepper,
    const typename StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION actLoc)
  { this->observe(sh, stepper, actLoc); }

public:

  /// Observe Trapezoidal Stepper.
  virtual void observe(
    Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<const StepperTrapezoidal<Scalar> > /* stepper */,
    const typename StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;

};

} // namespace Tempus

#endif // Tempus_StepperTrapezoidalObserverBase_hpp
