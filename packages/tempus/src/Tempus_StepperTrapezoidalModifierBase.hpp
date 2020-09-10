// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperTrapezoidalModifierBase_hpp
#define Tempus_StepperTrapezoidalModifierBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperTrapezoidalAppAction.hpp"

namespace Tempus {

/** \brief Base modifier for StepperTrapezoidal.
*
*  This class provides a means to modify values (e.g., solution variables
*  through SolutionHistory, and stepper member data through the Stepper),
*  and can be very powerful and easy to make changes to the stepper and
*  the solution.
*
*  Users deriving from this class can access a lot of data, and it is
*  expected that those users know what changes are allowable without
*  affecting the Stepper correctness, performance, accuracy and stability.
*  Thus the user should be careful when accessing data through classes
*  derived from the default modifier (i.e., USER BEWARE!!).
* 
*  \f{algorithm}{                                                             
*  \renewcommand{\thealgorithm}{}                                                          
*  \caption{Trapezoidal stepper with modfiy call locations indicated.}                    
*  \begin{algorithmic}[1]                                                                
*    \State {\it modifier.modify(solutionHistory, stepper, BEGIN\_STEP)}                 
*    \State Compute $y \leftarrow x_{n-1} + \frac{\Delta t_{n-1}}{2}f(x_{n-1},t_{n-1})$    
*    \State {\it modifier.modify(solutionHistory, stepper, BEFORE\_SOLVE)}                
*    \State Solve $x - y - \frac{\Delta t_{n-1}}{2}f(x,t_{n}) = 0$ for $x$    
*    \State {\it modifier.modify(solutionHistory, stepper, AFTER\_SOLVE)}                
*    \State $x_n \leftarrow x$ and $\dot x_n \leftarrow f(x,t_n)$                      
*    \State {\it modifier.modify(solutionHistory, stepper, END\_STEP)}           
*  \end{algorithmic}                                                                 
*  \f}                 
*/

template<class Scalar>
class StepperTrapezoidalModifierBase
  : virtual public Tempus::StepperTrapezoidalAppAction<Scalar>
{
private:

  /* \brief Adaptor execute function
   *
   *  This is an adaptor function to bridge between the AppAction
   *  interface and the Modifier interface.  It is meant to be private
   *  and non-virtual as deriving from this class should only need to
   *  implement the modify function.
   *
   *  For the Modifier interface, this adaptor is a "simple pass through".
   */
  void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperTrapezoidal<Scalar> > stepper,
    const typename StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION actLoc)
  { this->modify(sh, stepper, actLoc); }

public:

  /// Modify Trapezoidal Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperTrapezoidal<Scalar> > /* stepper */,
    const typename StepperTrapezoidalAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;

};

} // namespace Tempus

#endif // Tempus_StepperTrapezoidalModifierBase_hpp
