// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBDF2ModifierBase_hpp
#define Tempus_StepperBDF2ModifierBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBDF2AppAction.hpp"

namespace Tempus {

/** \brief Base modifier for StepperBDF2.
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
 *  \caption{BDF2 with the locations of the modifier modifications indicated.  Note                                     
 *  that the following algorithm in only applied after the first two steps (where                                    
 *  the appAction for the start-up stepper is used)}                                                                 
 *  \begin{algorithmic}[1]                                                                                           
 *    \State {\it modifier.modify(solutionHistory, stepper, BEGIN\_STEP)}                                          
 *    \State Set old values of $x_{n}$, $x_{n-1}$ to the new values of $x_{n-1}$, $x_{n-2}$ respectively.            
 *    \State {\it modifier.modify(solutionHistory, stepper, BEFORE\_SOLVE)}                                        
 *    \State Solve $F( (3 x_{n} - 4 x _{n-1} + x_{n-2})/(3\Delta t),x_{n},t_{n}) for $x_n$                           
 *    \State {\it modifier.modify(solutionHistory, stepper, AFTER\_SOLVE)}                                         
 *    \State $\dot{x}_{n} \leftarrow (3 x_{n} - 4 x_{n-1} + x_{n-2})/(3\Delta t)$                                    
 *    \State {\it modifier.modify(solutionHistory, stepper, END\_STEP)}                                            
 *  \end{algorithmic}                                                                                                
 *  \f}                                                                                                              
 */
template<class Scalar>
class StepperBDF2ModifierBase
  : virtual public Tempus::StepperBDF2AppAction<Scalar>
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
    Teuchos::RCP<StepperBDF2<Scalar> > stepper,
    const typename StepperBDF2AppAction<Scalar>::ACTION_LOCATION actLoc)
  { this->modify(sh, stepper, actLoc); }

public:

  /// Modify BDF2 Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperBDF2<Scalar> > /* stepper */,
    const typename StepperBDF2AppAction<Scalar>::ACTION_LOCATION actLoc) = 0;

};

} // namespace Tempus

#endif // Tempus_StepperBDF2ModifierBase_hpp
