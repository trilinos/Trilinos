// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBDF2ObserverBase_hpp
#define Tempus_StepperBDF2ObserverBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperBDF2AppAction.hpp"


namespace Tempus {

/** \brief Base observer for StepperBDF2.
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
 *  Below is the BDF2 algorithm with the locations of the observe calls
 *  italicized.
 *
 *  \f{algorithm}{                                                                               
 *  \renewcommand{\thealgorithm}{}                                                               
 *  \caption{BDF2 with the locations of the oberserver observations  indicated.   Note                             
 *  that the following algorithm in only applied after the first two steps (where                                    
 *  the appAction for the start-up stepper is used)}                                                                 
 *  \begin{algorithmic}[1]                                                                                           
 *    \State {\it observer.observe(solutionHistory, stepper, BEGIN\_STEP)}                                          
 *    \State Set old values of $x_{n}$, $x_{n-1}$ to the new values of $x_{n-1}$, $x_{n-2}$ respectively.            
 *    \State {\it oberserver.observe(solutionHistory, stepper, BEFORE\_SOLVE)}                                       
 *    \State Solve $F( (3 x_{n} - 4 x _{n-1} + x_{n-2})/(3\Delta t),x_{n},t_{n}) for $x_n$                           
 *    \State {\it observer.observe(solutionHistory, stepper, AFTER\_SOLVE)}                                         
 *    \State $\dot{x}_{n} \leftarrow (3 x_{n} - 4 x_{n-1} + x_{n-2})/(3\Delta t)$                                    
 *    \State {\it observer.observe(solutionHistory, stepper, END\_STEP)}                                            
 *  \end{algorithmic}   
 * \f}
 */ 
template<class Scalar>
class StepperBDF2ObserverBase
  : virtual public Tempus::StepperBDF2AppAction<Scalar>
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
    Teuchos::RCP<StepperBDF2<Scalar> > stepper,
    const typename StepperBDF2AppAction<Scalar>::ACTION_LOCATION actLoc)
  { this->observe(sh, stepper, actLoc); }

public:

  /// Observe BDF2 Stepper.
  virtual void observe(
    Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<const StepperBDF2<Scalar> > /* stepper */,
    const typename StepperBDF2AppAction<Scalar>::ACTION_LOCATION actLoc) = 0;

};

} // namespace Tempus

#endif // Tempus_StepperBDF2ObserverBase_hpp
