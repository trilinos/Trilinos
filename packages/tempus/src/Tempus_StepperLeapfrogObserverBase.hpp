// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperLeapfrogObserverBase_hpp
#define Tempus_StepperLeapfrogObserverBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperLeapfrogAppAction.hpp"


namespace Tempus {

/** \brief Base observer for StepperLeapfrog.
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
 *  Below is the Leapfrog algorithm with the locations of the observe calls
 *  italicized.
 *
 *  \f{algorithm}{                                                                               
 *  \renewcommand{\thealgorithm}{}                                                               
 *  \caption{Leapfrog with the locations of the application actions indicated.}             
 *  \begin{algorithmic}[1]
 *    \State \quad {\it observer.observe(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State Compute $\dot{x}_{n+1/2} = \dot{x}_n + 0.5\Delta t \ddot{x}_n$
 *    \State \quad {\it observer.observe(solutionHistory, stepper, BEFORE\_X\_UPDATE)}
 *    \State Compute $x_{n+1} = x_n + \Delta t \dot{x}_{n+1/2}$
 *    \State \quad {\it observer.observe(solutionHistory, stepper, BEFORE\_EXPLICIT\_EVAL)}
 *    \State Evaluate $\ddot{x}_{n+1} = f(x_{n+1},t_{n+1})$
 *    \State \quad {\it observer.observe(solutionHistory, stepper, BEFORE\_XDOT\_UPDATE)}
 *    \State Compute half-step sync $\dot{x}_{n+1} = \dot{x}_{n+1/2} + 0.5 \Delta t \ddot{x}_{n+1}$ or full step $\dot{x}_{n+3/2} = \dot{x}_{n+1/2} + \Delta t \ddot{x}_{n+1}$
 *  \end{algorithmic}                                                                            
 *  \f}                                                                                          
 */
template<class Scalar>
class StepperLeapfrogObserverBase
  : virtual public Tempus::StepperLeapfrogAppAction<Scalar>
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
    Teuchos::RCP<StepperLeapfrog<Scalar> > stepper,
    const typename StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION actLoc)
  { this->observe(sh, stepper, actLoc); }

public:

  /// Observe Leapfrog Stepper.
  virtual void observe(
    Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<const StepperLeapfrog<Scalar> > /* stepper */,
    const typename StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;

};

} // namespace Tempus

#endif // Tempus_StepperLeapfrogObserverBase_hpp
