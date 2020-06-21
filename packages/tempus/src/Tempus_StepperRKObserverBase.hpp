// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKObserverBase_hpp
#define Tempus_StepperRKObserverBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperRKAppAction.hpp"


namespace Tempus {

/** \brief Base observer for StepperRK.
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
 *  Below is the RK algorithm with the locations of the observe calls
 *  italicized.
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{DIRK Backward Euler with observe calls indicated.}
 *  \begin{algorithmic}[1]
 *    \State \quad {\it observer.observe(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State Compute the predictor (e.g., apply stepper to $x_n$).
 *    \State \quad {\it observer.observe(solutionHistory, stepper, BEFORE\_SOLVE)}
 *    \State Solve $\mathcal{F}_n(\dot{x}=(x_n-x_{n-1})/\Delta t_n, x_n, t_n)=0$ for $x_n$
 *    \State \quad {\it observer.observe(solutionHistory, stepper, AFTER\_SOLVE)}
 *    \State $\dot{x}_n \leftarrow (x_n-x_{n-1})/\Delta t_n$
 *    \State \quad {\it observer.observe(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperRKObserverBase
  : virtual public Tempus::StepperRKAppAction<Scalar>
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
    Teuchos::RCP<StepperRKBase<Scalar> > stepper,
    const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc)
  { this->observe(sh, stepper, actLoc); }

public:

  /// Observe RK Stepper.
  virtual void observe(
    Teuchos::RCP<const SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<const StepperRKBase<Scalar> > /* stepper */,
    const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;

};


} // namespace Tempus

#endif // Tempus_StepperRKObserverBase_hpp
