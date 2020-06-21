// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKModifierBase_hpp
#define Tempus_StepperRKModifierBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperRKBase.hpp"
#include "Tempus_StepperRKAppAction.hpp"


namespace Tempus {

/** \brief Base modifier for StepperRK.
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
 *  derived from the default modifier (i.e., USER RKWARE!!).
 *
 *  Below is the RK algorithm with the locations of the modify calls
 *  italicized.
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{Backward Euler with modify calls indicated.}
 *  \begin{algorithmic}[1]
 *    \State \quad {\it modifier.modify(solutionHistory, stepper, RKGIN\_STEP)}
 *    \State Compute the predictor (e.g., apply stepper to $x_n$).
 *    \State \quad {\it modifier.modify(solutionHistory, stepper, RKFORE\_SOLVE)}
 *    \State Solve $\mathcal{F}_n(\dot{x}=(x_n-x_{n-1})/\Delta t_n, x_n, t_n)=0$ for $x_n$
 *    \State \quad {\it modifier.modify(solutionHistory, stepper, AFTER\_SOLVE)}
 *    \State $\dot{x}_n \leftarrow (x_n-x_{n-1})/\Delta t_n$
 *    \State \quad {\it modifier.modify(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperRKModifierBase
  : virtual public Tempus::StepperRKAppAction<Scalar>
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
    Teuchos::RCP<StepperRKBase<Scalar> > stepper,
    const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc)
  { this->modify(sh, stepper, actLoc); }

public:

  /// Modify RK Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperRKBase<Scalar> > /* stepper */,
    const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;

};

} // namespace Tempus

#endif // Tempus_StepperRKModifierBase_hpp
