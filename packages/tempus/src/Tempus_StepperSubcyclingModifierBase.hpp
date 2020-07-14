// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperSubcyclingModifierBase_hpp
#define Tempus_StepperSubcyclingModifierBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperSubcyclingAppAction.hpp"

namespace Tempus {

/** \brief Base modifier for StepperSubcycling.
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
 *  \caption{Subcycling with the locations of the application actions indicated.}           
 *  \begin{algorithmic}[1]                                                                     
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}                    
 *    \State Compute $x_{n}$ from $x_{n-1}$, applying appActions from sub-steppers                       
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)} 
 *  \end{algorithmic}                                                                         
 *  \f}                                                                                       
 */

template<class Scalar>
class StepperSubcyclingModifierBase
  : virtual public Tempus::StepperSubcyclingAppAction<Scalar>
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
    Teuchos::RCP<StepperSubcycling<Scalar> > stepper,
    const typename StepperSubcyclingAppAction<Scalar>::ACTION_LOCATION actLoc)
  { this->modify(sh, stepper, actLoc); }

public:

  /// Modify Subcycling Stepper.
  virtual void modify(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Teuchos::RCP<StepperSubcycling<Scalar> > /* stepper */,
    const typename StepperSubcyclingAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;

};

} // namespace Tempus

#endif // Tempus_StepperSubcyclingModifierBase_hpp
