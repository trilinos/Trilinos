// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperSubcyclingAppAction_hpp
#define Tempus_StepperSubcyclingAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this AppAction <--> Stepper)
template<class Scalar> class StepperSubcycling;

/** \brief Application Action for StepperSubcycling.
 *
 *  This class provides a means to apply various actions with the Subcycling time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data (i.e., use but not change the
 *  data) to change any of it (USER BEWARE!).
 *
 *  Below is the Subcycling algorithm and includes the locations where the
 *  application can take actions (in italicized).
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{Subcyling with the locations of the application actions indicated.}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State Compute $x_{n}$ from $x_{n-1}$, applying appActions from sub-steppers
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperSubcyclingAppAction
{
public:

  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,     ///< At the beginning of the step.
    END_STEP        ///< At the end of the step.
  };

  /// Constructor
  StepperSubcyclingAppAction(){}

  /// Destructor
  virtual ~StepperSubcyclingAppAction(){}

  /// Execute application action for Subcycling Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperSubcycling<Scalar> > stepper,
    const typename StepperSubcyclingAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperSubcyclingAppAction_hpp
