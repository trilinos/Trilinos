// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBDF2Observer_hpp
#define Tempus_StepperBDF2Observer_hpp

#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this Observer <--> Stepper)
template<class Scalar> class StepperBDF2;

/** \brief StepperBDF2Observer class for StepperBDF2.
 *
 * This is a means for application developers to perform tasks
 * during the time steps, e.g.,
 *   - Compute specific quantities
 *   - Output information
 *   - "Massage" the working solution state
 *   - ...
 *
 * <b>Design Considerations</b>
 *   - StepperBDF2Observer is not stateless!  Developers may touch the
 *     solution state!  Developers need to be careful not to break the
 *     restart (checkpoint) capability.
 */
template<class Scalar>
class StepperBDF2Observer
 : virtual public Tempus::StepperObserver<Scalar>
{
public:

  /// Constructor
  StepperBDF2Observer(){}

  /// Destructor
  virtual ~StepperBDF2Observer(){}

  /// Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper){}

  /// Observe Stepper before nonlinear solve.
  virtual void observeBeforeSolve(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperBDF2<Scalar> & stepperBDF2){}

  /// Observe Stepper after nonlinear solve.
  virtual void observeAfterSolve(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperBDF2<Scalar> & stepperBDF2){}

  /// Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper){}
};
} // namespace Tempus
#endif // Tempus_StepperBDF2Observer_hpp
