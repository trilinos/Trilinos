// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperLeapfrogObserver_hpp
#define Tempus_StepperLeapfrogObserver_hpp

#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this Observer <--> Stepper)
template<class Scalar> class StepperLeapfrog;

/** \brief StepperLeapfrogObserver class for StepperLeapfrog.
 *
 * This is a means for application developers to perform tasks
 * during the time steps, e.g.,
 *   - Compute specific quantities
 *   - Output information
 *   - "Massage" the working solution state
 *   - ...
 *
 * <b>Design Considerations</b>
 *   - StepperLeapfrogObserver is not stateless!  Developers may touch the
 *     solution state!  Developers need to be careful not to break the
 *     restart (checkpoint) capability.
 */
template<class Scalar>
class StepperLeapfrogObserver
 : virtual public Tempus::StepperObserver<Scalar>
{
public:

  /// Constructor
  StepperLeapfrogObserver(){}

  /// Destructor
  virtual ~StepperLeapfrogObserver(){}

  /// Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper){}

  /// Observe Stepper before Explicit ME evaluation while initializing xDotDot
  virtual void observeBeforeExplicitInitialize(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperLeapfrog<Scalar> & stepperLF){}

  /// Observe Stepper before updating xDot while initializing xDotDot
  virtual void observeBeforeXDotUpdateInitialize(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperLeapfrog<Scalar> & stepperLF){}

  /// Observe Stepper before updating x
  virtual void observeBeforeXUpdate(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperLeapfrog<Scalar> & stepperLF){}

  /// Observe Stepper before Explicit ME evaluation.
  virtual void observeBeforeExplicit(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperLeapfrog<Scalar> & stepperLF){}

  /// Observe Stepper before updating xDot
  virtual void observeBeforeXDotUpdate(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperLeapfrog<Scalar> & stepperLF){}

  /// Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper){}
};
} // namespace Tempus
#endif // Tempus_StepperLeapfrogObserver_hpp
