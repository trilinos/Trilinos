// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRKObserver_hpp
#define Tempus_StepperExplicitRKObserver_hpp

#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this Observer <--> Stepper)
template<class Scalar> class StepperExplicitRK;

/** \brief StepperExplicitRKObserver class for StepperExplicitRK.
 *
 * This is a means for application developers to perform tasks
 * during the time steps, e.g.,
 *   - Compute specific quantities
 *   - Output information
 *   - "Massage" the working solution state
 *   - ...
 *
 * <b>Design Considerations</b>
 *   - StepperExplicitRKObserver is not stateless!  Developers may touch the
 *     solution state!  Developers need to be careful not to break the
 *     restart (checkpoint) capability.
 */
template<class Scalar>
class StepperExplicitRKObserver
 : virtual public Tempus::StepperObserver<Scalar>
{
public:

  /// Constructor
  StepperExplicitRKObserver(){}

  /// Destructor
  virtual ~StepperExplicitRKObserver(){}

  /// Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper){}

  /// Observe Stepper at beginning of each stage.
  virtual void observeBeginStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperExplicitRK<Scalar> & stepperExplicitRK){}

  /// Observe Stepper before Explicit evaluation of Implicit ODE ME.
  virtual void observeBeforeExplicit(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperExplicitRK<Scalar> & stepperExplicitRK){}

  /// Observe Stepper at end of each stage.
  virtual void observeEndStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperExplicitRK<Scalar> & stepperExplicitRK){}

  /// Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper){}
};
} // namespace Tempus
#endif // Tempus_StepperExplicitRKObserver_hpp
