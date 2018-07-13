// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperForwardEulerObserver_hpp
#define Tempus_StepperForwardEulerObserver_hpp

#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this Observer <--> Stepper)
template<class Scalar> class StepperForwardEuler;

/** \brief StepperForwardEulerObserver class for StepperForwardEuler.
 *
 * This is a means for application developers to perform tasks
 * during the time steps, e.g.,
 *   - Compute specific quantities
 *   - Output information
 *   - "Massage" the working solution state
 *   - ...
 *
 * <b>Design Considerations</b>
 *   - StepperForwardEulerObserver is not stateless!  Developers may touch the
 *     solution state!  Developers need to be careful not to break the
 *     restart (checkpoint) capability.
 */
template<class Scalar>
class StepperForwardEulerObserver
 : virtual public Tempus::StepperObserver<Scalar>
{
public:

  /// Constructor
  StepperForwardEulerObserver(){}

  /// Destructor
  virtual ~StepperForwardEulerObserver(){}

  /// Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper){}

  /// Observe Stepper before Explicit ME evaluation.
  virtual void observeBeforeExplicit(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperForwardEuler<Scalar> & stepperFE){}

  /// Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper){}
};
} // namespace Tempus
#endif // Tempus_StepperForwardEulerObserver_hpp
