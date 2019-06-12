// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKObserver_hpp
#define Tempus_StepperRKObserver_hpp

#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperObserver.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this Observer <--> Stepper)
//template<class Scalar> class Stepper;

/** \brief StepperRKObserver class for StepperRK.
 *
 * This is a means for application developers to perform tasks
 * during the time steps, e.g.,
 *   - Compute specific quantities
 *   - Output information
 *   - "Massage" the working solution state
 *   - ...
 *
 * <b>Design Considerations</b>
 *   - StepperRKObserver is not stateless!  Developers may touch the
 *     solution state!  Developers need to be careful not to break the
 *     restart (checkpoint) capability.
 */
template<class Scalar>
class StepperRKObserver
 : virtual public Tempus::StepperObserver<Scalar>
{
public:

  /// Constructor
  StepperRKObserver(){}

  /// Destructor
  virtual ~StepperRKObserver(){}

  /// 1.) Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}

  /// 2.) Observe Stepper at beginning of each stage.
  virtual void observeBeginStage(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}

  /// 3.) Observe Stepper before Explicit evaluation of Implicit ODE ME (IMEX).
  virtual void observeBeforeImplicitExplicitly(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}

  /// 4.) Observe Stepper before nonlinear solve (DIRK/IMEX).
  virtual void observeBeforeSolve(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}

  /// 5.) Observe Stepper after nonlinear solve (DIRK/IMEX).
  virtual void observeAfterSolve(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}

  /// 6.) Observe Stepper before Explicit evaluation of Implicit ODE ME (IMEX).
  virtual void observeBeforeExplicit(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}

  /// 7.) Observe Stepper at end of each stage.
  virtual void observeEndStage(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}

  /// 8.) Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}
};
} // namespace Tempus
#endif // Tempus_StepperRKObserver_hpp
