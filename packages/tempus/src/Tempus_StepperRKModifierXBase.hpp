// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKModifierXBase_hpp
#define Tempus_StepperRKModifierXBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperRKAppAction.hpp"


namespace Tempus {

/** \brief Base ModifierX for StepperRK.
 *
 *  This class provides a means to modify just the solution values
 *  (i.e., \f$x\f$ and \f$dot{x}\f$), and nothing else, but time and
 *  timestep are also provided.
 *
 *  Users deriving from this class can access and change the solution
 *  during the timestep (e.g., limiting the solution for monoticity).
 *  It is expected that the user knows what changes are allowable without
 *  affecting the Stepper correctness, performance, accuracy and stability
 *  (i.e., USER BEWARE!!).
 *
 *  Below is the RK algorithm with the locations of the ModifierX calls
 *  italicized.
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{Backward Euler with modify calls indicated.}
 *  \begin{algorithmic}[1]
 *    \State \quad {\it modifierX.modify(x, time, dt, X\_BEGIN\_STEP)}
 *    \State Compute the predictor (e.g., apply stepper to $x_n$).
 *    \State \quad {\it modifierX.modify(x, time, dt, X\_BEFORE\_SOLVE)}
 *    \State Solve $\mathcal{F}_n(\dot{x}=(x_n-x_{n-1})/\Delta t_n, x_n, t_n)=0$ for $x_n$
 *    \State \quad {\it modifierX.modify(x, time, dt, X\_AFTER\_SOLVE)}
 *    \State $\dot{x}_n \leftarrow (x_n-x_{n-1})/\Delta t_n$
 *    \State \quad {\it modifierX.modify(x, time, dt, XDOT\_END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 */
template<class Scalar>
class StepperRKModifierXBase
  : virtual public Tempus::StepperRKAppAction<Scalar>
{
private:

  /* \brief Adaptor execute function
   *
   *  This is an adaptor function to bridge between the AppAction
   *  interface and the ModifierX interface.  It is meant to be private
   *  and non-virtual as deriving from this class should only need to
   *  implement the modify function.
   *
   *  For the ModifierX interface, this adaptor maps the
   *  StepperRKAppAction::ACTION_LOCATION to the
   *  StepperRKModifierX::MODIFIERX_TYPE, and only pass the solution
   *  (\f$x\f$ and/or \f$\dot{x}\f$ and other parameters to the modify
   *  function.
   */
  void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperRKBase<Scalar> > stepper,
    const typename StepperRKAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    using Teuchos::RCP;

    MODIFIER_TYPE modType = X_BEGIN_STEP;
    RCP<SolutionState<Scalar> > workingState = sh->getWorkingState();
    const Scalar time = workingState->getTime();
    const Scalar dt   = workingState->getTimeStep();
    const int stageNumber = stepper->getStageNumber();
    RCP<Thyra::VectorBase<Scalar> > x;

    switch(actLoc) {
      case StepperRKAppAction<Scalar>::BEGIN_STEP:
      {
        modType = X_BEGIN_STEP;
        x = workingState->getX();
        break;
      }
      case StepperRKAppAction<Scalar>::BEGIN_STAGE:
      {
        modType = X_BEGIN_STAGE;
        x = workingState->getX();
        break;
      }
      case StepperRKAppAction<Scalar>::BEFORE_SOLVE:
      {
        modType = X_BEFORE_SOLVE;
        x = workingState->getX();
        break;
      }
      case StepperRKAppAction<Scalar>::AFTER_SOLVE:
      {
        modType = X_AFTER_SOLVE;
        x = workingState->getX();
        break;
      }
      case StepperRKAppAction<Scalar>::BEFORE_EXPLICIT_EVAL:
      {
        modType = X_BEFORE_EXPLICIT_EVAL;
        x = stepper->getStageX();
        break;
      }
      case StepperRKAppAction<Scalar>::END_STAGE:
      {
        modType = X_END_STAGE;
        x = stepper->getStageX();
        break;
      }
      case StepperRKAppAction<Scalar>::END_STEP:
      {
        modType = X_END_STEP;
        x = workingState->getX();
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown action location.\n");
    }

    this->modify(x, time, dt, stageNumber, modType);
  }

public:

  /// Indicates the location of application action (see algorithm).
  enum MODIFIER_TYPE {
    X_BEGIN_STEP,           ///< Modify \f$x\f$ at the beginning of the step.
    X_BEGIN_STAGE,          ///< Modify \f$x\f$ at the beginning of the stage.
    X_BEFORE_SOLVE,         ///< Modify \f$x\f$ before the implicit solve.
    X_AFTER_SOLVE,          ///< Modify \f$x\f$ after the implicit solve.
    X_BEFORE_EXPLICIT_EVAL, ///< Modify \f$x\f$ before the explicit evaluation.
    X_END_STAGE,            ///< Modify \f$x\f$ at the end of the stage.
    X_END_STEP              ///< Modify \f$x\f$ at the end of the step.
  };

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(
    Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
    const Scalar /* time */, const Scalar /* dt */,
    const int /* stageNumber */,
    const MODIFIER_TYPE modType) = 0;

};

} // namespace Tempus

#endif // Tempus_StepperRKModifierXBase_hpp
