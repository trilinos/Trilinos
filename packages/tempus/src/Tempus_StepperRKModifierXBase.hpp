//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperRKModifierXBase_hpp
#define Tempus_StepperRKModifierXBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperRKAppAction.hpp"

#include "Teuchos_SerialDenseVector.hpp"

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
 *  The locations of the StepperRKModifierXBase::MODIFIER_TYPE
 *  which correspond to the RK AppActions
 *  (StepperRKAppAction::ACTION_LOCATION)
 *  in takeStep are documented in each of the RK Algorithm sections:
 *  StepperExplicitRK, StepperDIRK and StepperIMEX_RK.
 */
template <class Scalar>
class StepperRKModifierXBase
  : virtual public Tempus::StepperRKAppAction<Scalar> {
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

    MODIFIER_TYPE modType                     = X_BEGIN_STEP;
    const int stageNumber                     = stepper->getStageNumber();
    Teuchos::SerialDenseVector<int, Scalar> c = stepper->getTableau()->c();
    RCP<SolutionState<Scalar> > workingState  = sh->getWorkingState();
    const Scalar dt                           = workingState->getTimeStep();
    Scalar time                               = sh->getCurrentState()->getTime();
    if (stageNumber >= 0) time += c(stageNumber) * dt;
    RCP<Thyra::VectorBase<Scalar> > x = workingState->getX();

    switch (actLoc) {
      case StepperRKAppAction<Scalar>::BEGIN_STEP: {
        modType = X_BEGIN_STEP;
        break;
      }
      case StepperRKAppAction<Scalar>::BEGIN_STAGE: {
        modType = X_BEGIN_STAGE;
        break;
      }
      case StepperRKAppAction<Scalar>::BEFORE_SOLVE: {
        modType = X_BEFORE_SOLVE;
        break;
      }
      case StepperRKAppAction<Scalar>::AFTER_SOLVE: {
        modType = X_AFTER_SOLVE;
        break;
      }
      case StepperRKAppAction<Scalar>::BEFORE_EXPLICIT_EVAL: {
        modType = X_BEFORE_EXPLICIT_EVAL;
        break;
      }
      case StepperRKAppAction<Scalar>::END_STAGE: {
        modType = X_END_STAGE;
        break;
      }
      case StepperRKAppAction<Scalar>::END_STEP: {
        modType = X_END_STEP;
        time    = workingState->getTime();
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
    X_BEGIN_STEP,            ///< Modify \f$x\f$ at the beginning of the step.
    X_BEGIN_STAGE,           ///< Modify \f$x\f$ at the beginning of the stage.
    X_BEFORE_SOLVE,          ///< Modify \f$x\f$ before the implicit solve.
    X_AFTER_SOLVE,           ///< Modify \f$x\f$ after the implicit solve.
    X_BEFORE_EXPLICIT_EVAL,  ///< Modify \f$x\f$ before the explicit evaluation.
    X_END_STAGE,             ///< Modify \f$x\f$ at the end of the stage.
    X_END_STEP               ///< Modify \f$x\f$ at the end of the step.
  };

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
                      const Scalar /* time */, const Scalar /* dt */,
                      const int /* stageNumber */,
                      const MODIFIER_TYPE modType) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperRKModifierXBase_hpp
