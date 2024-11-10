//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperLeapfrogModifierXBase_hpp
#define Tempus_StepperLeapfrogModifierXBase_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperLeapfrogAppAction.hpp"

namespace Tempus {

/** \brief Base ModifierX for StepperLeapfrog.
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
 *  The locations of the StepperLeapfrogModifierXBase::MODIFIER_TYPE
 *  which correspond to the AppAction calls
 *  (StepperLeapfrogAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperLeapfrog.
 */

template <class Scalar>
class StepperLeapfrogModifierXBase
  : virtual public Tempus::StepperLeapfrogAppAction<Scalar> {
 private:
  /* \brief Adaptor execute function
   *
   *  This is an adaptor function to bridge between the AppAction
   *  interface and the ModifierX interface.  It is meant to be private
   *  and non-virtual as deriving from this class should only need to
   *  implement the modify function.
   *
   *  For the ModifierX interface, this adaptor maps the
   *  StepperLeapfrogAppAction::ACTION_LOCATION to the
   *  StepperLeapfrogModifierX::MODIFIERX_TYPE, and only pass the solution
   *  (\f$x\f$ and/or \f$\dot{x}\f$ and other parameters to the modify
   *  function.
   */
  void execute(
      Teuchos::RCP<SolutionHistory<Scalar> > sh,
      Teuchos::RCP<StepperLeapfrog<Scalar> > stepper,
      const typename StepperLeapfrogAppAction<Scalar>::ACTION_LOCATION actLoc)
  {
    using Teuchos::RCP;

    MODIFIER_TYPE modType                    = X_BEGIN_STEP;
    RCP<SolutionState<Scalar> > workingState = sh->getWorkingState();
    const Scalar time                        = workingState->getTime();
    const Scalar dt                          = workingState->getTimeStep();
    RCP<Thyra::VectorBase<Scalar> > x;

    switch (actLoc) {
      case StepperLeapfrogAppAction<Scalar>::BEGIN_STEP: {
        modType = X_BEGIN_STEP;
        x       = workingState->getX();
        break;
      }
      case StepperLeapfrogAppAction<Scalar>::BEFORE_X_UPDATE: {
        modType = X_BEFORE_X_UPDATE;
        x       = workingState->getX();
        break;
      }
      case StepperLeapfrogAppAction<Scalar>::BEFORE_EXPLICIT_EVAL: {
        modType = X_BEFORE_EXPLICIT_EVAL;
        x       = workingState->getX();
        break;
      }
      case StepperLeapfrogAppAction<Scalar>::BEFORE_XDOT_UPDATE: {
        modType = X_BEFORE_XDOT_UPDATE;
        x       = workingState->getX();
        break;
      }
      case StepperLeapfrogAppAction<Scalar>::END_STEP: {
        modType = X_END_STEP;
        x       = workingState->getX();
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown action location.\n");
    }

    this->modify(x, time, dt, modType);
  }

 public:
  /// Indicates the location of application action (see algorithm).
  enum MODIFIER_TYPE {
    X_BEGIN_STEP,            ///< Modify \f$x\f$ at the beginning of the step.
    X_BEFORE_X_UPDATE,       ///< Modify \f$x\f$ before updating x
    X_BEFORE_EXPLICIT_EVAL,  ///< Modify \f$x\f$ before the explicit ME
                             ///< evaluation
    X_BEFORE_XDOT_UPDATE,    ///< Modify \f$x\f$ Before updating xDot
    X_END_STEP               ///< Modify \f$\dot{x}\f$ at the end of the step.
  };

  /// Modify solution based on the MODIFIER_TYPE.
  virtual void modify(Teuchos::RCP<Thyra::VectorBase<Scalar> > /* x */,
                      const Scalar /* time */, const Scalar /* dt */,
                      const MODIFIER_TYPE modType) = 0;
};

}  // namespace Tempus

#endif  // Tempus_StepperLeapfrogModifierXBase_hpp
