//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TutorialObserver_hpp
#define Tempus_TutorialObserver_hpp

#include <iomanip>
#include <iostream>

#include "Teuchos_RCP.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Tempus_SolutionHistory.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_TimeStepControl.hpp"

namespace Tempus_Tutorial {

/** \brief Tutorial-local observer used by Example 7.
 *
 *  This class mirrors the callback structure of \ref Tempus::IntegratorObserver,
 *  but it does not derive from that interface.  The production
 *  \ref Tempus::IntegratorObserver interface receives a
 *  \ref Tempus::Integrator, whereas this tutorial observer receives the
 *  objects owned by the application-managed loop directly:
 *
 *  - \ref Tempus::SolutionHistory
 *  - \ref Tempus::Stepper
 *  - \ref Tempus::TimeStepControl
 *  - \ref Tempus::Status
 *
 *  This keeps Example 7 focused on the observer pattern without introducing
 *  \ref Tempus::IntegratorBasic until the next tutorial step.
 */
template <class Scalar>
class TutorialObserver {
 public:
  using SolutionHistory = Tempus::SolutionHistory<Scalar>;
  using Stepper         = Tempus::Stepper<Scalar>;
  using TimeStepControl = Tempus::TimeStepControl<Scalar>;

  void observeStartIntegrator(
      const Teuchos::RCP<SolutionHistory>& solHistory,
      const Teuchos::RCP<Stepper>& stepper,
      const Teuchos::RCP<TimeStepControl>& timeStepControl,
      Tempus::Status& integratorStatus)
  {
    (void)solHistory;
    (void)integratorStatus;

    // Output
    std::cout << std::fixed;
    std::cout << std::setw(8)  << "index"
              << std::setw(10) << "time"
              << std::setw(12) << "x_0"
              << std::setw(12) << "x_1" << std::endl;

    auto currentState = solHistory->getCurrentState();
    std::cout << std::setw(8)  << currentState->getIndex()
              << std::setw(10) << std::setprecision(3) << currentState->getTime()
              << std::setw(12) << std::setprecision(4) << Thyra::get_ele(*(currentState->getX()), 0)
              << std::setw(12) << std::setprecision(4) << Thyra::get_ele(*(currentState->getX()), 1)
              << std::endl;
  }

  void observeStartTimeStep(
      const Teuchos::RCP<SolutionHistory>& solHistory,
      const Teuchos::RCP<Stepper>& stepper,
      const Teuchos::RCP<TimeStepControl>& timeStepControl,
      Tempus::Status& integratorStatus)
  {
    (void)solHistory;
    (void)stepper;
    (void)timeStepControl;
    (void)integratorStatus;
  }

  void observeNextTimeStep(
      const Teuchos::RCP<SolutionHistory>& solHistory,
      const Teuchos::RCP<Stepper>& stepper,
      const Teuchos::RCP<TimeStepControl>& timeStepControl,
      Tempus::Status& integratorStatus)
  {
    (void)solHistory;
    (void)stepper;
    (void)timeStepControl;
    (void)integratorStatus;
  }

  void observeBeforeTakeStep(
      const Teuchos::RCP<SolutionHistory>& solHistory,
      const Teuchos::RCP<Stepper>& stepper,
      const Teuchos::RCP<TimeStepControl>& timeStepControl,
      Tempus::Status& integratorStatus)
  {
    (void)solHistory;
    (void)stepper;
    (void)timeStepControl;
    (void)integratorStatus;
  }

  void observeAfterTakeStep(
      const Teuchos::RCP<SolutionHistory>& solHistory,
      const Teuchos::RCP<Stepper>& stepper,
      const Teuchos::RCP<TimeStepControl>& timeStepControl,
      Tempus::Status& integratorStatus)
  {
    (void)solHistory;
    (void)stepper;
    (void)timeStepControl;
    (void)integratorStatus;
  }

  void observeAfterCheckTimeStep(
      const Teuchos::RCP<SolutionHistory>& solHistory,
      const Teuchos::RCP<Stepper>& stepper,
      const Teuchos::RCP<TimeStepControl>& timeStepControl,
      Tempus::Status& integratorStatus)
  {
    (void)stepper;
    (void)timeStepControl;

    auto workingState = solHistory->getWorkingState();

    if (workingState->getSolutionStatus() == Tempus::Status::FAILED) {
      integratorStatus = Tempus::Status::FAILED;
    }
  }

  void observeEndTimeStep(
      const Teuchos::RCP<SolutionHistory>& solHistory,
      const Teuchos::RCP<Stepper>& stepper,
      const Teuchos::RCP<TimeStepControl>& timeStepControl,
      Tempus::Status& integratorStatus)
  {
    (void)stepper;
    (void)timeStepControl;
    (void)integratorStatus;

    auto currentState = solHistory->getCurrentState();

    // Output
    if (solHistory->getCurrentState()->getIndex() % 100 == 0) {
      currentState = solHistory->getCurrentState();
      std::cout << std::setw(8)  << currentState->getIndex()
                << std::setw(10) << std::setprecision(3) << currentState->getTime()
                << std::setw(12) << std::setprecision(4) << Thyra::get_ele(*(currentState->getX()), 0)
                << std::setw(12) << std::setprecision(4) << Thyra::get_ele(*(currentState->getX()), 1)
                << std::endl;
    }
  }

  void observeEndIntegrator(
      const Teuchos::RCP<SolutionHistory>& solHistory,
      const Teuchos::RCP<Stepper>& stepper,
      const Teuchos::RCP<TimeStepControl>& timeStepControl,
      Tempus::Status& integratorStatus)
  {
    (void)solHistory;
    (void)stepper;
    (void)timeStepControl;
    (void)integratorStatus;
  }
};

}  // namespace Tempus_Tutorial

#endif  // Tempus_TutorialObserver_hpp
