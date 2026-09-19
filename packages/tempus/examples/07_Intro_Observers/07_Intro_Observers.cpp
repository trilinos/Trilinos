//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "Teuchos_StandardCatchMacros.hpp"

#include "../00_Basic_Problem/Tutorial_Regression_Tester.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"

#include "../02_Use_ModelEvaluator/VanDerPol_ModelEvaluator_02.hpp"
#include "TutorialObserver.hpp"

#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_TimeStepControl.hpp"

using namespace std;
using Teuchos::RCP;

/** @file */

/** \page example-07 Example 7: Introduce Observers
 *
 *  This example introduces the observer pattern used to provide application
 *  access during time integration.  The van der Pol model is still provided
 *  through a \ref Thyra::ModelEvaluator, the Forward Euler step is still
 *  performed by \ref Tempus::StepperForwardEuler, timestep-size selection is
 *  still managed by \ref Tempus::TimeStepControl, and the application still
 *  manages the overall time loop and \ref Tempus::SolutionHistory.
 *
 *  The main purpose of this step is to separate application-level access
 *  points, such as output, diagnostics, and custom bookkeeping, from the core
 *  time-integration loop.
 *
 *  This example uses a tutorial-local `TutorialObserver` class, defined in
 *  `TutorialObserver.hpp`.  Its callback names mirror the
 *  \ref Tempus::IntegratorObserver lifecycle, but the arguments are adapted to
 *  the application-managed loop used in this tutorial.  In \ref example-08,
 *  these ideas are connected to \ref Tempus::IntegratorBasic and the
 *  production \ref Tempus::IntegratorObserverBasic interface.
 *
 *  Relative to \ref example-06:
 *  - a tutorial observer class is introduced
 *  - observer callbacks are added at key points in the time-step loop
 *  - step output is moved from the main loop into the observer
 *  - the stepper, timestep control, and solution-history logic remain
 *    essentially unchanged
 *
 *  The central idea behind observers is that applications often need access
 *  to integration events, solution states, and step-level activity without
 *  modifying the time-integration algorithm itself.  Observer callbacks provide
 *  structured access to those events.
 *
 *  This example uses only part of the full observer concept:
 *  - beginning and ending the integration
 *  - beginning and ending each time step
 *  - observing after timestep metadata is selected
 *  - observing before and after the stepper advances the solution
 *  - observing after the step status is checked
 *
 *  <hr>
 *  \par Transition notes
 *  See \ref tempus_tutorial_transition_06_07 for a detailed explanation
 *  of what changed from \ref example-06.
 *
 *  \htmlonly
 *  <div style="text-align:center;">
 *    <a href="example-06.html">← Previous Example</a> |
 *    <a href="tempus_tutorials_overview.html">Tutorial Overview</a> |
 *    <a href="example-08.html">Next Example →</a>
 *  </div>
 *  \endhtmlonly
 */
int main(int argc, char *argv[])
{
  bool verbose = true;
  bool success = false;
  try {
    // Construct ModelEvaluator
    Teuchos::RCP<const Thyra::ModelEvaluator<double> >
      model = Teuchos::rcp(new VanDerPol_ModelEvaluator_02<double>());

    // Setup initial condition SolutionState
    auto solState = Tempus::createSolutionStateX(
                      model->getNominalValues().get_x()->clone_v());
    solState->setIndex   (0);
    solState->setTime    (0.0);
    solState->setTimeStep(0.0);  // By convention, the IC has dt = 0.
    solState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are considered passed.

    // Create SolutionHistory
    auto solHistory = Tempus::createSolutionHistoryState<double>(solState);

    // Create and initialize StepperForwardEuler
    auto stepper = Teuchos::rcp(new Tempus::StepperForwardEuler<double>());
    stepper->setModel(model);
    stepper->initialize();
    stepper->setInitialConditions(solHistory);

    // Create and initialize TimeStepControl
    auto timeStepControl = Teuchos::rcp(new Tempus::TimeStepControl<double>());
    timeStepControl->setFinalTime(2.0);
    timeStepControl->setNumTimeSteps(2000);
    timeStepControl->initialize();

    Tempus::Status integratorStatus = Tempus::Status::WORKING;

    // Create the tutorial observer.
    Tempus_Tutorial::TutorialObserver<double> observer;

    observer.observeStartIntegrator(
      solHistory, stepper, timeStepControl, integratorStatus);

    // Advance the solution to the next timestep.
    while (solHistory->getCurrentState()->getSolutionStatus() == Tempus::Status::PASSED &&
           timeStepControl->timeInRange(solHistory->getCurrentTime()) &&
           timeStepControl->indexInRange(solHistory->getCurrentIndex())) {

      observer.observeStartTimeStep(
        solHistory, stepper, timeStepControl, integratorStatus);

      // Initialize next time step using SolutionHistory
      solHistory->initWorkingState();

      // Let TimeStepControl determine the next time-step metadata.
      timeStepControl->setNextTimeStep(solHistory, integratorStatus);

      observer.observeNextTimeStep(
        solHistory, stepper, timeStepControl, integratorStatus);

      observer.observeBeforeTakeStep(
        solHistory, stepper, timeStepControl, integratorStatus);

      // Take one Forward Euler step through Tempus::StepperForwardEuler
      stepper->takeStep(solHistory);

      observer.observeAfterTakeStep(
        solHistory, stepper, timeStepControl, integratorStatus);

      observer.observeAfterCheckTimeStep(
        solHistory, stepper, timeStepControl, integratorStatus);

      // Promote working state to current state
      solHistory->promoteWorkingState();

      observer.observeEndTimeStep(
        solHistory, stepper, timeStepControl, integratorStatus);
    }

    observer.observeEndIntegrator(
      solHistory, stepper, timeStepControl, integratorStatus);

    // Test for regression.
    auto finalState = solHistory->getCurrentState();
    bool passed = (finalState->getSolutionStatus() == Tempus::Status::PASSED);
    bool regressionPassed = tutorialRegressionTest(finalState);

    if (passed && regressionPassed) success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if(success)
    cout << "\nEnd Result: Test Passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
