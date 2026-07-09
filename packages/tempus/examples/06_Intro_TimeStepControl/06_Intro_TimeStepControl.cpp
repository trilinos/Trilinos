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

#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_TimeStepControl.hpp"

using namespace std;
using Teuchos::RCP;

/** @file */

/** \page example-06 Example 6: Introduce TimeStepControl
 *
 *  This example introduces \ref Tempus::TimeStepControl, which manages
 *  timestep-size selection independently of the stepper and application time
 *  loop. The van der Pol model is still provided through a
 *  \ref Thyra::ModelEvaluator, the Forward Euler step is still performed by
 *  \ref Tempus::StepperForwardEuler, and the application still manages the
 *  overall time loop and \ref Tempus::SolutionHistory. However, timestep-size
 *  selection is now delegated to a \ref Tempus::TimeStepControl object rather
 *  than being hardcoded directly in the application.
 *
 *  The main purpose of this step is to separate timestep-selection policy from
 *  both the stepping algorithm and the surrounding application logic.
 *  It should be noted that a \ref Tempus::Stepper may suggest a future
 *  timestep size, but the final timestep-size selection occurs in
 *  \ref Tempus::TimeStepControl.
 *
 *  Relative to \ref example-05:
 *  - timestep control is moved into \ref Tempus::TimeStepControl
 *  - the application initializes a dedicated timestep-control object
 *  - the time loop uses the control object to determine whether the current
 *    time and step index remain in range
 *  - timestep metadata for each new working state is assigned through the
 *    control object rather than directly in the application
 *  - the stepper and solution-history logic remain essentially unchanged
 *
 *  The central idea behind \ref Tempus::TimeStepControl is that timestep size
 *  should be managed by a dedicated object that can enforce integration
 *  limits, policies, and strategies independently of the stepper itself.
 *
 *  This example uses only part of the full \ref Tempus::TimeStepControl
 *  capability:
 *  - setting the final integration time
 *  - setting the number of time steps for constant-timestep control
 *  - using the control object to determine timestep metadata for each step
 *  - using the control object to test whether time and step index remain in range
 *
 *  The example continues to print the evolving solution in a simple table with
 *  columns for step index, time, \f$x_0\f$, and \f$x_1\f$. As in the previous
 *  examples, overall tutorial success should require both a successful
 *  integration status and a passing regression comparison of the final
 *  solution.
 *
 *  <hr>
 *  \par Transition notes
 *  See \ref tempus_tutorial_transition_05_06 for a detailed explanation
 *  of what changed from \ref example-05.
 *
 *  \htmlonly
 *  <div style="text-align:center;">
 *    <a href="example-05.html">← Previous Example</a> |
 *    <a href="tempus_tutorials_overview.html">Tutorial Overview</a> |
 *    <a href="example-07.html">Next Example →</a>
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

    // Output
    cout << std::fixed;
    cout << std::setw(8)  << "index"
         << std::setw(10) << "time"
         << std::setw(12) << "x_0"
         << std::setw(12) << "x_1" << endl;

    auto currentState = solHistory->getCurrentState();
    cout << std::setw(8)  << currentState->getIndex()
         << std::setw(10) << std::setprecision(3) << currentState->getTime()
         << std::setw(12) << std::setprecision(4) << Thyra::get_ele(*(currentState->getX()), 0)
         << std::setw(12) << std::setprecision(4) << Thyra::get_ele(*(currentState->getX()), 1)
         << endl;

    // Advance the solution to the next timestep.
    while (solHistory->getCurrentState()->getSolutionStatus() == Tempus::Status::PASSED &&
           timeStepControl->timeInRange(solHistory->getCurrentTime()) &&
           timeStepControl->indexInRange(solHistory->getCurrentIndex())) {

      // Initialize next time step using SolutionHistory
      solHistory->initWorkingState();

      // Let TimeStepControl determine the next time-step metadata.
      timeStepControl->setNextTimeStep(solHistory, integratorStatus);

      // Take one Forward Euler step through Tempus::StepperForwardEuler
      stepper->takeStep(solHistory);

      // Promote working state to current state
      solHistory->promoteWorkingState();

      // Output
      if (solHistory->getCurrentState()->getIndex() % 100 == 0) {
        currentState = solHistory->getCurrentState();
        cout << std::setw(8)  << currentState->getIndex()
             << std::setw(10) << std::setprecision(3) << currentState->getTime()
             << std::setw(12) << std::setprecision(4) << Thyra::get_ele(*(currentState->getX()), 0)
             << std::setw(12) << std::setprecision(4) << Thyra::get_ele(*(currentState->getX()), 1)
             << endl;
      }
    }

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
