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

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"

#include "../02_Use_ModelEvaluator/VanDerPol_ModelEvaluator_02.hpp"

#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_StepperForwardEuler.hpp"


using namespace std;
using Teuchos::RCP;

/** @file */

/** \page example-05 Example 5: Introduce Stepper
 *
 *  This example introduces \ref Tempus::Stepper through
 *  \ref Tempus::StepperForwardEuler.  The van der Pol model is still provided
 *  through a \ref Thyra::ModelEvaluator, and the application still manages the
 *  overall time loop and \ref Tempus::SolutionHistory.  However, the Forward
 *  Euler step itself is now delegated to a Tempus stepper object rather than
 *  being written explicitly in the application code.
 *
 *  The main purpose of this step is to separate the stepping algorithm from
 *  the surrounding application logic.
 *
 *  Relative to \ref example-04:
 *  - the Forward Euler update is performed by a \ref Tempus::Stepper
 *  - the specific stepper used is \ref Tempus::StepperForwardEuler
 *  - the application still manages the time loop and solution history
 *  - the stepper advances the working state using the model and timestep
 *    information already stored in the history
 *
 *  The central idea behind \ref Tempus::Stepper is that a time integration
 *  algorithm should be encapsulated in its own object.  This allows the
 *  application to reuse the same surrounding logic while changing only the
 *  stepping algorithm.
 *
 *  This example uses only part of the full \ref Tempus::Stepper capability:
 *  - construction and initialization of a \ref Tempus::StepperForwardEuler
 *  - attaching a \ref Thyra::ModelEvaluator to the stepper
 *  - setting initial conditions for the stepper
 *  - advancing one step at a time through the stepper
 *  - continued use of \ref Tempus::SolutionHistory under application control
 *
 *  <hr>
 *  \par Transition notes
 *  See \ref tempus_tutorial_transition_04_05 for a detailed explanation
 *  of what changed from \ref example-04.
 *
 *  \htmlonly
 *  <div style="text-align:center;">
 *    <a href="example-04.html">← Previous Example</a> |
 *    <a href="tempus_tutorials_overview.html">Tutorial Overview</a> |
 *    <a href="example-06.html">Next Example →</a>
 *  </div>
 *  \endhtmlonly
 */
int main(int  /*argc*/, char * /*argv*/[])
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

    // Timestep size
    double finalTime = 2.0;
    int nTimeSteps = 2000;
    const double constDT = finalTime/nTimeSteps;

    // Advance the solution to the next timestep.
    while (solHistory->getCurrentState()->getSolutionStatus() == Tempus::Status::PASSED &&
           solHistory->getCurrentState()->getTime() < finalTime &&
           solHistory->getCurrentState()->getIndex() < nTimeSteps) {

      // Initialize next time step using SolutionHistory
      solHistory->initWorkingState();
      auto workingState = solHistory->getWorkingState();

      // Set the timestep and time for the working solution, i.e., n+1.
      int index = workingState->getIndex();  // Already incremented by initWorkingState()
      double dt = constDT;
      double time = index*dt;
      workingState->setTime(time);
      workingState->setTimeStep(dt);

      // Take one Forward Euler step through Tempus::StepperForwardEuler
      stepper->takeStep(solHistory);

      // Promote working state to current state
      solHistory->promoteWorkingState();

      // Output
      if (solHistory->getCurrentState()->getIndex() % 100 == 0) {
        auto currentState = solHistory->getCurrentState();
        cout << currentState->getIndex() << "  " << currentState->getTime()
             << "  " << Thyra::get_ele(*(currentState->getX()), 0)
             << "  " << Thyra::get_ele(*(currentState->getX()), 1) << endl;
      }
    }

    // Test for regression.
    auto finalState = solHistory->getCurrentState();
    RCP<Thyra::VectorBase<double> > x_n = finalState->getX();
    RCP<Thyra::VectorBase<double> > x_regress = x_n->clone_v();
    {
      Thyra::DetachedVectorView<double> x_regress_view(*x_regress);
      x_regress_view[0] = -1.59496108218721311;
      x_regress_view[1] =  0.96359412806611255;
    }

    RCP<Thyra::VectorBase<double> > x_error = x_n->clone_v();
    Thyra::V_VmV(x_error.ptr(), *x_n, *x_regress);
    double x_L2norm_error   = Thyra::norm_2(*x_error  );
    double x_L2norm_regress = Thyra::norm_2(*x_regress);

    cout << "Relative L2 Norm of the error (regression) = "
         << x_L2norm_error/x_L2norm_regress << endl;
    if ( x_L2norm_error > 1.0e-08*x_L2norm_regress) {
      finalState->setSolutionStatus(Tempus::Status::FAILED);
      cout << "FAILED regression constraint!" << endl;
    }
    if (finalState->getSolutionStatus() == Tempus::Status::PASSED) success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if(success)
    cout << "\nEnd Result: Test Passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
