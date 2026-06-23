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


using namespace std;
using Teuchos::RCP;

/** @file */

/** \page example-04 Example 4: Add SolutionHistory
 *
 *  This example introduces \ref Tempus::SolutionHistory, which manages
 *  multiple \ref Tempus::SolutionState objects during time integration.
 *  The van der Pol model is still provided through a
 *  \ref Thyra::ModelEvaluator, and the Forward Euler update is still written
 *  explicitly in the application code, but the evolving solution is now
 *  organized through a history object rather than a single state.
 *
 *  The main purpose of this step is to move from managing one solution state
 *  at a time to managing a sequence of states through a core \ref Tempus
 *  container.
 *
 *  Relative to \ref example-03:
 *  - the primary state is now stored in a \ref Tempus::SolutionHistory
 *  - the current and working states are accessed through the history object
 *  - the history object initializes the working state for each new step
 *  - successful steps are accepted by promoting the working state
 *  - the example begins to follow the state-history management pattern used
 *    by \ref Tempus steppers and integrators
 *
 *  The central idea behind \ref Tempus::SolutionHistory is that time
 *  integration algorithms often need access to more than one solution state.
 *  A history object provides a structured way to manage the current state,
 *  tentative working states, and previously accepted states in support of
 *  stepping algorithms, interpolation, restart, and error recovery.
 *
 *  This example uses only part of the full \ref Tempus::SolutionHistory
 *  capability:
 *  - storing the current and working solution states
 *  - initializing a working state through `initWorkingState()`
 *  - accepting a step through `promoteWorkingState()`
 *
 *  For additional details, see \ref SolutionHistory_Description.
 *
 *  <hr>
 *  \par Transition notes
 *  See \ref tempus_tutorial_transition_03_04 for a detailed explanation
 *  of what changed from \ref example-03.
 *
 *  \htmlonly
 *  <div style="text-align:center;">
 *    <a href="example-03.html">← Previous Example</a> |
 *    <a href="tempus_tutorials_overview.html">Tutorial Overview</a> |
 *    <a href="example-05.html">Next Example →</a>
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
    RCP<Thyra::VectorBase<double> > xDot_n =
      model->getNominalValues().get_x_dot()->clone_v();

    // Create SolutionHistory
    auto solHistory = Tempus::createSolutionHistoryState<double>(solState);

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
      auto currentState = solHistory->getCurrentState();
      auto workingState = solHistory->getWorkingState();
      RCP<Thyra::VectorBase<double> > x_n    = currentState->getX();
      RCP<Thyra::VectorBase<double> > x_np1  = workingState->getX();

      // Set the timestep and time for the working solution, i.e., n+1.
      int index = workingState->getIndex();  // Already incremented by initWorkingState()
      double dt = constDT;
      double time = index*dt;
      workingState->setTime(time);
      workingState->setTimeStep(dt);
      
      // For explicit ODE formulation, xDot = f(x, t),
      // xDot is part of the outArgs.
      auto inArgs  = model->createInArgs();
      auto outArgs = model->createOutArgs();
      inArgs.set_t(time);
      inArgs.set_x(x_n);
      inArgs.set_x_dot(Teuchos::null);
      outArgs.set_f(xDot_n);

      // Righthand side evaluation and time-derivative at n.
      model->evalModel(inArgs, outArgs);

      // Take the timestep - Forward Euler
      Thyra::V_VpStV(x_np1.ptr(), *x_n, dt, *xDot_n);

      // Test if solution has passed.
      if ( std::isnan(Thyra::norm(*x_np1)) ) {
        workingState->setSolutionStatus(Tempus::Status::FAILED);
      } else {
        workingState->setSolutionStatus(Tempus::Status::PASSED);
      }
      // Promote working state to current state
      solHistory->promoteWorkingState();

      // Output
      if (solHistory->getCurrentState()->getIndex() % 100 == 0) {
        currentState = solHistory->getCurrentState();
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
