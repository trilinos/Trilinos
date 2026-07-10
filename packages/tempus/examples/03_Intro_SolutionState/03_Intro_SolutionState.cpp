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


using namespace std;
using Teuchos::RCP;

/** @file */

/** \page example-03 Example 3: Introduce SolutionState
 *
 *  This example introduces \ref Tempus::SolutionState, which stores the
 *  solution and selected metadata at a particular time. The model is still
 *  provided through a \ref Thyra::ModelEvaluator, and the Forward Euler update
 *  is still written explicitly in the application code, but the evolving state
 *  of the integration is now organized using a core \ref Tempus data
 *  structure.
 *
 *  The main purpose of this step is to move from ad hoc scalar bookkeeping to
 *  a structured representation of the solution and its integration metadata.
 *
 *  Relative to \ref example-02:
 *  - the primary solution is stored in a \ref Tempus::SolutionState
 *  - metadata such as index, time, timestep size, and status are tracked
 *    through that object
 *  - pass/fail logic is expressed through \ref Tempus::Status
 *  - the example begins to follow \ref Tempus conventions for managing
 *    solution state during time integration
 *
 *  The central idea behind \ref Tempus::SolutionState is that it should
 *  contain the information needed to represent a solution at a specific time
 *  in a form suitable for restart, checkpointing, diagnostics, and recovery
 *  from failed timesteps.
 *
 *  This example uses only part of the full \ref Tempus::SolutionState
 *  capability:
 *  - the solution vector \f$x(t)\f$
 *  - selected metadata from \ref Tempus::SolutionStateMetaData
 *
 *  For additional details, see \ref SolutionState_Description.
 *
 *  The example continues to print the evolving solution in a simple table with
 *  columns for step index, time, \f$x_0\f$, and \f$x_1\f$. In this example,
 *  the local variable `passed` used at the end reflects whether the final
 *  \ref Tempus::SolutionState status is `PASSED`. Overall tutorial success
 *  additionally requires that the final solution satisfy the regression check.
 *
 *  <hr>
 *  \par Transition notes
 *  See \ref tempus_tutorial_transition_02_03 for a detailed explanation
 *  of what changed from \ref example-02.
 *
 *  \htmlonly
 *  <div style="text-align:center;">
 *    <a href="example-02.html">← Previous Example</a> |
 *    <a href="tempus_tutorials_overview.html">Tutorial Overview</a> |
 *    <a href="example-04.html">Next Example →</a>
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

    // Setup initial condition SolutionState --------------------
    auto solState = Tempus::createSolutionStateX(
                      model->getNominalValues().get_x()->clone_v());
    solState->setIndex   (0);
    solState->setTime    (0.0);
    solState->setTimeStep(0.0);  // By convention, the IC has dt=0.
    solState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are considered passed.
    RCP<Thyra::VectorBase<double> > xDot_n =
      model->getNominalValues().get_x_dot()->clone_v();


    // Timestep size
    double finalTime = 2.0;
    int nTimeSteps = 2000;
    const double constDT = finalTime/nTimeSteps;

    // Output
    cout << std::fixed;
    cout << std::setw(8)  << "index"
         << std::setw(10) << "time"
         << std::setw(12) << "x_0"
         << std::setw(12) << "x_1" << endl;

    cout << std::setw(8)  << solState->getIndex()
         << std::setw(10) << std::setprecision(3) << solState->getTime()
         << std::setw(12) << std::setprecision(4) << get_ele(*(solState->getX()), 0)
         << std::setw(12) << std::setprecision(4) << get_ele(*(solState->getX()), 1)
         << endl;

    // Advance the solution to the next timestep.
    while (solState->getSolutionStatus() == Tempus::Status::PASSED &&
           solState->getTime() < finalTime &&
           solState->getIndex() < nTimeSteps) {

      // Initialize next time step
      RCP<Thyra::VectorBase<double> > x_n    = solState->getX();
      RCP<Thyra::VectorBase<double> > x_np1  = solState->getX()->clone_v(); // at time index n+1
      solState->setSolutionStatus(Tempus::Status::WORKING);

      // Set the timestep and time for the working solution i.e., n+1.
      int index = solState->getIndex()+1;
      double dt = constDT;
      double time = index*dt;

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
        solState->setSolutionStatus(Tempus::Status::FAILED);
      } else {
        // Promote to next step (n <- n+1).
        Thyra::V_V(x_n.ptr(), *x_np1);
        solState->setIndex   (index);
        solState->setTime    (time);
        solState->setTimeStep(constDT);
        solState->setSolutionStatus(Tempus::Status::PASSED);
      }

      // Output
      if (solState->getIndex() % 100 == 0)
        cout << std::setw(8)  << solState->getIndex()
             << std::setw(10) << std::setprecision(3) << solState->getTime()
             << std::setw(12) << std::setprecision(4) << get_ele(*(solState->getX()), 0)
             << std::setw(12) << std::setprecision(4) << get_ele(*(solState->getX()), 1)
             << endl;
    }

    // Test for regression.
    bool passed = (solState->getSolutionStatus() == Tempus::Status::PASSED);
    bool regressionPassed = tutorialRegressionTest(solState);

    if (passed && regressionPassed) success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if(success)
    cout << "\nEnd Result: Test Passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
