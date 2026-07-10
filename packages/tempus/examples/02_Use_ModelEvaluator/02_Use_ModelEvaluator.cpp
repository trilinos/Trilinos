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

#include "VanDerPol_ModelEvaluator_02.hpp"


using namespace std;
using Teuchos::RCP;

/** @file */

/** \page example-02 Example 2: Use ModelEvaluator
 *
 *  This example moves the van der Pol model into a
 *  \ref Thyra::ModelEvaluator. The time integration loop still performs a
 *  hand-written Forward Euler update, but the right-hand-side evaluation and
 *  nominal initial conditions are now obtained through the model interface.
 *
 *  The main purpose of this step is to separate the problem physics from the
 *  time integration algorithm.
 *
 *  Relative to \ref example-01:
 *  - the problem physics are encapsulated in a \ref Thyra::ModelEvaluator,
 *    `VanDerPol_ModelEvaluator_02`
 *  - initial conditions are obtained from `getNominalValues()`
 *  - right-hand-side evaluation is performed through `evalModel()`
 *  - the application begins to separate the model definition from the time
 *    integration logic
 *
 *  \ref Thyra::ModelEvaluator provides a common Trilinos interface between
 *  application models and
 *  <a href="https://www.osti.gov/servlets/purl/1264635">abstract numerical algorithms</a>.
 *  For \ref Tempus, this interface is the standard mechanism for exposing
 *  governing equations to steppers and integrators.
 *
 *  This example focuses only on the subset of \ref Thyra::ModelEvaluator
 *  behavior needed for explicit first-order ODEs.
 *
 *  As in the previous examples, the code prints the evolving solution in a
 *  simple table with columns for the step index, time, \f$x_0\f$, and
 *  \f$x_1\f$. The variable `passed` still indicates whether the hand-written
 *  stepping loop succeeds, while overall tutorial success additionally
 *  requires the final regression comparison to pass.
 *
 *  <hr>
 *  \par Transition notes
 *  See \ref tempus_tutorial_transition_01_02 for a detailed explanation
 *  of what changed from \ref example-01.
 *
 *  \htmlonly
 *  <div style="text-align:center;">
 *    <a href="example-01.html">← Previous Example</a> |
 *    <a href="tempus_tutorials_overview.html">Tutorial Overview</a> |
 *    <a href="example-03.html">Next Example →</a>
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

    // Get initial conditions from ModelEvaluator::getNominalValues.
    int n = 0;
    double time = 0.0;
    bool passed = true;   // ICs are considered passed.
    RCP<Thyra::VectorBase<double> > x_n    =
      model->getNominalValues().get_x()->clone_v();
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

    cout << std::setw(8)  << n
         << std::setw(10) << std::setprecision(3) << time
         << std::setw(12) << std::setprecision(4) << get_ele(*x_n, 0)
         << std::setw(12) << std::setprecision(4) << get_ele(*x_n, 1) << endl;

    // Advance the solution to the next timestep.
    while (passed && time < finalTime && n < nTimeSteps) {

      // Initialize next time step
      RCP<Thyra::VectorBase<double> > x_np1 = x_n->clone_v(); // at time index n+1

      // Set the timestep and time.
      double dt = constDT;
      time = (n+1)*dt;

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
        passed = false;
      } else {
        // Promote to next step (n <- n+1).
        Thyra::V_V(x_n.ptr(), *x_np1);
        n++;
      }

      // Output
      if (n % 100 == 0)
        cout << std::setw(8)  << n
             << std::setw(10) << std::setprecision(3) << time
             << std::setw(12) << std::setprecision(4) << get_ele(*x_n, 0)
             << std::setw(12) << std::setprecision(4) << get_ele(*x_n, 1) << endl;
    }

    // Test for regression.
    bool regressionPassed = tutorialRegressionTest(x_n);

    if (passed && regressionPassed) success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if(success)
    cout << "\nEnd Result: Test Passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
