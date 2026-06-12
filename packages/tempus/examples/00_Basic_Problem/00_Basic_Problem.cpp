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

using namespace std;

/** @file */

/** \page example-00 Example 0: Basic Problem
 *
 *  This example provides a minimal application code for integrating the
 *  van der Pol problem using a hand-written Forward Euler time loop.
 *  It does not yet use Tempus, Thyra, or other Trilinos abstractions for
 *  the state or time integration algorithm.
 *
 *  The purpose of this example is to establish the baseline structure used
 *  throughout the tutorial sequence.  Later examples introduce \ref Tempus
 *  and Trilinos capabilities one step at a time while preserving the same
 *  basic problem setup whenever possible.
 *
 *  \section vanderpol van der Pol Problem
 *
 *  The scaled explicit first-order ODE is
 *  \f{eqnarray*}{
 *    \dot{x}_0(t) & = & x_1(t) \\
 *    \dot{x}_1(t) & = & \left[(1-x_0^2)x_1-x_0\right]/\epsilon
 *  \f}
 *  with initial conditions
 *  \f{eqnarray*}{
 *    x_0(0) & = & 2 \\
 *    x_1(0) & = & 0.
 *  \f}
 *
 *  For the model definition and additional details, 
 *  see \ref Tempus_Test::VanDerPolModel
 *  and \ref tempus_vanderpol_model "van der Pol Model".
 *
 *  \section example-00-outline Outline
 *
 *  This example demonstrates the basic structure of an application-level
 *  time integration loop:
 *  - declare the state and its time derivative
 *  - set the initial conditions
 *  - advance the solution from the initial time to the final time
 *    - select the timestep size
 *    - evaluate the right-hand side of the governing equations
 *    - apply the Forward Euler update
 *    - check a simple pass/fail criterion
 *    - accept the step and promote the solution
 *
 *  The regression check at the end of the run is secondary to the tutorial
 *  discussion and can be ignored on a first reading.
 *
 *  In this version:
 *  - the solution is stored in raw C++ arrays,
 *  - the right-hand side is evaluated directly in the application code,
 *  - timestep management is handled manually,
 *  - solution status is tracked with local scalar variables.
 *
 *  The next example replaces raw arrays with \ref Thyra vectors while
 *  preserving the same overall algorithmic structure.
 *
 *  \htmlonly
 *  <div style="text-align:center;">
 *    <a href="tempus_tutorials_overview.html">Tutorial Overview</a> |
 *    <a href="example-01.html">Next Example →</a>
 *  </div>
 *  \endhtmlonly
 */
int main(int argc, char *argv[])
{
  bool verbose = true;
  bool success = false;
  try {
    // Solution and its time-derivative.
    double x_n[2];      // at time index n
    double xDot_n[2];   // at time index n

    // Initial Conditions
    int n = 0;
    double time = 0.0;
    double epsilon = 1.0e-1;
    bool passed = true;   // ICs are considered passed.
    x_n   [0] = 2.0;
    x_n   [1] = 0.0;
    xDot_n[0] = 0.0;
    xDot_n[1] = -2.0/epsilon;

    // Timestep size
    double finalTime = 2.0;
    int nTimeSteps = 2001;
    const double constDT = finalTime/(nTimeSteps-1);

    // Advance the solution to the next timestep.
    cout << n << "  " << time << "  " << x_n[0] << "  " << x_n[1] << endl;
    while (passed && time < finalTime && n < nTimeSteps) {

      // Initialize next time step
      double x_np1[2];    // at time index n+1
      x_np1[0] = x_n[0];
      x_np1[1] = x_n[1];

      // Set the timestep and time.
      double dt = constDT;
      time = (n+1)*dt;

      // Righthand side evaluation and time-derivative at n.
      xDot_n[0] = x_n[1];
      xDot_n[1] = ((1.0 - x_n[0]*x_n[0])*x_n[1] - x_n[0])/epsilon;

      // Take the timestep - Forward Euler
      x_np1[0] = x_n[0] + dt*xDot_n[0];
      x_np1[1] = x_n[1] + dt*xDot_n[1];

      // Test if solution has passed.
      if ( std::isnan(x_n[0]) || std::isnan(x_n[1]) ) {
        passed = false;
      } else {
        // Promote to next step (n <- n+1).
        x_n[0] = x_np1[0];
        x_n[1] = x_np1[1];
        n++;
      }

      // Output
      if ( n%100 == 0 )
        cout << n << "  " << time << "  " << x_n[0] << "  " << x_n[1] << endl;

    }

    // Test for regression.
    double x_regress[2];      // Regression results for nTimeSteps = 2001
    x_regress[0] = -1.59496108218721311;
    x_regress[1] =  0.96359412806611255;
    double x_L2norm_error = 0.0;
    double x_L2norm_regress = 0.0;
    for (int i=0; i < 2; i++) {
      x_L2norm_error += (x_n[i]-x_regress[i])*(x_n[i]-x_regress[i]);
      x_L2norm_regress += x_regress[1]*x_regress[1];
    }
    x_L2norm_error   = sqrt(x_L2norm_error  );
    x_L2norm_regress = sqrt(x_L2norm_regress);
    cout << "Relative L2 Norm of the error (regression) = "
         << x_L2norm_error/x_L2norm_regress << endl;
    if ( x_L2norm_error > 1.0e-08*x_L2norm_regress) {
      passed = false;
      cout << "FAILED regression constraint!" << endl;
    }

    double x_best[2];         // Best resolution with nTimeSteps = 2000000001
    x_best[0]    = -1.58184083624543947;
    x_best[1]    =  0.97844890081968072;
    x_L2norm_error = 0.0;
    double x_L2norm_best = 0.0;
    for (int i=0; i < 2; i++) {
      x_L2norm_error += (x_n[i]-x_best[i])*(x_n[i]-x_best[i]);
      x_L2norm_best += x_best[1]*x_best[1];
    }
    x_L2norm_error = sqrt(x_L2norm_error);
    x_L2norm_best  = sqrt(x_L2norm_best );
    cout << "Relative L2 Norm of the error (best)       = "
         << x_L2norm_error/x_L2norm_best << endl;
    if ( x_L2norm_error > 0.02*x_L2norm_best) {
      passed = false;
      cout << "FAILED best constraint!" << endl;
    }
    if (passed) success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if(success)
    cout << "\nEnd Result: Test Passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
