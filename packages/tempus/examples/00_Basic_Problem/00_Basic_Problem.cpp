// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "Teuchos_StandardCatchMacros.hpp"

using namespace std;

/** @file */

/** \brief Description:
 *
 *  \section example-00 Example 0: Basic Problem
 *
 *  This problem is an example of a typical application code before
 *  utilizing Tempus, thus it has no Tempus or other Trilinos
 *  functionality at this point.  It is the starting point to help
 *  illustrate the various stages an application can do to take
 *  advantage of Tempus's capabilities (additionally see the other
 *  Tempus examples).
 *
 *  We have kept the problem itself fairly simple, Van Der Pol problem,
 *  with a simple time stepper, Forward Euler, to allow us to focus
 *  on the time integration structure and Tempus features as we add them.
 *  With this, the implementation is not the most compact or efficient
 *  for this simple problem and stepper, but it does have features
 *  one would want in a more capable time integrator.
 *
 *  \subsection vanderpol van der Pol Problem
 *
 *  van der Pol problem is a nonlinear electrical circuit, and is a
 *  canonical equation of a nonlinear oscillator (Hairer, Norsett,
 *  and Wanner, pp. 111-115, and Hairer and Wanner, pp. 4-5) for an
 *  electrical circuit.  The explicit ODE form, \f$ \dot{x} = f(x,t) \f$,
 *  of the scaled problem is
 *  \f{eqnarray*}{
 *    \dot{x}_0(t) & = & f_0 = x_1(t) \\
 *    \dot{x}_1(t) & = & f_1 = [(1-x_0^2)x_1-x_0]/\epsilon
 *  \f}
 *  where the initial conditions are
 *  \f{eqnarray*}{
 *    x_0(t_0=0) & = & 2 \\
 *    x_1(t_0=0) & = & 0
 *  \f}
 *  and the initial time derivatives are
 *  \f{eqnarray*}{
 *    \dot{x}_0(t_0=0) & = & x_1(t_0=0) = 0 \\
 *    \dot{x}_1(t_0=0) & = & [(1-x_0^2)x_1-x_0]/\epsilon = -2/\epsilon
 *  \f}
 *
 *  \section basicCodeStructure Basic Code Structure
 *
 *  We will review this program's code structure, 00_Basic_Problem.cpp,
 *  and go over its basic details.
 *
 *  \subsection solutionIC Solution and Initial Conditions
 *
 *  One of the first things an application does is to define the
 *  solution and set its initial values.
 *  @code
 *    // Solution and its time-derivative.
 *    double x_n[2];      // at time index n
 *    double xDot_n[2];   // at time index n
 *  @endcode
 *  Here we are using an array of two doubles for the solution and its
 *  time derivative, and set their initial values (see \ref vanderpol).
 *  @code
 *    // Initial Conditions
 *    double time = 0.0;
 *    double epsilon = 1.0e-1;
 *    x_n   [0] = 2.0;
 *    x_n   [1] = 0.0;
 *    xDot_n[0] = 0.0;
 *    xDot_n[1] = -2.0/epsilon;
 *  @endcode
 *  Note that xDot_n is the time derivative, \f$ \dot{x} \f$, but is also
 *  the right-hand side evaluation, \f$ f(x,t) \f$.
 *
 *  \subsection timeStepParameters Timestep Parameters
 *
 *  The next piece of code sets parameters associated with
 *  time-integration, e.g., length of time integration, number of timesteps
 *  and the timestep size.
 *  @code
 *    // Timestep size
 *    double finalTime = 2.0;
 *    int nTimeSteps = 2001;
 *    const double constDT = finalTime/(nTimeSteps-1);
 *  @endcode
 *
 *  \subsection timeLoop Time-Integration Loop
 *
 *  The time-integration loop simply advances the solution until the
 *  solution is not "passing", the time has reached the final time, or
 *  the number timesteps is reached.
 *  @code
 *    // Advance the solution to the next timestep.
 *    int n = 0;
 *    bool passing = true;
 *    cout << n << "  " << time << "  " << x_n[0] << "  " << x_n[1] << endl;
 *    while (passing && time < finalTime && n < nTimeSteps) {
 *
 *      ...
 *
 *    }
 *  @endcode
 *
 *  \subsubsection initNextTimeStep Initialize the Next Timestep
 *
 *  Within the time loop, we initialize the next timestep, \f$ x^{n+1} \f$,
 *  to the previous timestep, \f$ x^n \f$.  Although this is not needed
 *  in this simple case, if one wants to recover from a failed
 *  timestep, the previous timestep, \f$ x^n \f$, needs to be perserved
 *  until the next timestep is accepted.
 *  @code
 *    // Initialize next time step
 *    double x_np1[2];    // at time index n+1
 *    x_np1[0] = x_n[0];
 *    x_np1[1] = x_n[1];
 *  @endcode
 *
 *  \subsubsection setTimeStepTime Set the Timestep and Time
 *
 *  Although this a constant timestep example, for variable timestepping,
 *  the timestep size needs to be determined and the related time for the
 *  next timestep set.
 *  @code
 *    // Set the timestep and time.
 *    double dt = constDT;
 *    time = (n+1)*dt;
 *  @endcode
 *
 *  \subsubsection evalModel Evaluate Righthand Side Function
 *
 *  This is a simple evaluation of the van der Pol problem and
 *  setting it to \f$ \dot{x} \f$.
 *  @code
 *    // Righthand side evaluation and time-derivative at n.
 *    xDot_n[0] = x_n[1];
 *    xDot_n[1] = ((1.0 - x_n[0]*x_n[0])*x_n[1] - x_n[0])/epsilon;
 *  @endcode
 *
 *  \subsubsection takeStep Take the Timestep
 *
 *  For Forward Euler, the timestep is a simple daxpy.
 *  @code
 *    // Take the timestep - Forward Euler
 *    x_np1[0] = x_n[0] + dt*xDot_n[0];
 *    x_np1[1] = x_n[1] + dt*xDot_n[1];
 *  @endcode
 *
 *  \subsubsection acceptStep Test If Solution is Passing.
 *
 *  Here we have a simple test if the solution is passing, i.e.,
 *  the solution does not have a NAN.  Otherwise "promote" the solution
 *  to the next timestep and increment the timestep index.
 *  @code
 *    // Test if solution is passing.
 *    if ( std::isnan(x_n[0]) || std::isnan(x_n[1]) ) {
 *      passing = false;
 *    } else {
 *      // Promote to next step (n -> n+1).
 *      n++;
 *      x_n[0] = x_np1[0];
 *      x_n[1] = x_np1[1];
 *    }
 *  @endcode
 *
 *  \subsubsection observe Output Solution
 *
 *  Of course, we would like to know that progress of the solution
 *  during the simulation.
 *  @code
 *    // Output
 *    if ( n%100 == 0 )
 *      cout << n << "  " << time << "  " << x_n[0] << "  " << x_n[1] << endl;
 *  @endcode
 *
 *  \subsection aftertimeLoop Following the Time-Integration Loop
 *
 *  Following the time loop is some basic regression testing, but
 *  serves to show tasks that one may want to complete after
 *  completion of time integration.
 *
 *  \section example-00-summary Example 0: Basic Problem Summary
 *
 *  So this outlines the basic time integration of the van der Pol
 *  problem with a few features.  We will use this example
 *  and introduce Tempus and Trilinos capabilities in the following
 *  examples.
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
    double time = 0.0;
    double epsilon = 1.0e-1;
    x_n   [0] = 2.0;
    x_n   [1] = 0.0;
    xDot_n[0] = 0.0;
    xDot_n[1] = -2.0/epsilon;

    // Timestep size
    double finalTime = 2.0;
    int nTimeSteps = 2001;
    const double constDT = finalTime/(nTimeSteps-1);

    // Advance the solution to the next timestep.
    int n = 0;
    bool passing = true;
    cout << n << "  " << time << "  " << x_n[0] << "  " << x_n[1] << endl;
    while (passing && time < finalTime && n < nTimeSteps) {

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

      // Test if solution is passing.
      if ( std::isnan(x_n[0]) || std::isnan(x_n[1]) ) {
        passing = false;
      } else {
        // Promote to next step (n -> n+1).
        n++;
        x_n[0] = x_np1[0];
        x_n[1] = x_np1[1];
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
      passing = false;
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
      passing = false;
      cout << "FAILED best constraint!" << endl;
    }
    if (passing) success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if(success)
    cout << "\nEnd Result: Test Passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
