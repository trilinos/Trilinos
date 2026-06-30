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


using namespace std;
using Teuchos::RCP;

/** @file */

/** \page example-01 Example 1: Utilize Thyra
 *
 *  This example replaces the raw C++ arrays from \ref example-00 with
 *  Thyra vectors and vector spaces.  The mathematical problem, timestepper,
 *  and overall time-integration structure remain essentially unchanged,
 *  The main purpose is to introduce the abstract vector interfaces used
 *  throughout \ref Tempus and other Trilinos packages.
 *
 *  Changes relative to \ref example-00:
 *  - `Teuchos::RCP` is used for memory management.
 *  - state storage moves from raw arrays to `Thyra::VectorBase`,
 *  - a `Thyra::VectorSpaceBase` is introduced to define the state space,
 *  - vector initialization and element access use `Thyra::DetachedVectorView`,
 *  - vector algebra uses Thyra helper routines, `Thyra::V_VpStV`, `Thyra::norm`, `Thyra::V_V`
 *
 *  This example shows how the same application logic can be expressed in
 *  terms of <a href="https://www.osti.gov/servlets/purl/1264635">abstract numerical objects</a> rather than concrete array storage.
 *  That abstraction is a prerequisite for later examples that introduce
 *  `Thyra::ModelEvaluator` and \ref Tempus solution-management objects.
 *
 *  <hr>
 *  \par Transition notes
 *  See \ref tempus_tutorial_transition_00_01 for a detailed explanation
 *  of what changed from \ref example-00.
 *
 *  \htmlonly
 *  <div style="text-align:center;">
 *    <a href="example-00.html">← Previous Example</a> |
 *    <a href="tempus_tutorials_overview.html">Tutorial Overview</a> |
 *    <a href="example-02.html">Next Example →</a>
 *  </div>
 *  \endhtmlonly
 */
int main(int argc, char *argv[])
{
  bool verbose = true;
  bool success = false;
  try {
    // Solution and its time-derivative.
    int vectorLength = 2;  // number state unknowns
    RCP<const Thyra::VectorSpaceBase<double> >  xSpace =
      Thyra::defaultSpmdVectorSpace<double>(vectorLength);

    RCP<Thyra::VectorBase<double> > x_n    = Thyra::createMember(xSpace);
    RCP<Thyra::VectorBase<double> > xDot_n = Thyra::createMember(xSpace);

    // Initial Conditions
    int n = 0;
    double time = 0.0;
    double epsilon = 1.0e-1;
    bool passed = true;   // ICs are considered passed.
    { // Scope to delete DetachedVectorViews
      Thyra::DetachedVectorView<double> x_n_view(*x_n);
      x_n_view[0] = 2.0;
      x_n_view[1] = 0.0;
      Thyra::DetachedVectorView<double> xDot_n_view(*xDot_n);
      xDot_n_view[0] = 0.0;
      xDot_n_view[1] = -2.0/epsilon;
    }

    // Timestep size
    double finalTime = 2.0;
    int nTimeSteps = 2000;
    const double constDT = finalTime/nTimeSteps;

    // Advance the solution to the next timestep.
    cout << n << "  " << time << "  " << get_ele(*(x_n), 0)
                              << "  " << get_ele(*(x_n), 1) << endl;
    while (passed && time < finalTime && n < nTimeSteps) {

      // Initialize next time step
      RCP<Thyra::VectorBase<double> > x_np1 = x_n->clone_v(); // at time index n+1

      // Set the timestep and time.
      double dt = constDT;
      time = (n+1)*dt;

      // Righthand side evaluation and time-derivative at n.
      {
        Thyra::ConstDetachedVectorView<double> x_n_view(*x_n);
        Thyra::DetachedVectorView<double> xDot_n_view(*xDot_n);
        xDot_n_view[0] = x_n_view[1];
        xDot_n_view[1] =
          ((1.0-x_n_view[0]*x_n_view[0])*x_n_view[1]-x_n_view[0])/epsilon;
      }

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
      if ( n%100 == 0 )
        cout << n << "  " << time << "  " << get_ele(*(x_n), 0)
                                  << "  " << get_ele(*(x_n), 1) << endl;
    }

    // Test for regression.
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
      passed = false;
      cout << "FAILED regression constraint!" << endl;
    }
    if (passed) success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if(success)
    cout << "\nEnd Result: Test Passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
