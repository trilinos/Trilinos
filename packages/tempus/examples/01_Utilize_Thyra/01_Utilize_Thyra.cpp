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

/** \brief Description:
 *
 *  \section example-01 Example 1: Utilize Thyra
 *
 *  This problem takes \ref 00_Basic_Problem.cpp "Basic Problem"
 *  and utilizes Thyra vectors, replacing the C++ double arrays,
 *  01_Utilize_Thyra.cpp.  Tempus uses Thyra for its Abstract Numerical
 *  Algorithms (ANAs), which are the mathematical concepts of
 *  vectors, vector spaces, and linear operators. All other ANA
 *  interfaces and support software are built on these fundamental
 *  operator/vector interfaces (including ModelEvaluators).
 *
 *  In the following table, code snippets from the \ref 00_Basic_Problem.cpp
 *  "Basic Problem" tutorial are replaced with snippets using Thyra to create
 *  01_Utilize_Thyra.cpp for the \ref 01_Utilize_Thyra.cpp "Utilize Thyra"
 *  tutorial.  This is similar to performing a diff between
 *  00_Basic_Problem.cpp and 01_Utilize_Thyra.cpp, but the first column
 *  provides comments related to the changes.  You may want to to do a
 *  diff (e.g., vimdiff or tkdiff) to see these changes within main
 *  (e.g., vimdiff 00_Basic_Problem/00_Basic_Problem.cpp 01_Utilize_Thyra/01_Utilize_Thyra.cpp).
 *
 *  <table>
 *    <tr> <th> Comments <th> Original "Basic Problem" Code Snippet
 *                        <th> Replacement "Utilize Thyra" Code Snippet
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Setup Thyra vectors.</b> We first need to replace the C++
 *      double arrays with a vector space to construct the Thyra::Vector.
 *      We additionally have introduced the use of Teuchos
 *      Reference-Counted Pointers (Teuchos:RCP), which are Trilinos's
 *      smart pointers.  Details on RCP can be found at
 *      https://www.osti.gov/servlets/purl/919177.
 *    <td>
 *      @code
 *        // Solution and its time-derivative.
 *        double x_n[2];      // at time index n
 *        double xDot_n[2];   // at time index n
 *      @endcode
 *    <td>
 *      @code
 *        // Solution and its time-derivative.
 *        int vectorLength = 2;  // number state unknowns
 *        RCP<const Thyra::VectorSpaceBase<double> >  xSpace =
 *          Thyra::defaultSpmdVectorSpace<double>(vectorLength);
 *
 *        RCP<Thyra::VectorBase<double> > x_n    = Thyra::createMember(xSpace);
 *        RCP<Thyra::VectorBase<double> > xDot_n = Thyra::createMember(xSpace);
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Initialize Thyra vectors.</b>
 *      The initialization can be achieved with the
 *      Thyra::DetachedVectorView's.  The scoping ensures they are deleted
 *      after use.
 *    <td>
 *      @code
 *        // Initial Conditions
 *        double time = 0.0;
 *        double epsilon = 1.0e-1;
 *        x_n   [0] = 2.0;
 *        x_n   [1] = 0.0;
 *        xDot_n[0] = 0.0;
 *        xDot_n[1] = -2.0/epsilon;
 *      @endcode
 *    <td>
 *      @code
 *        // Initial Conditions
 *        double time = 0.0;
 *        double epsilon = 1.0e-1;
 *        { // Scope to delete DetachedVectorViews
 *          Thyra::DetachedVectorView<double> x_n_view(*x_n);
 *          x_n_view[0] = 2.0;
 *          x_n_view[1] = 0.0;
 *          Thyra::DetachedVectorView<double> xDot_n_view(*xDot_n);
 *          xDot_n_view[0] = 0.0;
 *          xDot_n_view[1] = -2.0/epsilon;
 *        }
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Access Thyra vectors.</b>
 *      Elements of the Thyra::Vector can be quickly accessed with get_ele.
 *    <td>
 *      @code
 *        cout << n << "  " << time << "  " << x_n[0] << "  " << x_n[1] << endl;
 *      @endcode
 *    <td>
 *      @code
 *        cout << n << "  " << time << "  " << get_ele(*(x_n), 0)
 *                                  << "  " << get_ele(*(x_n), 1) << endl;
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Clone Thyra vectors.</b>
 *      To initialize the solution at the next time step, we can simply
 *      clone the current timestep.
 *    <td>
 *      @code
 *        // Initialize next time step
 *        double x_np1[2];    // at time index n+1
 *        x_np1[0] = x_n[0];
 *        x_np1[1] = x_n[1];
 *      @endcode
 *    <td>
 *      @code
 *        // Initialize next time step
 *        RCP<Thyra::VectorBase<double> > x_np1 = x_n->clone_v(); // at time index n+1
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Accessing elements of Thyra vectors.</b>
 *      The model evaluation is achieved through Thyra::DetachedVectorView.
 *    <td>
 *      @code
 *        // Righthand side evaluation and time-derivative at n.
 *        xDot_n[0] = x_n[1];
 *        xDot_n[1] = ((1.0 - x_n[0]*x_n[0])*x_n[1] - x_n[0])/epsilon;
 *      @endcode
 *    <td>
 *      @code
 *        // Righthand side evaluation and time-derivative at n.
 *        {
 *          Thyra::ConstDetachedVectorView<double> x_n_view(*x_n);
 *          Thyra::DetachedVectorView<double> xDot_n_view(*xDot_n);
 *          xDot_n_view[0] = x_n_view[1];
 *          xDot_n_view[1] =
 *            ((1.0-x_n_view[0]*x_n_view[0])*x_n_view[1]-x_n_view[0])/epsilon;
 *        }
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Compute with Thyra vectors.</b>
 *      The Forward Euler time stepping is achieved by using
 *      Thyra::V_VpStV, which performs an axpy.
 *    <td>
 *      @code
 *        // Take the timestep - Forward Euler
 *        x_np1[0] = x_n[0] + dt*xDot_n[0];
 *        x_np1[1] = x_n[1] + dt*xDot_n[1];
 *      @endcode
 *    <td>
 *      @code
 *        // Take the timestep - Forward Euler
 *        Thyra::V_VpStV(x_np1.ptr(), *x_n, dt, *xDot_n);
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Other Thyra vector utilities.</b>
 *      We can also take advantage other Thyra features, Thyra::norm, to check
 *      if the solution has passed.
 *    <td>
 *      @code
 *        // Test if solution has passed.
 *        if ( std::isnan(x_n[0]) || std::isnan(x_n[1]) ) {
 *      @endcode
 *    <td>
 *      @code
 *        // Test if solution has passed.
 *        if ( std::isnan(Thyra::norm(*x_np1) ) {
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Use Thyra vector assignment.</b>
 *      The solution update is achieved via Thyra::V_V.
 *    <td>
 *      @code
 *          // Promote to next step (n <- n+1).
 *          x_n[0] = x_np1[0];
 *          x_n[1] = x_np1[1];
 *          n++;
 *      @endcode
 *    <td>
 *      @code
 *          // Promote to next step (n <- n+1).
 *          Thyra::V_V(x_n.ptr(), *x_np1);
 *          n++;
 *      @endcode
 *  </table>
 *
 *  The remainder of 01_Utilize_Thyra.cpp (e.g., regression testing) has
 *  similar changes to those from above.
 *
 *  \subsection example-01_Next Links
 *
 *  - Back to: \ref tutorials
 *  - Previous: \ref example-00
 *  - Next: \ref example-02
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
    int nTimeSteps = 2001;
    const double constDT = finalTime/(nTimeSteps-1);

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

    RCP<Thyra::VectorBase<double> > x_best = x_n->clone_v();
    {
      Thyra::DetachedVectorView<double> x_best_view(*x_best);
      x_best_view[0] = -1.59496108218721311;
      x_best_view[1] =  0.96359412806611255;
    }

    Thyra::V_VmV(x_error.ptr(), *x_n, *x_best);
           x_L2norm_error = Thyra::norm_2(*x_error);
    double x_L2norm_best  = Thyra::norm_2(*x_best );

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
