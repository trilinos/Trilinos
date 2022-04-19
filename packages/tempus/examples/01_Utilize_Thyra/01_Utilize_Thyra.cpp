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
 *  and utilizes Thyra vectors, replacing the standard double vectors,
 *  01_Utilize_Thyra.cpp.  Tempus uses Thyra for its Abstract Numerical
 *  Algorithms (ANAs), which are the mathematical concepts of
 *  vectors, vector spaces, and linear operators. All other ANA
 *  interfaces and support software are built on these fundamental
 *  operator/vector interfaces (including ModelEvaluators).
 *
 *  \subsection declaration Declaration and Initialization
 *
 *  We first need to create a vector space to construct the Thyra::Vector
 *  we need.  Thus we replace
 *  @code
 *    // Solution and its time-derivative.
 *    double x_n[2];      // at time index n
 *    double xDot_n[2];   // at time index n
 *  @endcode
 *  with
 *  @code
 *    // Solution and its time-derivative.
 *    int vectorLength = 2;  // number state unknowns
 *    RCP<const Thyra::VectorSpaceBase<double> >  xSpace =
 *      Thyra::defaultSpmdVectorSpace<double>(vectorLength);
 *
 *    RCP<Thyra::VectorBase<double> > x_n    = Thyra::createMember(xSpace);
 *    RCP<Thyra::VectorBase<double> > xDot_n = Thyra::createMember(xSpace);
 *  @endcode
 *
 *  The initialization can be achieved by replacing
 *  @code
 *    // Initial Conditions
 *    double time = 0.0;
 *    double epsilon = 1.0e-1;
 *    x_n   [0] = 2.0;
 *    x_n   [1] = 0.0;
 *    xDot_n[0] = 0.0;
 *    xDot_n[1] = -2.0/epsilon;
 *  @endcode
 *  with
 *  @code
 *    // Initial Conditions
 *    double time = 0.0;
 *    double epsilon = 1.0e-1;
 *    { // Scope to delete DetachedVectorViews
 *      Thyra::DetachedVectorView<double> x_n_view(*x_n);
 *      x_n_view[0] = 2.0;
 *      x_n_view[1] = 0.0;
 *      Thyra::DetachedVectorView<double> xDot_n_view(*xDot_n);
 *      xDot_n_view[0] = 0.0;
 *      xDot_n_view[1] = -2.0/epsilon;
 *    }
 *  @endcode
 *  and the scoped Thyra::DetachedVectorView's are used to delete them
 *  after use.
 *
 *  Elements of the thyra::Vector can be quickly accessed with get_ele by
 *  replacing
 *  @code
 *    cout << n << "  " << time << "  " << x_n[0] << "  " << x_n[1] << endl;
 *  @endcode
 *  with
 *  @code
 *    cout << n << "  " << time << "  " << get_ele(*(x_n), 0)
 *                              << "  " << get_ele(*(x_n), 1) << endl;
 *  @endcode
 *
 *  To initialize the solution at the next time step, we can simply
 *  replace
 *  @code
 *    // Initialize next time step
 *    double x_np1[2];    // at time index n+1
 *    x_np1[0] = x_n[0];
 *    x_np1[1] = x_n[1];
 *  @endcode
 *  with
 *  @code
 *    // Initialize next time step
 *    RCP<Thyra::VectorBase<double> > x_np1 = x_n->clone_v(); // at time index n+1
 *  @endcode
 *
 *  The model evaluation is again achieved through Thyra::DetachedVectorView
 *  by replacing
 *  @code
 *    // Righthand side evaluation and time-derivative at n.
 *    xDot_n[0] = x_n[1];
 *    xDot_n[1] = ((1.0 - x_n[0]*x_n[0])*x_n[1] - x_n[0])/epsilon;
 *  @endcode
 *  with
 *  @code
 *    // Righthand side evaluation and time-derivative at n.
 *    {
 *      Thyra::ConstDetachedVectorView<double> x_n_view(*x_n);
 *      Thyra::DetachedVectorView<double> xDot_n_view(*xDot_n);
 *      xDot_n_view[0] = x_n_view[1];
 *      xDot_n_view[1] =
 *        ((1.0-x_n_view[0]*x_n_view[0])*x_n_view[1]-x_n_view[0])/epsilon;
 *    }
 *  @endcode
 *
 *  The Forward Euler time step
 *  @code
 *    // Take the timestep - Forward Euler
 *    x_np1[0] = x_n[0] + dt*xDot_n[0];
 *    x_np1[1] = x_n[1] + dt*xDot_n[1];
 *  @endcode
 *  can be replaced with
 *  @code
 *    // Take the timestep - Forward Euler
 *    Thyra::V_VpStV(x_np1.ptr(), *x_n, dt, *xDot_n);
 *  @endcode
 *  where Thyra::V_VpStV performs an axpy.
 *
 *  We can also take advantage other Thyra features, Thyra::norm, to check
 *  if the solution is passing, i.e.,
 *  @code
 *    // Test if solution is passing.
 *    if ( std::isnan(x_n[0]) || std::isnan(x_n[1]) ) {
 *  @endcode
 *  is replaced with
 *  @code
 *    // Test if solution is passing.
 *    if ( std::isnan(Thyra::norm(*x_np1) ) {
 *  @endcode
 *
 *  The solution update is changed from
 *  @code
 *      // Promote to next step (n <- n+1).
 *      n++;
 *      x_n[0] = x_np1[0];
 *      x_n[1] = x_np1[1];
 *  @endcode
 *  to
 *  @code
 *      // Promote to next step (n <- n+1).
 *      n++;
 *      Thyra::V_V(x_n.ptr(), *x_np1);
 *  @endcode
 *
 *  The remainder of 01_Utilize_Thyra.cpp (e.g., regression testing) has
 *  similar changes to those from above.
 *
 *  \subsection example-01_Next Links
 *
 *  - Back to: \ref tutorials
 *  - Previous: \ref 00_Basic_Problem.cpp "Basic Problem"
 *  - Next: Use ModelEvaluator
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
    double time = 0.0;
    double epsilon = 1.0e-1;
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
    int n = 0;
    bool passing = true;
    cout << n << "  " << time << "  " << get_ele(*(x_n), 0)
                              << "  " << get_ele(*(x_n), 1) << endl;
    while (passing && time < finalTime && n < nTimeSteps) {

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

      // Test if solution is passing.
      if ( std::isnan(Thyra::norm(*x_np1)) ) {
        passing = false;
      } else {
        // Promote to next step (n <- n+1).
        n++;
        Thyra::V_V(x_n.ptr(), *x_np1);
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
      passing = false;
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
