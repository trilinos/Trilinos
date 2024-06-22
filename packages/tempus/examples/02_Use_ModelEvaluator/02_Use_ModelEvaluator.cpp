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

#include "VanDerPol_ModelEvaluator_02.hpp"


using namespace std;
using Teuchos::RCP;

/** @file */

/** \brief Description:
 *
 *  \section example-02 Example 2: Use ModelEvaluator
 *
 *  The next step in our progression is to use a Thyra::ModelEvaluator
 *  for the application problem, i.e., \ref vanderpol, and will
 *  provide a common interface to Abstract Numerical Algorithms
 *  (ANAs), e.g., nonlinear equations, NOX; stability analysis,
 *  LOCA; transient nonlinear ODEs, Tempus; and sensitivity analysis,
 *  Sacado.  For this example, we will be focusing on setting up a
 *  ModelEvaluator for transient nonlinear equations (e.g., explicit
 *  and implicit ODEs), and therefore this is not meant to be an
 *  indepth tutorial on ModelEvaluators.  Please refer to
 *  https://docs.trilinos.org/dev/packages/thyra/doc/html/classThyra_1_1ModelEvaluator.html
 *  and https://bartlettroscoe.github.io/publications/ModelEvaluator_HPCSW2008.pdf.
 *  for additional details on ModelEvaluators.
 *
 *  \subsection MEbasics ModelEvaluator Basics
 *
 *  Because the ModelEvaluator is a generic interface to all the ANAs,
 *  the primary function to evaluate the various quantities is simply
 *  @code
 *    Thyra::ModelEvaluator< Scalar >::evalModel(
 *      const ModelEvaluatorBase::InArgs< Scalar > & inArgs,
 *      const ModelEvaluatorBase::OutArgs< Scalar > & outArgs ) const
 *  @endcode
 *  where InArgs contains all the required information for the
 *  ModelEvaluator to compute the request from the ANA (i.e.,
 *  Tempus), and OutArgs will contain the computed quantities.
 *  In Tempus's case of explicit first-order ODEs (\f$ \dot{x} = f(x,t)\f$),
 *  the InArgs is just \f$x\f$ and \f$t\f$ and the OutArgs is
 *  \f$\dot{x}\f$.
 *  For implicit first-order ODEs (\f$ f(\dot{x},x,t) \f$),
 *  the InArgs is \f$ \dot{x} \f$, \f$x\f$ and \f$t\f$ and the OutArgs is
 *  the function evaluation, \f$ f(\dot{x},x,t) \f$.  Thus through the
 *  InArgs and OutArgs, the ModelEvaluator can determine which evaluation
 *  is being requested.
 *
 *  One important attribute of ModelEvaluators is they are
 *  <b>stateless</b>!  This means they do not storage any information
 *  that might change over the duration of the simulation, i.e., the
 *  ModelEvaluator does <b>not</b> contain any state information.
 *  Any state information must be passed in through the inputs (i.e.,
 *  InArgs).  A simple test for this attribute is the ModelEvaluator
 *  will return the exact same evaluation (i.e., OutArgs) for two
 *  consecutive evaluations, i.e.,
 *  @code
 *    evalModel(inArgs, outArgs1);
 *    evalModel(inArgs, outArgs2);
 *  @endcode
 *  and outArgs1 = outArgs2!
 *  The stateless attribute is very important for Tempus, because
 *  timesteps may fail and need to be retried with a new timestep size.
 *  If the ModelEvaluator is not stateless, then the re-evaluation will
 *  change!  The solvers (e.g., NOX) have a similar issue with multiple
 *  evaluations, and also require the stateless attribute.
 *
 *  For details on the ModelEvaluator for the \ref vanderpol implemented
 *  in this example, please review VanDerPol_ModelEvaluator_02.
 *
 *  \subsection changesToUseME Changes to Use ModelEvaluator
 *
 *  In the following table are code snippets of the changes from the
 *  \ref 01_Utilize_Thyra.cpp to 02_Use_ModelEvaluator.cpp.  This is
 *  similar taking a difference between them for the main() function.
 *
 *  <table>
 *    <tr> <th> Comments <th> "Utilize Thyra" Code Snippet
 *                       <th> "Use ModelEvaluator" Code Snippet
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Construct the ModelEvaluator.</b>
 *      The ModelEvaluator, VanDerPol_ModelEvaluator_02, now has all
 *      information to construct the problem, including the problem
 *      size, vectorLength, parameters, and initial conditions.
 *    <td>
 *      @code
 *        // Solution and its time-derivative.
 *        int vectorLength = 2;  // number state unknowns
 *        RCP<const Thyra::VectorSpaceBase<double> >  xSpace =
 *          Thyra::defaultSpmdVectorSpace<double>(vectorLength);
 *      @endcode
 *    <td>
 *      @code
 *        // Construct ModelEvaluator
 *        Teuchos::RCP<const Thyra::ModelEvaluator<double> >
 *          model = Teuchos::rcp(new VanDerPol_ModelEvaluator_02<double>());
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Get ICs from the ModelEvaluator.</b>
 *      The initial conditions (i.e., nominal values) have moved to
 *      the VanDerPol_ModelEvaluator_02 constructor, and can retrieved
 *      via getNominalValues().
 *    <td>
 *      @code
 *        RCP<Thyra::VectorBase<double> > x_n    = Thyra::createMember(xSpace);
 *        RCP<Thyra::VectorBase<double> > xDot_n = Thyra::createMember(xSpace);
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
 *    <td>
 *      @code
 *        // Get initial conditions from ModelEvaluator::getNominalValues.
 *        RCP<Thyra::VectorBase<double> > x_n    =
 *          model->getNominalValues().get_x()->clone_v();
 *        RCP<Thyra::VectorBase<double> > xDot_n =
 *          model->getNominalValues().get_x_dot()->clone_v();
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Setup InArgs and OutArgs for the ModelEvaluator.</b>
 *      Before calling evalModel(InArgs, OutArgs), the InArgs and
 *      OutArgs need to be setup.  For an explicit ODE, we need to
 *      set the time, the solution vector, x, and the time derivative,
 *      x_dot. In InArgs, x_dot needs to be set to null to indicate
 *      it is an explicit ODE evaluation.  Finally, we set the outArgs
 *      f evaluation to the time derivative, xDot, so it is filled
 *      with the result.
 *    <td>
 *      @code
 *      @endcode
 *    <td>
 *      @code
 *        // For explicit ODE formulation, xDot = f(x, t),
 *        // xDot is part of the outArgs.
 *        auto inArgs  = model->createInArgs();
 *        auto outArgs = model->createOutArgs();
 *        inArgs.set_t(time);
 *        inArgs.set_x(x_n);
 *        inArgs.set_x_dot(Teuchos::null);
 *        outArgs.set_f(xDot_n);
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Evaluate the ModelEvaluator.</b>
 *      The righthand-side evaluation is now in the ModelEvaluator,
 *      and is called with evalModel(InArgs, OutArgs).  Note that
 *      xDot has been filled in and can be accessed for the Forward-Euler
 *      update.
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
 *    <td>
 *      @code
 *        // Righthand side evaluation and time-derivative at n.
 *        model->evalModel(inArgs, outArgs);
 *      @endcode
 *  </table>
 *
 *  The remainder of 02_Use_ModelEvaluator.cpp is regression testing, etc.
 *
 *  \subsection example-02_Next Links
 *
 *  - Back to: \ref tutorials
 *  - Previous: \ref example-01
 *  - Next: \ref example-03
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
