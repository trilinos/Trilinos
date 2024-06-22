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


using namespace std;
using Teuchos::RCP;

/** @file */

/** \brief Description:
 *
 *  \section example-03 Example 3: Introduce SolutionState
 *
 *
 *  The basic idea behind the SolutionState object is that it
 *  contains all the information about the solution at a particular
 *  time, i.e., \f$x(t)\f$.  The primary requirement for SolutionState
 *  is that it contains all the information needed to restart the
 *  solution from scratch, e.g., check-pointing and recover from a
 *  failed timestep.  The four primary components of the
 *  Tempus::SolutionState are
 *  - The state - \f$x(t)\f$, \f$\dot{x}(t)\f$ (optional), and \f$\ddot{x}(t)\f$ (optional)
 *  - The metastate - timestep, time, order, error, etc.
 *  - The StepperState - any data needed by the stepper
 *  - The PhysicsState - any application-specific data
 *
 *  <b>For additional details about the Tempus::SolutionState, please see
 *  \ref SolutionState_Description!</b>
 *
 *  For this example, we will only utilize and demonstrate the state
 *  and the metastate (Tempus::SolutionStateMetaData).
 *
 *  \subsection changesToUseSolutionState Changes to Use SolutionState
 *
 *  In the following table are code snippets of the changes from the
 *  \ref 02_Use_ModelEvaluator.cpp to 03_Intro_SolutionState.cpp.  This is
 *  similar taking a difference between them for the main() function.
 *
 *  <table>
 *    <tr> <th> Comments <th> "Use ModelEvaluator" Code Snippet
 *                       <th> "Intro SolutionState" Code Snippet
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Initial the SolutionState.</b>
 *      Create a SolutionState object with the initial-condition
 *      solution, \f$x\f$, from the ModelEvaluator, using the one
 *      of the non-member constructors.  These non-member constructors
 *      set the remaining member data to default values (see
 *      Tempus::createSolutionStateX).  Member data can also be set
 *      through the basic set functions, e.g., setIndex(0).
 *
 *      <b>Note:</b> we did not include its time derivative, \f$\dot{x}\f$,
 *      as part of the SolutionState because it is not needed for
 *      output or diagnostics.  However we still require memory for it
 *      in the Forward Euler timestepping.  Foreshadowing: the time
 *      derivative can be handled and maintained by all the
 *      time steppers internally.
 *    <td>
 *      @code
 *        // Get initial conditions from ModelEvaluator::getNominalValues.
 *        int n = 0;
 *        double time = 0.0;
 *        RCP<Thyra::VectorBase<double> > x_n    =
 *          model->getNominalValues().get_x()->clone_v();
 *        RCP<Thyra::VectorBase<double> > xDot_n =
 *          model->getNominalValues().get_x_dot()->clone_v();
 *      @endcode
 *    <td>
 *      @code
 *        // Setup initial condition SolutionState --------------------
 *        auto solState = Tempus::createSolutionStateX(
 *                          model->getNominalValues().get_x()->clone_v());
 *        solState->setIndex   (0);
 *        solState->setTime    (0.0);
 *        solState->setTimeStep(0.0);  // By convention, the IC has dt=0.
 *        solState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are considered passed.
 *        RCP<Thyra::VectorBase<double> > xDot_n =
 *          model->getNominalValues().get_x_dot()->clone_v();
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Access member data of the SolutionState.</b>
 *      All the information about the solution and state are within the
 *      SolutionState object, and can be accessed through get functions,
 *      e.g., the index, time, solution, and solution status.
 *    <td>
 *      @code
 *        cout << n << "  " << time << "  " << get_ele(*(x_n), 0)
 *                                  << "  " << get_ele(*(x_n), 1) << endl;
 *        while (passed && time < finalTime && n < nTimeSteps) {
 *      @endcode
 *    <td>
 *      @code
 *        cout << solState->getIndex() << "  " << solState->getTime()
 *             << "  " << get_ele(*(solState->getX()), 0)
 *             << "  " << get_ele(*(solState->getX()), 1) << endl;
 *        while (solState->getSolutionStatus() == Tempus::Status::PASSED &&
 *               solState->getTime() < finalTime &&
 *               solState->getIndex() < nTimeSteps) {
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Setup the next timestep.</b>
 *      To initialize the next time step, we set some Referenced
 *      Counted Pointers (RCPs), and create a scratch vector.
 *      Additionally, we set the index, time and timestep for the
 *      next ("working") time step solution.  Although, not evident
 *      at this point, we will use the label "working" to help
 *      manage solutions that have not met the "passing" criteria,
 *      and we can try again with new settings, e.g., timestep size.
 *    <td>
 *      @code
 *        // Initialize next time step
 *        RCP<Thyra::VectorBase<double> > x_np1 = x_n->clone_v(); // at time index n+1

 *        // Set the timestep and time.
 *        double dt = constDT;
 *        time = (n+1)*dt;
 *      @endcode
 *    <td>
 *      @code
 *        // Initialize next time step
 *        RCP<Thyra::VectorBase<double> > x_n    = solState->getX();
 *        RCP<Thyra::VectorBase<double> > x_np1  = solState->getX()->clone_v(); // at time index n+1
 *        solState->setSolutionStatus(Tempus::Status::WORKING);

 *        // Set the timestep and time for the working solution i.e., n+1.
 *        int index = solState->getIndex()+1;
 *        double dt = constDT;
 *        double time = index*dt;

 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Set status of SolutionState.</b>
 *        Need to set SolutionState status to PASSED or FAILED, increment
 *        the index and time, and set the timestep.
 *    <td>
 *      @code
 *        // Test if solution has passed.
 *        if ( std::isnan(Thyra::norm(*x_np1)) ) {
 *          passed = false;
 *        } else {
 *          // Promote to next step (n <- n+1).
 *          Thyra::V_V(x_n.ptr(), *x_np1);
 *          n++;
 *        }
 *      @endcode
 *    <td>
 *      @code
 *        // Test if solution has passed.
 *        if ( std::isnan(Thyra::norm(*x_np1)) ) {
 *          solState->setSolutionStatus(Tempus::Status::FAILED);
 *        } else {
 *          // Promote to next step (n <- n+1).
 *          Thyra::V_V(x_n.ptr(), *x_np1);
 *          solState->setIndex   (index);
 *          solState->setTime    (time);
 *          solState->setTimeStep(constDT);
 *          solState->setSolutionStatus(Tempus::Status::PASSED);
 *        }
 *      @endcode
 *    <tr VALIGN=TOP>
 *    <td>
 *      <b>Access the SolutionState.</b>
 *      Need to access the SolutionState for the solution index.
 *      Note the RCP x_n still has access to the solution vector.
 *    <td>
 *      @code
 *        // Output
 *        if ( n%100 == 0 )
 *          cout << n << "  " << time << "  " << get_ele(*(x_n), 0)
 *                                    << "  " << get_ele(*(x_n), 1) << endl;
 *      @endcode
 *    <td>
 *      @code
 *        // Output
 *        if ( solState->getIndex()%100 == 0 )
 *          cout << solState->getIndex() << "  " << time
 *               << "  " << get_ele(*(x_n), 0)
 *               << "  " << get_ele(*(x_n), 1) << endl;
 *      @endcode
 *  </table>
 *
 *  The remainder of 03_Intro_SolutionState.cpp is regression testing, etc.
 *
 *  \subsection example-03_Next Links
 *
 *  - Back to: \ref tutorials
 *  - Previous: \ref example-02
 *  - Next: Adding SolutionHistory
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
    int nTimeSteps = 2001;
    const double constDT = finalTime/(nTimeSteps-1);

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
      if ( solState->getIndex()%100 == 0 )
        cout << solState->getIndex() << "  " << time
             << "  " << get_ele(*(x_n), 0)
             << "  " << get_ele(*(x_n), 1) << endl;
    }

    // Test for regression.
    RCP<Thyra::VectorBase<double> > x_n    = solState->getX();
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
      solState->setSolutionStatus(Tempus::Status::FAILED);
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
      solState->setSolutionStatus(Tempus::Status::FAILED);
      cout << "FAILED best constraint!" << endl;
    }
    if (solState->getSolutionStatus() == Tempus::Status::PASSED) success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if(success)
    cout << "\nEnd Result: Test Passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
