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

/** \page example-03 Example 3: Introduce SolutionState
 *
 *  This example introduces \ref Tempus::SolutionState, which stores the
 *  solution and selected metadata at a particular time.  The model is still
 *  provided through a \ref Thyra::ModelEvaluator, and the Forward Euler update
 *  is still written explicitly in the application code, but the evolving state
 *  of the integration is now organized using a core \ref Tempus data structure.
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
 *  The central idea behind \ref Tempus::SolutionState is that it should contain
 *  the information needed to represent a solution at a specific time in a form
 *  suitable for restart, checkpointing, diagnostics, and recovery from failed
 *  timesteps.
 *
 *  This example uses only part of the full \ref Tempus::SolutionState
 *  capability:
 *  - the solution vector \f$x(t)\f$
 *  - selected metadata from \ref Tempus::SolutionStateMetaData
 *  - For additional details, see \ref SolutionState_Description.
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
    if (solState->getSolutionStatus() == Tempus::Status::PASSED) success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  if(success)
    cout << "\nEnd Result: Test Passed!" << std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
