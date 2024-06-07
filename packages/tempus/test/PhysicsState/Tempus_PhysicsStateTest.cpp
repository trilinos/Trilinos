//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_PhysicsStateTest_StepperForwardEuler.hpp"
#include "Tempus_PhysicsStateCounter.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::getParametersFromXmlFile;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::sublist;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(PhysicsState, SinCos)
{
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 7;
  double dt                = 0.2;
  double order             = 0.0;
  for (int n = 0; n < nTimeStepSizes; n++) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_PhysicsState_SinCos.xml");

    // std::ofstream ftmp("PL.txt");
    // pList->print(ftmp);
    // ftmp.close();

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    // RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
    RCP<SinCosModel<double> > model =
        Teuchos::rcp(new SinCosModel<double>(scm_pl));

    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::createIntegratorBasic<double>(pl, model);

    // Replace Tempus::StepperForwardEuler with
    // Tempus_Test::StepperPhysicsStateTest
    Teuchos::RCP<Tempus::Stepper<double> > physicsStepper =
        Teuchos::rcp(new StepperPhysicsStateTest<double>(model));
    integrator->setStepper(physicsStepper);
    order = integrator->getStepper()->getOrder();

    // Initial Conditions
    // During the Integrator construction, the initial SolutionState
    // is set by default to model->getNominalVales().get_x().  However,
    // the application can set it also by integrator->initializeSolutionHistory.
    RCP<Thyra::VectorBase<double> > x0 =
        model->getNominalValues().get_x()->clone_v();
    integrator->initializeSolutionHistory(0.0, x0);

    // Replace Tempus::PhysicsState with
    // Tempus_Test::PhysicsStateCounter
    RCP<PhysicsStateCounter<double> > pSC =
        Teuchos::rcp(new PhysicsStateCounter<double>("PhysicsStateTest", 0));
    integrator->getSolutionHistory()->getCurrentState()->setPhysicsState(pSC);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test PhysicsState
    Teuchos::RCP<Tempus::PhysicsState<double> > pS =
        integrator->getSolutionHistory()->getCurrentState()->getPhysicsState();
    TEST_EQUALITY(pS->getName(), "PhysicsStateTest");
    pSC = Teuchos::rcp_dynamic_cast<PhysicsStateCounter<double> >(pS);
    // out << " Counter = " << pSC->getCounter() << std::endl;
    TEST_EQUALITY(pSC->getCounter(), (int)(1.0 / dt));

    // Test if at 'Final Time'
    double time      = integrator->getTime();
    double timeFinal = pl->sublist("Demo Integrator")
                           .sublist("Time Step Control")
                           .get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::VectorBase<double> > x_exact =
        model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution
    if (n == 0) {
      std::ofstream ftmp("Tempus_ForwardEuler_SinCos.dat");
      RCP<const SolutionHistory<double> > solutionHistory =
          integrator->getSolutionHistory();
      RCP<const Thyra::VectorBase<double> > x_exact_plot;
      for (int i = 0; i < solutionHistory->getNumStates(); i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        double time_i                                   = solutionState->getTime();
        RCP<const Thyra::VectorBase<double> > x_plot    = solutionState->getX();
        x_exact_plot                                    = model->getExactSolution(time_i).get_x();
        ftmp << time_i << "   " << Thyra::get_ele(*(x_plot), 0) << "   "
             << Thyra::get_ele(*(x_plot), 1) << "   "
             << Thyra::get_ele(*(x_exact_plot), 0) << "   "
             << Thyra::get_ele(*(x_exact_plot), 1) << std::endl;
      }
      ftmp.close();
    }

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
    StepSize.push_back(dt);
    const double L2norm = Thyra::norm_2(*xdiff);
    ErrorNorm.push_back(L2norm);
  }

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
  out << "  Stepper = ForwardEuler" << std::endl;
  out << "  =========================" << std::endl;
  out << "  Expected order: " << order << std::endl;
  out << "  Observed order: " << slope << std::endl;
  out << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(slope, order, 0.01);
  TEST_FLOATING_EQUALITY(ErrorNorm[0], 0.051123, 1.0e-4);

  std::ofstream ftmp("Tempus_ForwardEuler_SinCos-Error.dat");
  double error0 = 0.8 * ErrorNorm[0];
  for (int n = 0; n < nTimeStepSizes; n++) {
    ftmp << StepSize[n] << "   " << ErrorNorm[n] << "   "
         << error0 * (pow(StepSize[n] / StepSize[0], order)) << std::endl;
  }
  ftmp.close();

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
