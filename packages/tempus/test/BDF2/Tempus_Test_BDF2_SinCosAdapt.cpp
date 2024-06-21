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
#include "Teuchos_DefaultComm.hpp"

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_StepperBDF2.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::getParametersFromXmlFile;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::sublist;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BDF2, SinCosAdapt)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BDF2_SinCos_AdaptDt.xml");
  // Set initial time step = 2*dt specified in input file (for convergence
  // study)
  double dt = pList->sublist("Tempus")
                  .sublist("Default Integrator")
                  .sublist("Time Step Control")
                  .get<double>("Initial Time Step");
  dt *= 2.0;

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  const int nTimeStepSizes  = scm_pl->get<int>("Number of Time Step Sizes", 7);
  std::string output_file_string =
      scm_pl->get<std::string>("Output File Name", "Tempus_BDF2_SinCos");
  std::string output_file_name  = output_file_string + ".dat";
  std::string err_out_file_name = output_file_string + "-Error.dat";
  double time                   = 0.0;
  for (int n = 0; n < nTimeStepSizes; n++) {
    auto model = rcp(new SinCosModel<double>(scm_pl));

    dt /= 2;

    RCP<ParameterList> tempusPL =
        getParametersFromXmlFile("Tempus_BDF2_SinCos_AdaptDt.xml");
    RCP<ParameterList> pl = sublist(tempusPL, "Tempus", true);

    // Setup the Integrator and reset initial time step
    pl->sublist("Default Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt / 4.0);
    // Ensure time step does not get larger than the initial time step size,
    // as that would mess up the convergence rates.
    pl->sublist("Default Integrator")
        .sublist("Time Step Control")
        .set("Maximum Time Step", dt);
    // Ensure time step does not get too small and therefore too many steps.
    pl->sublist("Default Integrator")
        .sublist("Time Step Control")
        .set("Minimum Time Step", dt / 4.0);
    // For the SinCos problem eta is directly related to dt
    pl->sublist("Default Integrator")
        .sublist("Time Step Control")
        .sublist("Time Step Control Strategy")
        .set("Minimum Value Monitoring Function", dt * 0.99);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);

    // Initial Conditions
    // During the Integrator construction, the initial SolutionState
    // is set by default to model->getNominalVales().get_x().  However,
    // the application can set it also by integrator->initializeSolutionHistory.
    RCP<Thyra::VectorBase<double>> x0 =
        model->getNominalValues().get_x()->clone_v();
    integrator->initializeSolutionHistory(0.0, x0);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    time             = integrator->getTime();
    double timeFinal = pl->sublist("Default Integrator")
                           .sublist("Time Step Control")
                           .get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double>> x = integrator->getX();
    RCP<const Thyra::VectorBase<double>> x_exact =
        model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution
    if (n == 0) {
      std::ofstream ftmp(output_file_name);
      // Warning: the following assumes serial run
      FILE *gold_file = fopen("Tempus_BDF2_SinCos_AdaptDt_gold.dat", "r");
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      RCP<const Thyra::VectorBase<double>> x_exact_plot;
      for (int i = 0; i < solutionHistory->getNumStates(); i++) {
        char time_gold_char[100];
        fgets(time_gold_char, 100, gold_file);
        double time_gold;
        sscanf(time_gold_char, "%lf", &time_gold);
        RCP<const SolutionState<double>> solutionState = (*solutionHistory)[i];
        double time_i                                  = solutionState->getTime();
        // Throw error if time does not match time in gold file to specified
        // tolerance
        TEST_FLOATING_EQUALITY(time_i, time_gold, 1.0e-5);
        RCP<const Thyra::VectorBase<double>> x_plot = solutionState->getX();
        x_exact_plot                                = model->getExactSolution(time_i).get_x();
        ftmp << time_i << "   " << get_ele(*(x_plot), 0) << "   "
             << get_ele(*(x_plot), 1) << "   " << get_ele(*(x_exact_plot), 0)
             << "   " << get_ele(*(x_exact_plot), 1) << std::endl;
      }
      ftmp.close();
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()), solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()), solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes - 1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solutionExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),
                  solutionExact.ptr());
      solutions.push_back(solutionExact);
      auto solutionDotExact = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDotExact.ptr());
      solutionsDot.push_back(solutionDotExact);
    }
  }

  // Check the order and intercept
  if (nTimeStepSizes > 1) {
    double xSlope    = 0.0;
    double xDotSlope = 0.0;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;
    RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
    // double order = stepper->getOrder();
    writeOrderError("Tempus_BDF2_SinCos-Error.dat", stepper, StepSize,
                    solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                    xDotSlope, out);

    TEST_FLOATING_EQUALITY(xSlope, 1.932, 0.01);
    TEST_FLOATING_EQUALITY(xDotSlope, 1.932, 0.01);
    TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.000192591, 1.0e-4);
    TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 0.000192591, 1.0e-4);
  }

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
