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
TEUCHOS_UNIT_TEST(BDF2, SinCos)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;

  // Read params from .xml file
  RCP<ParameterList> pList = getParametersFromXmlFile("Tempus_BDF2_SinCos.xml");
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
  std::string ref_out_file_name = output_file_string + "-Ref.dat";
  std::string err_out_file_name = output_file_string + "-Error.dat";
  double time                   = 0.0;
  for (int n = 0; n < nTimeStepSizes; n++) {
    auto model = rcp(new SinCosModel<double>(scm_pl));

    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> tempusPL =
        getParametersFromXmlFile("Tempus_BDF2_SinCos.xml");
    RCP<ParameterList> pl = sublist(tempusPL, "Tempus", true);
    pl->sublist("Default Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
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

    // Plot sample solution and exact solution
    if (n == 0) {
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      writeSolution(output_file_name, solutionHistory);

      auto solnHistExact = rcp(new Tempus::SolutionHistory<double>());
      for (int i = 0; i < solutionHistory->getNumStates(); i++) {
        double time_i = (*solutionHistory)[i]->getTime();
        auto state    = Tempus::createSolutionStateX(
               rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x()),
               rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution(ref_out_file_name, solnHistExact);
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
    double order                         = stepper->getOrder();
    writeOrderError(err_out_file_name, stepper, StepSize, solutions, xErrorNorm,
                    xSlope, solutionsDot, xDotErrorNorm, xDotSlope, out);

    TEST_FLOATING_EQUALITY(xSlope, order, 0.01);
    TEST_FLOATING_EQUALITY(xDotSlope, order, 0.01);
    TEST_FLOATING_EQUALITY(xErrorNorm[0], 5.13425e-05, 1.0e-4);
    TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 5.13425e-05, 1.0e-4);
  }

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
