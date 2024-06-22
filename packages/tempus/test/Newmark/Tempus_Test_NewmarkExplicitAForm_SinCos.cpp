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

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"

#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperNewmarkImplicitAForm.hpp"
#include "Tempus_StepperNewmarkImplicitDForm.hpp"

#include "../TestModels/HarmonicOscillatorModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::getParametersFromXmlFile;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// ************************************************************
TEUCHOS_UNIT_TEST(NewmarkExplicitAForm, SinCos)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 9;
  double time              = 0.0;

  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_Test_NewmarkExplicitAForm_SinCos.xml");

  // Setup the HarmonicOscillatorModel
  RCP<ParameterList> hom_pl = sublist(pList, "HarmonicOscillatorModel", true);
  RCP<HarmonicOscillatorModel<double>> model =
      Teuchos::rcp(new HarmonicOscillatorModel<double>(hom_pl));

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl        = sublist(pList, "Tempus", true);
  RCP<ParameterList> stepperPL = sublist(pl, "Default Stepper", true);
  stepperPL->remove("Zero Initial Guess");

  // Set initial time step = 2*dt specified in input file (for convergence
  // study)
  //
  double dt = pl->sublist("Default Integrator")
                  .sublist("Time Step Control")
                  .get<double>("Initial Time Step");
  dt *= 2.0;

  for (int n = 0; n < nTimeStepSizes; n++) {
    // Perform time-step refinement
    dt /= 2;
    out << "\n \n time step #" << n << " (out of " << nTimeStepSizes - 1
        << "), dt = " << dt << "\n";
    pl->sublist("Default Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);

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
      writeSolution("Tempus_Test_NewmarkExplicitAForm_SinCos.dat",
                    solutionHistory);

      RCP<Tempus::SolutionHistory<double>> solnHistExact =
          Teuchos::rcp(new Tempus::SolutionHistory<double>());
      for (int i = 0; i < solutionHistory->getNumStates(); i++) {
        double time_i                            = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double>> state = Tempus::createSolutionStateX(
            rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x()),
            rcp_const_cast<Thyra::VectorBase<double>>(
                model->getExactSolution(time_i).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_Test_NewmarkExplicitAForm_SinCos-Ref.dat",
                    solnHistExact);
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
  double xSlope                        = 0.0;
  double xDotSlope                     = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  double order                         = stepper->getOrder();
  writeOrderError("Tempus_Test_NewmarkExplicitAForm_SinCos-Error.dat", stepper,
                  StepSize, solutions, xErrorNorm, xSlope, solutionsDot,
                  xDotErrorNorm, xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.02);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.0157928, 1.0e-4);
  TEST_FLOATING_EQUALITY(xDotSlope, order, 0.02);
  TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 0.233045, 1.0e-4);

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
