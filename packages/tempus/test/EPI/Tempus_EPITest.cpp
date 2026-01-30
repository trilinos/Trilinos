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
#include "Thyra_DetachedVectorView.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include "Tempus_StepperEPI.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
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
// TEUCHOS_UNIT_TEST(EPI, ParameterList)
// {
//   // Read params from .xml file
//   RCP<ParameterList> pList =
//       getParametersFromXmlFile("Tempus_EPI_SinCos.xml");

//   // Setup the SinCosModel
//   RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
//   auto model                = rcp(new SinCosModel<double>(scm_pl));

//   RCP<ParameterList> tempusPL = sublist(pList, "Tempus", true);

//   // Test constructor IntegratorBasic(tempusPL, model)
//   {
//     RCP<Tempus::IntegratorBasic<double>> integrator =
//         Tempus::createIntegratorBasic<double>(tempusPL, model);

//     RCP<ParameterList> stepperPL = sublist(tempusPL, "Demo Stepper", true);
//     RCP<const ParameterList> defaultPL =
//         integrator->getStepper()->getValidParameters();

//     bool pass = haveSameValuesSorted(*stepperPL, *defaultPL, true);
//     if (!pass) {
//       out << std::endl;
//       out << "stepperPL -------------- \n"
//           << *stepperPL << std::endl;
//       out << "defaultPL -------------- \n"
//           << *defaultPL << std::endl;
//     }
//     TEST_ASSERT(pass)
//   }

//   // Test constructor IntegratorBasic(model, stepperType)
//   {
//     RCP<Tempus::IntegratorBasic<double>> integrator =
//         Tempus::createIntegratorBasic<double>(model,
//                                               std::string("EPI"));

//     RCP<ParameterList> stepperPL = sublist(tempusPL, "Demo Stepper", true);
//     RCP<const ParameterList> defaultPL =
//         integrator->getStepper()->getValidParameters();

//     bool pass = haveSameValuesSorted(*stepperPL, *defaultPL, true);
//     if (!pass) {
//       out << std::endl;
//       out << "stepperPL -------------- \n"
//           << *stepperPL << std::endl;
//       out << "defaultPL -------------- \n"
//           << *defaultPL << std::endl;
//     }
//     TEST_ASSERT(pass)
//   }
// }

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, SinCos)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 7;
  double dt                = 0.2;
  double time              = 0.0;
  const int expected_order = 3;
  for (int n = 0; n < nTimeStepSizes; n++) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_EPI_SinCos.xml");

    // std::ofstream ftmp("PL.txt");
    // pList->print(ftmp);
    // ftmp.close();

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    // RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
    auto model = rcp(new SinCosModel<double>(scm_pl));

    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
    // The reason why this is set to 3 in here rather than xml file is that
    // xml file values are compared with the default parameter list values.
    pl->sublist("Demo Stepper")
        .set("Taylor Expansion Order", expected_order);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);
    // Initial Conditions
    // During the Integrator construction, the initial SolutionState
    // is set by default to model->getNominalVales().get_x().  However,
    // the application can set it also by integrator->initializeSolutionHistory.
    RCP<Thyra::VectorBase<double>> x0 =
        model->getNominalValues().get_x()->clone_v();

    integrator->initializeSolutionHistory(0.0, x0);
    integrator->initialize();

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test PhysicsState
    RCP<Tempus::PhysicsState<double>> physicsState =
        integrator->getSolutionHistory()->getCurrentState()->getPhysicsState();
    TEST_EQUALITY(physicsState->getName(), "Tempus::PhysicsState");

    // Test if at 'Final Time'
    time             = integrator->getTime();
    double timeFinal = pl->sublist("Demo Integrator")
                           .sublist("Time Step Control")
                           .get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double>> x = integrator->getX();
    RCP<const Thyra::VectorBase<double>> x_exact =
        model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution
    if (n == 0) {
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      writeSolution("Tempus_EPI_SinCos.dat", solutionHistory);

    //   solutionHistory->printHistory("high");

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
      writeSolution("Tempus_EPI_SinCos-Ref.dat", solnHistExact);
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
  // This is a linear problem, so the expected order is 3 due to 
  // Taylor expansion order regardless of the integration order.
  double order                         = (double)expected_order;
  writeOrderError("Tempus_EPI_SinCos-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                  xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.01);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 4.16603e-05, 1.0e-5);

  Teuchos::TimeMonitor::summarize();
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(EPI, VanDerPol)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 7;  // 8 for Error plot
  double dt                = 0.2;
  for (int n = 0; n < nTimeStepSizes; n++) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_EPI_VanDerPol.xml");

    // Setup the VanDerPolModel
    RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
    auto model                 = rcp(new VanDerPolModel<double>(vdpm_pl));

    // Set the step size
    dt /= 2;
    if (n == nTimeStepSizes - 1) dt /= 10.0;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time      = integrator->getTime();
    double timeFinal = pl->sublist("Demo Integrator")
                           .sublist("Time Step Control")
                           .get<double>("Final Time");
    double tol = 100.0 * std::numeric_limits<double>::epsilon();
    TEST_FLOATING_EQUALITY(time, timeFinal, tol);

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()), solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXDot()), solutionDot.ptr());
    solutionsDot.push_back(solutionDot);

    // Output finest temporal solution for plotting
    // This only works for ONE MPI process
    if ((n == 0) || (n == nTimeStepSizes - 1)) {
      std::string fname = "Tempus_EPI_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_EPI_VanDerPol.dat";
      RCP<const SolutionHistory<double>> solutionHistory =
          integrator->getSolutionHistory();
      writeSolution(fname, solutionHistory);
    }
  }

  // Check the order and intercept
  double xSlope                        = 0.0;
  double xDotSlope                     = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  double order                         = stepper->getOrder();
  writeOrderError("Tempus_EPI_VanDerPol-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                  xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.15);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.00210467, 1.0e-5);

  Teuchos::TimeMonitor::summarize();
}

// // ************************************************************
// // ************************************************************
// TEUCHOS_UNIT_TEST(EPI, NumberTimeSteps)
// {
//   std::vector<double> StepSize;
//   std::vector<double> ErrorNorm;
//   // const int nTimeStepSizes = 7;
//   // double dt = 0.2;
//   // double order = 0.0;

//   // Read params from .xml file
//   RCP<ParameterList> pList =
//       getParametersFromXmlFile("Tempus_EPI_NumberOfTimeSteps.xml");

//   // Setup the VanDerPolModel
//   RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
//   auto model                 = rcp(new VanDerPolModel<double>(vdpm_pl));

//   // Setup the Integrator and reset initial time step
//   RCP<ParameterList> pl = sublist(pList, "Tempus", true);

//   // dt = pl->sublist("Demo Integrator")
//   //         .sublist("Time Step Control")
//   //         .get<double>("Initial Time Step");
//   const int numTimeSteps = pl->sublist("Demo Integrator")
//                                .sublist("Time Step Control")
//                                .get<int>("Number of Time Steps");

//   RCP<Tempus::IntegratorBasic<double>> integrator =
//       Tempus::createIntegratorBasic<double>(pl, model);

//   // Integrate to timeMax
//   bool integratorStatus = integrator->advanceTime();
//   TEST_ASSERT(integratorStatus)

//   // check that the number of time steps taken is whats is set
//   // in the parameter list
//   TEST_EQUALITY(numTimeSteps, integrator->getIndex());
// }

// // ************************************************************
// // ************************************************************
// TEUCHOS_UNIT_TEST(EPI, Variable_TimeSteps)
// {
//   // Read params from .xml file
//   RCP<ParameterList> pList =
//       getParametersFromXmlFile("Tempus_EPI_VanDerPol.xml");

//   // Setup the VanDerPolModel
//   RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
//   auto model                 = rcp(new VanDerPolModel<double>(vdpm_pl));

//   // Setup the Integrator and reset initial time step
//   RCP<ParameterList> pl = sublist(pList, "Tempus", true);

//   // Set parameters for this test.
//   pl->sublist("Demo Integrator")
//       .sublist("Time Step Control")
//       .set("Initial Time Step", 0.01);

//   pl->sublist("Demo Integrator")
//       .sublist("Time Step Control")
//       .sublist("Time Step Control Strategy")
//       .set("Reduction Factor", 0.9);
//   pl->sublist("Demo Integrator")
//       .sublist("Time Step Control")
//       .sublist("Time Step Control Strategy")
//       .set("Amplification Factor", 1.15);
//   pl->sublist("Demo Integrator")
//       .sublist("Time Step Control")
//       .sublist("Time Step Control Strategy")
//       .set("Minimum Value Monitoring Function", 0.05);
//   pl->sublist("Demo Integrator")
//       .sublist("Time Step Control")
//       .sublist("Time Step Control Strategy")
//       .set("Maximum Value Monitoring Function", 0.1);

//   pl->sublist("Demo Integrator")
//       .sublist("Solution History")
//       .set("Storage Type", "Static");
//   pl->sublist("Demo Integrator")
//       .sublist("Solution History")
//       .set("Storage Limit", 3);

//   RCP<Tempus::IntegratorBasic<double>> integrator =
//       Tempus::createIntegratorBasic<double>(pl, model);

//   // Integrate to timeMax
//   bool integratorStatus = integrator->advanceTime();
//   TEST_ASSERT(integratorStatus)

//   // Check 'Final Time'
//   double time      = integrator->getTime();
//   double timeFinal = pl->sublist("Demo Integrator")
//                          .sublist("Time Step Control")
//                          .get<double>("Final Time");
//   TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

//   // Check TimeStep size
//   auto state = integrator->getCurrentState();
//   double dt  = state->getTimeStep();
//   TEST_FLOATING_EQUALITY(dt, 0.008310677297208358, 1.0e-12);

//   // Check number of time steps taken
//   const int numTimeSteps = 60;
//   TEST_EQUALITY(numTimeSteps, integrator->getIndex());

//   // Time-integrated solution and the reference solution
//   RCP<Thyra::VectorBase<double>> x     = integrator->getX();
//   RCP<Thyra::VectorBase<double>> x_ref = x->clone_v();
//   {
//     Thyra::DetachedVectorView<double> x_ref_view(*x_ref);
//     x_ref_view[0] = -1.931946840284863;
//     x_ref_view[1] = 0.645346748303107;
//   }

//   // Calculate the error
//   RCP<Thyra::VectorBase<double>> xdiff = x->clone_v();
//   Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_ref, -1.0, *(x));

//   // Check the solution
//   out << "  Stepper = EPI" << std::endl;
//   out << "  =========================" << std::endl;
//   out << "  Reference solution: " << get_ele(*(x_ref), 0) << "   "
//       << get_ele(*(x_ref), 1) << std::endl;
//   out << "  Computed solution : " << get_ele(*(x), 0) << "   "
//       << get_ele(*(x), 1) << std::endl;
//   out << "  Difference        : " << get_ele(*(xdiff), 0) << "   "
//       << get_ele(*(xdiff), 1) << std::endl;
//   out << "  =========================" << std::endl;
//   TEST_FLOATING_EQUALITY(get_ele(*(x), 0), get_ele(*(x_ref), 0), 1.0e-12);
//   TEST_FLOATING_EQUALITY(get_ele(*(x), 1), get_ele(*(x_ref), 1), 1.0e-12);
// }

}  // namespace Tempus_Test
