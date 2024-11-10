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
#include "Tempus_StepperTrapezoidal.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <vector>
#include <fstream>
#include <sstream>
#include <limits>

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
TEUCHOS_UNIT_TEST(Trapezoidal, ParameterList)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_Trapezoidal_SinCos.xml");

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  auto model                = rcp(new SinCosModel<double>(scm_pl));

  RCP<ParameterList> tempusPL = sublist(pList, "Tempus", true);

  // Test constructor IntegratorBasic(tempusPL, model)
  {
    RCP<Tempus::IntegratorBasic<double>> integrator =
        Tempus::createIntegratorBasic<double>(tempusPL, model);

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);

    RCP<const ParameterList> defaultPL =
        integrator->getStepper()->getValidParameters();
    bool pass = haveSameValuesSorted(*stepperPL, *defaultPL, true);
    if (!pass) {
      out << std::endl;
      out << "stepperPL -------------- \n"
          << *stepperPL << std::endl;
      out << "defaultPL -------------- \n"
          << *defaultPL << std::endl;
    }
    TEST_ASSERT(pass)
  }

  // Test constructor IntegratorBasic(model, stepperType)
  {
    RCP<Tempus::IntegratorBasic<double>> integrator =
        Tempus::createIntegratorBasic<double>(
            model, std::string("Trapezoidal Method"));

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
    RCP<const ParameterList> defaultPL =
        integrator->getStepper()->getValidParameters();

    bool pass = haveSameValuesSorted(*stepperPL, *defaultPL, true);
    if (!pass) {
      out << std::endl;
      out << "stepperPL -------------- \n"
          << *stepperPL << std::endl;
      out << "defaultPL -------------- \n"
          << *defaultPL << std::endl;
    }
    TEST_ASSERT(pass)
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Trapezoidal, ConstructingFromDefaults)
{
  double dt = 0.1;
  std::vector<std::string> options;
  options.push_back("Default Parameters");
  options.push_back("ICConsistency and Check");

  for (const auto& option : options) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_Trapezoidal_SinCos.xml");
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    // RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
    auto model = rcp(new SinCosModel<double>(scm_pl));

    // Setup Stepper for field solve ----------------------------
    auto stepper = rcp(new Tempus::StepperTrapezoidal<double>());
    //{
    //  // Turn on NOX output.
    //  RCP<ParameterList> sPL = stepper->getNonconstParameterList();
    //  std::string solverName = sPL->get<std::string>("Solver Name");
    //  RCP<ParameterList> solverPL = Teuchos::sublist(sPL, solverName, true);
    //  solverPL->sublist("NOX").sublist("Printing")
    //           .sublist("Output Information").set("Outer Iteration", true);
    //  solverPL->sublist("NOX").sublist("Printing")
    //           .sublist("Output Information").set("Parameters", true);
    //  solverPL->sublist("NOX").sublist("Printing")
    //           .sublist("Output Information").set("Details", true);
    //  stepper->setSolver(solverPL);
    //}
    stepper->setModel(model);
    if (option == "ICConsistency and Check") {
      stepper->setICConsistency("Consistent");
      stepper->setICConsistencyCheck(true);
    }
    stepper->initialize();

    // Setup TimeStepControl ------------------------------------
    auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
    ParameterList tscPL =
        pl->sublist("Default Integrator").sublist("Time Step Control");
    timeStepControl->setInitIndex(tscPL.get<int>("Initial Time Index"));
    timeStepControl->setInitTime(tscPL.get<double>("Initial Time"));
    timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
    timeStepControl->setInitTimeStep(dt);
    timeStepControl->initialize();

    // Setup initial condition SolutionState --------------------
    auto inArgsIC = model->getNominalValues();
    auto icSoln   = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
    auto icSolnDot =
        rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x_dot());
    auto icState = Tempus::createSolutionStateX(icSoln, icSolnDot);
    icState->setTime(timeStepControl->getInitTime());
    icState->setIndex(timeStepControl->getInitIndex());
    icState->setTimeStep(0.0);
    icState->setOrder(stepper->getOrder());
    icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

    // Setup SolutionHistory ------------------------------------
    auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
    solutionHistory->setName("Forward States");
    solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
    solutionHistory->setStorageLimit(2);
    solutionHistory->addState(icState);

    // Setup Integrator -----------------------------------------
    RCP<Tempus::IntegratorBasic<double>> integrator =
        Tempus::createIntegratorBasic<double>();
    integrator->setStepper(stepper);
    integrator->setTimeStepControl(timeStepControl);
    integrator->setSolutionHistory(solutionHistory);
    // integrator->setObserver(...);
    integrator->initialize();

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time      = integrator->getTime();
    double timeFinal = pl->sublist("Default Integrator")
                           .sublist("Time Step Control")
                           .get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double>> x = integrator->getX();
    RCP<const Thyra::VectorBase<double>> x_exact =
        model->getExactSolution(time).get_x();

    // Calculate the error
    RCP<Thyra::VectorBase<double>> xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));

    // Check the order and intercept
    out << "  Stepper = " << stepper->description() << "\n            with "
        << option << std::endl;
    out << "  =========================" << std::endl;
    out << "  Exact solution   : " << get_ele(*(x_exact), 0) << "   "
        << get_ele(*(x_exact), 1) << std::endl;
    out << "  Computed solution: " << get_ele(*(x), 0) << "   "
        << get_ele(*(x), 1) << std::endl;
    out << "  Difference       : " << get_ele(*(xdiff), 0) << "   "
        << get_ele(*(xdiff), 1) << std::endl;
    out << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 0.841021, 1.0e-4);
    TEST_FLOATING_EQUALITY(get_ele(*(x), 1), 0.541002, 1.0e-4);
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Trapezoidal, SinCos)
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
  for (int n = 0; n < nTimeStepSizes; n++) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_Trapezoidal_SinCos.xml");

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
    pl->sublist("Default Integrator")
        .sublist("Time Step Control")
        .set("Initial Time Step", dt);
    integrator = Tempus::createIntegratorBasic<double>(pl, model);

    // Initial Conditions
    // During the Integrator construction, the initial SolutionState
    // is set by default to model->getNominalVales().get_x().  However,
    // the application can set it also by integrator->initializeSolutionHistory.
    {
      RCP<Thyra::VectorBase<double>> x0 =
          model->getNominalValues().get_x()->clone_v();
      RCP<Thyra::VectorBase<double>> xdot0 =
          model->getNominalValues().get_x_dot()->clone_v();
      integrator->initializeSolutionHistory(0.0, x0, xdot0);
    }

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
      writeSolution("Tempus_Trapezoidal_SinCos.dat", solutionHistory);

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
      writeSolution("Tempus_Trapezoidal_SinCos-Ref.dat", solnHistExact);
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

  double xSlope                        = 0.0;
  double xDotSlope                     = 0.0;
  RCP<Tempus::Stepper<double>> stepper = integrator->getStepper();
  double order                         = stepper->getOrder();
  writeOrderError("Tempus_Trapezoidal_SinCos-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                  xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.01);
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.000832086, 1.0e-4);
  TEST_FLOATING_EQUALITY(xDotSlope, order, 0.01);
  TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 0.000832086, 1.0e-4);

  Teuchos::TimeMonitor::summarize();
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Trapezoidal, VanDerPol)
{
  RCP<Tempus::IntegratorBasic<double>> integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 4;
  double dt                = 0.05;
  double time              = 0.0;
  for (int n = 0; n < nTimeStepSizes; n++) {
    // Read params from .xml file
    RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_Trapezoidal_VanDerPol.xml");

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
    time             = integrator->getTime();
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
      std::string fname = "Tempus_Trapezoidal_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_Trapezoidal_VanDerPol.dat";
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
  writeOrderError("Tempus_Trapezoidal_VanDerPol-Error.dat", stepper, StepSize,
                  solutions, xErrorNorm, xSlope, solutionsDot, xDotErrorNorm,
                  xDotSlope, out);

  TEST_FLOATING_EQUALITY(xSlope, order, 0.10);
  TEST_FLOATING_EQUALITY(xDotSlope, order, 0.10);  //=order at samll dt
  TEST_FLOATING_EQUALITY(xErrorNorm[0], 0.00085855, 1.0e-4);
  TEST_FLOATING_EQUALITY(xDotErrorNorm[0], 0.00120695, 1.0e-4);

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
