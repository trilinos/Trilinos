// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include "Tempus_StepperForwardEuler.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;


// Comment out any of the following tests to exclude from build/run.
#define TEST_PARAMETERLIST
#define TEST_CONSTRUCTING_FROM_DEFAULTS
#define TEST_SINCOS
#define TEST_VANDERPOL
#define TEST_NUMBER_TIMESTEPS


#ifdef TEST_PARAMETERLIST
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, ParameterList)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_ForwardEuler_SinCos.xml");

  //std::ofstream ftmp("PL.txt");
  //pList->print(ftmp);
  //ftmp.close();

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  RCP<SinCosModel<double> > model =
    Teuchos::rcp(new SinCosModel<double> (scm_pl));

  RCP<ParameterList> tempusPL  = sublist(pList, "Tempus", true);

  // Test constructor IntegratorBasic(tempusPL, model)
  {
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(tempusPL, model);

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Demo Stepper", true);
    RCP<ParameterList> defaultPL =
      integrator->getStepper()->getDefaultParameters();
    TEST_ASSERT(haveSameValues(*stepperPL, *defaultPL, true))
  }

  // Test constructor IntegratorBasic(model, stepperType)
  {
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(model, "Forward Euler");

    RCP<ParameterList> stepperPL = sublist(tempusPL, "Demo Stepper", true);
    RCP<ParameterList> defaultPL =
      integrator->getStepper()->getDefaultParameters();

    TEST_ASSERT(haveSameValues(*stepperPL, *defaultPL, true))
  }
}
#endif // TEST_PARAMETERLIST


#ifdef TEST_CONSTRUCTING_FROM_DEFAULTS
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, ConstructingFromDefaults)
{
  double dt = 0.1;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_ForwardEuler_SinCos.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  //RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
  RCP<SinCosModel<double> > model =
    Teuchos::rcp(new SinCosModel<double>(scm_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperForwardEuler<double> > stepper =
    Teuchos::rcp(new Tempus::StepperForwardEuler<double>(model));

  // Setup TimeStepControl ------------------------------------
  RCP<Tempus::TimeStepControl<double> > timeStepControl =
    Teuchos::rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL = pl->sublist("Demo Integrator")
                           .sublist("Time Step Control");
  timeStepControl->setStepType (tscPL.get<std::string>("Integrator Step Type"));
  timeStepControl->setInitIndex(tscPL.get<int>   ("Initial Time Index"));
  timeStepControl->setInitTime (tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
  timeStepControl->setInitTimeStep(dt);
  timeStepControl->initialize();

  // Setup initial condition SolutionState --------------------
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =
    stepper->getModel()->getNominalValues();
  RCP<Thyra::VectorBase<double> > icSolution =
    Teuchos::rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  RCP<Tempus::SolutionState<double> > icState =
      Teuchos::rcp(new Tempus::SolutionState<double>(icSolution));
  icState->setTime    (timeStepControl->getInitTime());
  icState->setIndex   (timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  RCP<Tempus::SolutionHistory<double> > solutionHistory =
    Teuchos::rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->setName("Forward States");
  solutionHistory->setStorageType(Tempus::STORAGE_TYPE_STATIC);
  solutionHistory->setStorageLimit(2);
  solutionHistory->addState(icState);

  // Setup Integrator -----------------------------------------
  RCP<Tempus::IntegratorBasic<double> > integrator =
    Tempus::integratorBasic<double>();
  integrator->setStepperWStepper(stepper);
  integrator->setTimeStepControl(timeStepControl);
  integrator->setSolutionHistory(solutionHistory);
  //integrator->setObserver(...);
  integrator->initialize();


  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)


  // Test if at 'Final Time'
  double time = integrator->getTime();
  double timeFinal =pl->sublist("Demo Integrator")
     .sublist("Time Step Control").get<double>("Final Time");
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Time-integrated solution and the exact solution
  RCP<Thyra::VectorBase<double> > x = integrator->getX();
  RCP<const Thyra::VectorBase<double> > x_exact =
    model->getExactSolution(time).get_x();

  // Calculate the error
  RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
  Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));

  // Check the order and intercept
  std::cout << "  Stepper = ForwardEuler" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Exact solution   : " << get_ele(*(x_exact), 0) << "   "
                                       << get_ele(*(x_exact), 1) << std::endl;
  std::cout << "  Computed solution: " << get_ele(*(x      ), 0) << "   "
                                       << get_ele(*(x      ), 1) << std::endl;
  std::cout << "  Difference       : " << get_ele(*(xdiff  ), 0) << "   "
                                       << get_ele(*(xdiff  ), 1) << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 0.882508, 1.0e-4 );
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1), 0.570790, 1.0e-4 );
}
#endif // TEST_CONSTRUCTING_FROM_DEFAULTS


#ifdef TEST_SINCOS
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, SinCos)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 7;
  double dt = 0.2;
  double time = 0.0;
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_ForwardEuler_SinCos.xml");

    //std::ofstream ftmp("PL.txt");
    //pList->print(ftmp);
    //ftmp.close();

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    //RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
    RCP<SinCosModel<double> > model =
      Teuchos::rcp(new SinCosModel<double> (scm_pl));

    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::integratorBasic<double>(pl, model);

    // Initial Conditions
    // During the Integrator construction, the initial SolutionState
    // is set by default to model->getNominalVales().get_x().  However,
    // the application can set it also by integrator->setInitialState.
    RCP<Thyra::VectorBase<double> > x0 =
      model->getNominalValues().get_x()->clone_v();
    integrator->setInitialState(0.0, x0);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test PhysicsState
    Teuchos::RCP<Tempus::PhysicsState<double> > physicsState =
      integrator->getSolutionHistory()->getCurrentState()->getPhysicsState();
    TEST_EQUALITY(physicsState->getName(), "Tempus::PhysicsState");

    // Test if at 'Final Time'
    time = integrator->getTime();
    double timeFinal = pl->sublist("Demo Integrator")
      .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    RCP<Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::VectorBase<double> > x_exact =
      model->getExactSolution(time).get_x();

    // Plot sample solution and exact solution
    if (n == 0) {
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      writeSolution("Tempus_ForwardEuler_SinCos.dat", solutionHistory);

      RCP<Tempus::SolutionHistory<double> > solnHistExact =
        Teuchos::rcp(new Tempus::SolutionHistory<double>());
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        double time = (*solutionHistory)[i]->getTime();
        RCP<Tempus::SolutionState<double> > state =
          Teuchos::rcp(new Tempus::SolutionState<double>(
            model->getExactSolution(time).get_x(),
            model->getExactSolution(time).get_x_dot()));
        state->setTime((*solutionHistory)[i]->getTime());
        solnHistExact->addState(state);
      }
      writeSolution("Tempus_ForwardEuler_SinCos-Ref.dat", solnHistExact);
    }

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXdot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);
    if (n == nTimeStepSizes-1) {  // Add exact solution last in vector.
      StepSize.push_back(0.0);
      auto solution = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x()),solution.ptr());
      solutions.push_back(solution);
      auto solutionDot = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                  solutionDot.ptr());
      solutionsDot.push_back(solutionDot);
    }
  }

  // Check the order and intercept
  double xSlope = 0.0;
  double xDotSlope = 0.0;
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  double order = stepper->getOrder();
  writeOrderError("Tempus_ForwardEuler_SinCos-Error.dat",
                  stepper, StepSize,
                  solutions,    xErrorNorm,    xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope);

  TEST_FLOATING_EQUALITY( xSlope,               order, 0.01   );
  TEST_FLOATING_EQUALITY( xErrorNorm[0],     0.051123, 1.0e-4 );
  // xDot not yet available for Forward Euler.
  //TEST_FLOATING_EQUALITY( xDotSlope,            order, 0.01   );
  //TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 0.0486418, 1.0e-4 );

  Teuchos::TimeMonitor::summarize();
}
#endif // TEST_SINCOS


#ifdef TEST_VANDERPOL
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, VanDerPol)
{
  RCP<Tempus::IntegratorBasic<double> > integrator;
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
  std::vector<double> StepSize;
  std::vector<double> xErrorNorm;
  std::vector<double> xDotErrorNorm;
  const int nTimeStepSizes = 7;  // 8 for Error plot
  double dt = 0.2;
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_ForwardEuler_VanDerPol.xml");

    // Setup the VanDerPolModel
    RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
    RCP<VanDerPolModel<double> > model =
      Teuchos::rcp(new VanDerPolModel<double>(vdpm_pl));

    // Set the step size
    dt /= 2;
    if (n == nTimeStepSizes-1) dt /= 10.0;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Demo Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    integrator = Tempus::integratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Demo Integrator")
      .sublist("Time Step Control").get<double>("Final Time");
    double tol = 100.0 * std::numeric_limits<double>::epsilon();
    TEST_FLOATING_EQUALITY(time, timeFinal, tol);

    // Store off the final solution and step size
    StepSize.push_back(dt);
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    auto solutionDot = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getXdot()),solutionDot.ptr());
    solutionsDot.push_back(solutionDot);

    // Output finest temporal solution for plotting
    // This only works for ONE MPI process
    if ((n == 0) or (n == nTimeStepSizes-1)) {
      std::string fname = "Tempus_ForwardEuler_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_ForwardEuler_VanDerPol.dat";
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      writeSolution(fname, solutionHistory);
    }
  }

  // Check the order and intercept
  double xSlope = 0.0;
  double xDotSlope = 0.0;
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  double order = stepper->getOrder();
  writeOrderError("Tempus_ForwardEuler_VanDerPol-Error.dat",
                  stepper, StepSize,
                  solutions,    xErrorNorm,    xSlope,
                  solutionsDot, xDotErrorNorm, xDotSlope);

  TEST_FLOATING_EQUALITY( xSlope,            order, 0.10   );
  TEST_FLOATING_EQUALITY( xErrorNorm[0],  0.387476, 1.0e-4 );
  // xDot not yet available for Forward Euler.
  //TEST_FLOATING_EQUALITY( xDotSlope,       1.74898, 0.10   );
  //TEST_FLOATING_EQUALITY( xDotErrorNorm[0], 1.0038, 1.0e-4 );

  // Calculate the error - use the most temporally refined mesh for
  // the reference solution.
  auto ref_solution = solutions[solutions.size()-1];
  std::vector<double> StepSizeCheck;
  for (std::size_t i=0; i < (solutions.size()-1); ++i) {
    auto tmp = solutions[i];
    Thyra::Vp_StV(tmp.ptr(), -1.0, *ref_solution);
    const double L2norm = Thyra::norm_2(*tmp);
    StepSizeCheck.push_back(StepSize[i]);
    xErrorNorm.push_back(L2norm);
  }

  Teuchos::TimeMonitor::summarize();
}
#endif // TEST_VANDERPOL


#ifdef TEST_NUMBER_TIMESTEPS
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ForwardEuler, NumberTimeSteps)
{

  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  //const int nTimeStepSizes = 7;
  //double dt = 0.2;
  //double order = 0.0;

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_ForwardEuler_NumberOfTimeSteps.xml");

    // Setup the VanDerPolModel
    RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
    RCP<VanDerPolModel<double> > model =
      Teuchos::rcp(new VanDerPolModel<double>(vdpm_pl));

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);

    //dt = pl->sublist("Demo Integrator")
    //        .sublist("Time Step Control")
    //        .get<double>("Initial Time Step");
    const int numTimeSteps = pl->sublist("Demo Integrator")
                                .sublist("Time Step Control")
                                .get<int>("Number of Time Steps");
    const std::string integratorStepperType =
      pl->sublist("Demo Integrator")
         .sublist("Time Step Control")
         .get<std::string>("Integrator Step Type");

    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // check that the number of time steps taken is whats is set
    // in the parameter list
    TEST_EQUALITY(numTimeSteps, integrator->getIndex());
}
#endif // TEST_NUMBER_TIMESTEPS


} // namespace Tempus_Test
