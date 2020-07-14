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
#include "Teuchos_DefaultComm.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_StepperFactory.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ExplicitRK, ParameterList)
{
  std::vector<std::string> RKMethods;
  RKMethods.push_back("General ERK");
  RKMethods.push_back("RK Forward Euler");
  RKMethods.push_back("RK Explicit 4 Stage");
  RKMethods.push_back("RK Explicit 3/8 Rule");
  RKMethods.push_back("RK Explicit 4 Stage 3rd order by Runge");
  RKMethods.push_back("RK Explicit 5 Stage 3rd order by Kinnmark and Gray");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order TVD");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order by Heun");
  RKMethods.push_back("RK Explicit Midpoint");
  RKMethods.push_back("RK Explicit Trapezoidal");
  RKMethods.push_back("Heuns Method");
  RKMethods.push_back("Bogacki-Shampine 3(2) Pair");
  RKMethods.push_back("Merson 4(5) Pair");

  for(std::vector<std::string>::size_type m = 0; m != RKMethods.size(); m++) {

    std::string RKMethod = RKMethods[m];
    std::replace(RKMethod.begin(), RKMethod.end(), ' ', '_');
    std::replace(RKMethod.begin(), RKMethod.end(), '/', '.');

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_ExplicitRK_SinCos.xml");

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    auto model = rcp(new SinCosModel<double>(scm_pl));

    // Set the Stepper
    RCP<ParameterList> tempusPL  = sublist(pList, "Tempus", true);
    if (RKMethods[m] == "General ERK") {
      tempusPL->sublist("Demo Integrator").set("Stepper Name", "Demo Stepper 2");
    } else {
      tempusPL->sublist("Demo Stepper").set("Stepper Type", RKMethods[m]);
    }

    // Set IC consistency to default value.
    tempusPL->sublist("Demo Stepper")
                 .set("Initial Condition Consistency", "None");

    // Test constructor IntegratorBasic(tempusPL, model)
    {
      RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::integratorBasic<double>(tempusPL, model);

      RCP<ParameterList> stepperPL = sublist(tempusPL, "Demo Stepper", true);
      if (RKMethods[m] == "General ERK")
        stepperPL = sublist(tempusPL, "Demo Stepper 2", true);
      RCP<ParameterList> defaultPL =
        Teuchos::rcp_const_cast<Teuchos::ParameterList>(
          integrator->getStepper()->getValidParameters());
      defaultPL->remove("Description");

      bool pass = haveSameValues(*stepperPL, *defaultPL, true);
      if (!pass) {
        std::cout << std::endl;
        std::cout << "stepperPL -------------- \n" << *stepperPL << std::endl;
        std::cout << "defaultPL -------------- \n" << *defaultPL << std::endl;
      }
      TEST_ASSERT(pass)
    }

    // Test constructor IntegratorBasic(model, stepperType)
    {
      RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::integratorBasic<double>(model, RKMethods[m]);

      RCP<ParameterList> stepperPL = sublist(tempusPL, "Demo Stepper", true);
      if (RKMethods[m] == "General ERK")
        stepperPL = sublist(tempusPL, "Demo Stepper 2", true);
      RCP<ParameterList> defaultPL =
        Teuchos::rcp_const_cast<Teuchos::ParameterList>(
          integrator->getStepper()->getValidParameters());
      defaultPL->remove("Description");

      bool pass = haveSameValues(*stepperPL, *defaultPL, true);
      if (!pass) {
        std::cout << std::endl;
        std::cout << "stepperPL -------------- \n" << *stepperPL << std::endl;
        std::cout << "defaultPL -------------- \n" << *defaultPL << std::endl;
      }
      TEST_ASSERT(pass)
    }
  }
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ExplicitRK, ConstructingFromDefaults)
{
  double dt = 0.1;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_ExplicitRK_SinCos.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  //RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
  auto model = rcp(new SinCosModel<double>(scm_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperFactory<double> > sf =
    Teuchos::rcp(new Tempus::StepperFactory<double>());
  RCP<Tempus::Stepper<double> > stepper =
    sf->createStepper("RK Explicit 4 Stage");
  stepper->setModel(model);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
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
  auto icSolution = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime    (timeStepControl->getInitTime());
  icState->setIndex   (timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);
  icState->setOrder   (stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
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
  std::cout << "  Stepper = RK Explicit 4 Stage" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Exact solution   : " << get_ele(*(x_exact), 0) << "   "
                                       << get_ele(*(x_exact), 1) << std::endl;
  std::cout << "  Computed solution: " << get_ele(*(x      ), 0) << "   "
                                       << get_ele(*(x      ), 1) << std::endl;
  std::cout << "  Difference       : " << get_ele(*(xdiff  ), 0) << "   "
                                       << get_ele(*(xdiff  ), 1) << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 0.841470, 1.0e-4 );
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1), 0.540303, 1.0e-4 );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ExplicitRK, SinCos)
{
  std::vector<std::string> RKMethods;
  RKMethods.push_back("General ERK");
  RKMethods.push_back("RK Forward Euler");
  RKMethods.push_back("RK Explicit 4 Stage");
  RKMethods.push_back("RK Explicit 3/8 Rule");
  RKMethods.push_back("RK Explicit 4 Stage 3rd order by Runge");
  RKMethods.push_back("RK Explicit 5 Stage 3rd order by Kinnmark and Gray");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order TVD");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order by Heun");
  RKMethods.push_back("RK Explicit Midpoint");
  RKMethods.push_back("RK Explicit Trapezoidal");
  RKMethods.push_back("Heuns Method");
  RKMethods.push_back("Bogacki-Shampine 3(2) Pair");
  RKMethods.push_back("Merson 4(5) Pair"); // slope = 3.87816
  RKMethods.push_back("General ERK Embedded");
  RKMethods.push_back("SSPERK22");
  RKMethods.push_back("SSPERK33");
  RKMethods.push_back("SSPERK54");  // slope = 3.94129
  RKMethods.push_back("RK2"); 

  std::vector<double> RKMethodErrors;
  RKMethodErrors.push_back(8.33251e-07);
  RKMethodErrors.push_back(0.051123);
  RKMethodErrors.push_back(8.33251e-07);
  RKMethodErrors.push_back(8.33251e-07);
  RKMethodErrors.push_back(4.16897e-05);
  RKMethodErrors.push_back(8.32108e-06);
  RKMethodErrors.push_back(4.16603e-05);
  RKMethodErrors.push_back(4.16603e-05);
  RKMethodErrors.push_back(4.16603e-05);
  RKMethodErrors.push_back(0.00166645);
  RKMethodErrors.push_back(0.00166645);
  RKMethodErrors.push_back(0.00166645);
  RKMethodErrors.push_back(4.16603e-05);
  RKMethodErrors.push_back(1.39383e-07);
  RKMethodErrors.push_back(4.16603e-05);
  RKMethodErrors.push_back(0.00166645);
  RKMethodErrors.push_back(4.16603e-05);
  RKMethodErrors.push_back(3.85613e-07); // SSPERK54
  RKMethodErrors.push_back(0.00166644); 


  TEST_ASSERT(RKMethods.size() == RKMethodErrors.size() );

  for(std::vector<std::string>::size_type m = 0; m != RKMethods.size(); m++) {

    std::string RKMethod = RKMethods[m];
    std::replace(RKMethod.begin(), RKMethod.end(), ' ', '_');
    std::replace(RKMethod.begin(), RKMethod.end(), '/', '.');

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
        getParametersFromXmlFile("Tempus_ExplicitRK_SinCos.xml");

      // Setup the SinCosModel
      RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
      //RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
      auto model = rcp(new SinCosModel<double>(scm_pl));

      // Set the Stepper
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);
      if (RKMethods[m] == "General ERK") {
         pl->sublist("Demo Integrator").set("Stepper Name", "Demo Stepper 2");
      } else if (RKMethods[m] == "General ERK Embedded"){
         pl->sublist("Demo Integrator").set("Stepper Name", "General ERK Embedded Stepper");
      } else {
         pl->sublist("Demo Stepper").set("Stepper Type", RKMethods[m]);
      }


      dt /= 2;

      // Setup the Integrator and reset initial time step
      pl->sublist("Demo Integrator")
        .sublist("Time Step Control").set("Initial Time Step", dt);
      integrator = Tempus::integratorBasic<double>(pl, model);

      // Initial Conditions
      // During the Integrator construction, the initial SolutionState
      // is set by default to model->getNominalVales().get_x().  However,
      // the application can set it also by integrator->initializeSolutionHistory.
      RCP<Thyra::VectorBase<double> > x0 =
        model->getNominalValues().get_x()->clone_v();
      integrator->initializeSolutionHistory(0.0, x0);

      // Integrate to timeMax
      bool integratorStatus = integrator->advanceTime();
      TEST_ASSERT(integratorStatus)

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
        writeSolution("Tempus_"+RKMethod+"_SinCos.dat", solutionHistory);

        auto solnHistExact = rcp(new Tempus::SolutionHistory<double>());
        for (int i=0; i<solutionHistory->getNumStates(); i++) {
          double time_i = (*solutionHistory)[i]->getTime();
          auto state = Tempus::createSolutionStateX(
            rcp_const_cast<Thyra::VectorBase<double> > (
              model->getExactSolution(time_i).get_x()),
            rcp_const_cast<Thyra::VectorBase<double> > (
              model->getExactSolution(time_i).get_x_dot()));
          state->setTime((*solutionHistory)[i]->getTime());
          solnHistExact->addState(state);
        }
        writeSolution("Tempus_"+RKMethod+"_SinCos-Ref.dat", solnHistExact);
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
        auto solutionExact = Thyra::createMember(model->get_x_space());
        Thyra::copy(*(model->getExactSolution(time).get_x()),solutionExact.ptr());
        solutions.push_back(solutionExact);
        auto solutionDotExact = Thyra::createMember(model->get_x_space());
        Thyra::copy(*(model->getExactSolution(time).get_x_dot()),
                    solutionDotExact.ptr());
        solutionsDot.push_back(solutionDotExact);
      }
    }

    // Check the order and intercept
    double xSlope = 0.0;
    double xDotSlope = 0.0;
    RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
    double order = stepper->getOrder();
    writeOrderError("Tempus_"+RKMethod+"_SinCos-Error.dat",
                    stepper, StepSize,
                    solutions,    xErrorNorm,    xSlope,
                    solutionsDot, xDotErrorNorm, xDotSlope);

    double order_tol = 0.01;
    if (RKMethods[m] == "Merson 4(5) Pair") order_tol = 0.04;
    if (RKMethods[m] == "SSPERK54")         order_tol = 0.06;

    TEST_FLOATING_EQUALITY( xSlope,                    order, order_tol );
    TEST_FLOATING_EQUALITY( xErrorNorm[0], RKMethodErrors[m],    1.0e-4 );
    // xDot not yet available for ExplicitRK methods.
    //TEST_FLOATING_EQUALITY( xDotSlope,                 order, 0.01   );
    //TEST_FLOATING_EQUALITY( xDotErrorNorm[0],      0.0486418, 1.0e-4 );

  }
  //Teuchos::TimeMonitor::summarize();
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ExplicitRK, EmbeddedVanDerPol)
{

  std::vector<std::string> IntegratorList;
  IntegratorList.push_back("Embedded_Integrator_PID");
  IntegratorList.push_back("Demo_Integrator");
  IntegratorList.push_back("Embedded_Integrator");
  IntegratorList.push_back("General_Embedded_Integrator");
  IntegratorList.push_back("Embedded_Integrator_PID_General");

  // the embedded solution will test the following:
  // using the starting stepsize routine, this has now decreased
  const int refIstep = 45;

  for(auto integratorChoice : IntegratorList){

     std::cout << "Using Integrator: " << integratorChoice << " !!!" << std::endl;

     // Read params from .xml file
     RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_ExplicitRK_VanDerPol.xml");


     // Setup the VanDerPolModel
     RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
     auto model = rcp(new VanDerPolModel<double>(vdpm_pl));


     // Set the Integrator and Stepper
     RCP<ParameterList> pl = sublist(pList, "Tempus", true);
     pl->set("Integrator Name", integratorChoice);

     // Setup the Integrator
     RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::integratorBasic<double>(pl, model);

     const std::string RKMethod =
        pl->sublist(integratorChoice).get<std::string>("Stepper Name");

     // Integrate to timeMax
     bool integratorStatus = integrator->advanceTime();
     TEST_ASSERT(integratorStatus);

     // Test if at 'Final Time'
     double time = integrator->getTime();
     double timeFinal = pl->sublist(integratorChoice)
        .sublist("Time Step Control").get<double>("Final Time");
     TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);


     // Numerical reference solution at timeFinal (for \epsilon = 0.1)
     RCP<Thyra::VectorBase<double> > x = integrator->getX();
     RCP<Thyra::VectorBase<double> > xref = x->clone_v();
     Thyra::set_ele(0, -1.5484458614405929, xref.ptr());
     Thyra::set_ele(1,  1.0181127316101317, xref.ptr());

     // Calculate the error
     RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
     Thyra::V_StVpStV(xdiff.ptr(), 1.0, *xref, -1.0, *(x));
     const double L2norm = Thyra::norm_2(*xdiff);

     // Test number of steps, failures, and accuracy
     if ((integratorChoice == "Embedded_Integrator_PID") or
           (integratorChoice == "Embedded_Integrator_PID_General")) {

        const double absTol = pl->sublist(integratorChoice).
           sublist("Time Step Control").get<double>("Maximum Absolute Error");
        const double relTol = pl->sublist(integratorChoice).
           sublist("Time Step Control").get<double>("Maximum Relative Error");

        // Should be close to the prescribed tolerance (magnitude)
        TEST_COMPARE(std::log10(L2norm), <, std::log10(absTol));
        TEST_COMPARE(std::log10(L2norm), <, std::log10(relTol));

        // get the number of time steps and number of step failure
        //const int nFailure_c = integrator->getSolutionHistory()->
        //getCurrentState()->getMetaData()->getNFailures();
        const int iStep = integrator->getSolutionHistory()->
           getCurrentState()->getIndex();
        const int nFail = integrator->getSolutionHistory()->
           getCurrentState()->getMetaData()->getNRunningFailures();

        // test for number of steps
        TEST_EQUALITY(iStep, refIstep);
        std::cout << "Tolerance = " << absTol
           << " L2norm = "   << L2norm
           << " iStep = "    << iStep
           << " nFail = "    << nFail << std::endl;
     }

     // Plot sample solution and exact solution
     std::ofstream ftmp("Tempus_"+integratorChoice+RKMethod+"_VDP_Example.dat");
     RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
     int nStates = solutionHistory->getNumStates();
     //RCP<const Thyra::VectorBase<double> > x_exact_plot;
     for (int i=0; i<nStates; i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        double time_i = solutionState->getTime();
        RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
        //x_exact_plot = model->getExactSolution(time_i).get_x();
        ftmp << time_i << "   "
           << Thyra::get_ele(*(x_plot), 0) << "   "
           << Thyra::get_ele(*(x_plot), 1) << "   " << std::endl;
     }
     ftmp.close();
  }

  Teuchos::TimeMonitor::summarize();
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(ExplicitRK, stage_number)
{
  double dt = 0.1;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_ExplicitRK_SinCos.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  //RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
  auto model = rcp(new SinCosModel<double>(scm_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperFactory<double> > sf =
    Teuchos::rcp(new Tempus::StepperFactory<double>());
  RCP<Tempus::Stepper<double> > stepper =
    sf->createStepper("RK Explicit 4 Stage");
  stepper->setModel(model);
  stepper->initialize();

  // Setup TimeStepControl ------------------------------------
  auto timeStepControl = rcp(new Tempus::TimeStepControl<double>());
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
  auto icSolution = rcp_const_cast<Thyra::VectorBase<double> > (inArgsIC.get_x());
  auto icState = Tempus::createSolutionStateX(icSolution);
  icState->setTime    (timeStepControl->getInitTime());
  icState->setIndex   (timeStepControl->getInitIndex());
  icState->setTimeStep(0.0);
  icState->setOrder   (stepper->getOrder());
  icState->setSolutionStatus(Tempus::Status::PASSED);  // ICs are passing.

  // Setup SolutionHistory ------------------------------------
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
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

  // get the RK stepper
  auto erk_stepper = Teuchos::rcp_dynamic_cast<Tempus::StepperExplicitRK<double>  >(stepper,true);

  TEST_EQUALITY( -1 , erk_stepper->getStageNumber());
  const std::vector<int> ref_stageNumber = { 1, 4, 8, 10, 11, -1, 5};
  for(auto stage_number : ref_stageNumber) {
    // set the stage number
    erk_stepper->setStageNumber(stage_number);
    // make sure we are getting the correct stage number
    TEST_EQUALITY( stage_number, erk_stepper->getStageNumber());
  }
}


} // namespace Tempus_Test
