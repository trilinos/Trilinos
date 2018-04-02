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
#include "Tempus_StepperDIRK.hpp"

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


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK, ParameterList)
{
  std::vector<std::string> RKMethods;
  RKMethods.push_back("RK Backward Euler");
  RKMethods.push_back("IRK 1 Stage Theta Method");
  RKMethods.push_back("SDIRK 1 Stage 1st order");
  RKMethods.push_back("SDIRK 2 Stage 2nd order");
  RKMethods.push_back("SDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage Theta Method");
  RKMethods.push_back("SDIRK 3 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 5th order");

  for(std::vector<std::string>::size_type m = 0; m != RKMethods.size(); m++) {

    std::string RKMethod_ = RKMethods[m];

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_DIRK_SinCos.xml");

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    RCP<SinCosModel<double> > model =
      Teuchos::rcp(new SinCosModel<double>(scm_pl));

    RCP<ParameterList> tempusPL  = sublist(pList, "Tempus", true);
    tempusPL->sublist("Default Stepper").set("Stepper Type", RKMethods[m]);

    if (RKMethods[m] == "IRK 1 Stage Theta Method" ||
        RKMethods[m] == "EDIRK 2 Stage Theta Method") {
      // Construct in the same order as default.
      RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
      RCP<ParameterList> solverPL = Teuchos::parameterList();
      *solverPL  = *(sublist(stepperPL, "Default Solver", true));
      tempusPL->sublist("Default Stepper").remove("Default Solver");
      tempusPL->sublist("Default Stepper")
           .set<double>("theta", 0.5);
      tempusPL->sublist("Default Stepper").remove("Zero Initial Guess");
      tempusPL->sublist("Default Stepper").set<bool>("Zero Initial Guess", 0);
      tempusPL->sublist("Default Stepper").set("Default Solver", *solverPL);
    } else if (RKMethods[m] == "SDIRK 2 Stage 2nd order") {
      // Construct in the same order as default.
      RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
      RCP<ParameterList> solverPL = Teuchos::parameterList();
      *solverPL  = *(sublist(stepperPL, "Default Solver", true));
      tempusPL->sublist("Default Stepper").remove("Default Solver");
      tempusPL->sublist("Default Stepper")
           .set<double>("gamma", 0.2928932188134524);
      tempusPL->sublist("Default Stepper").remove("Zero Initial Guess");
      tempusPL->sublist("Default Stepper").set<bool>("Zero Initial Guess", 0);
      tempusPL->sublist("Default Stepper").set("Default Solver", *solverPL);
    } else if (RKMethods[m] == "SDIRK 2 Stage 3rd order") {
      // Construct in the same order as default.
      RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
      RCP<ParameterList> solverPL = Teuchos::parameterList();
      *solverPL  = *(sublist(stepperPL, "Default Solver", true));
      tempusPL->sublist("Default Stepper").remove("Default Solver");
      tempusPL->sublist("Default Stepper").set("3rd Order A-stable", true);
      tempusPL->sublist("Default Stepper").set("2nd Order L-stable", false);
      tempusPL->sublist("Default Stepper")
           .set<double>("gamma", 0.7886751345948128);
      tempusPL->sublist("Default Stepper").remove("Zero Initial Guess");
      tempusPL->sublist("Default Stepper").set<bool>("Zero Initial Guess", 0);
      tempusPL->sublist("Default Stepper").set("Default Solver", *solverPL);
    }

    // Test constructor IntegratorBasic(tempusPL, model)
    {
      RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::integratorBasic<double>(tempusPL, model);

      RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
      RCP<ParameterList> defaultPL =
        integrator->getStepper()->getDefaultParameters();
      defaultPL->remove("Description");

      TEST_ASSERT(haveSameValues(*stepperPL, *defaultPL, true))
    }

    // Test constructor IntegratorBasic(model, stepperType)
    {
      RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::integratorBasic<double>(model, RKMethods[m]);

      RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
      RCP<ParameterList> defaultPL =
        integrator->getStepper()->getDefaultParameters();
      defaultPL->remove("Description");

    //std::cout << std::endl;
    //std::cout << "stepperPL ----------------- \n" << *stepperPL << std::endl;
    //std::cout << "defaultPL ----------------- \n" << *defaultPL << std::endl;
      TEST_ASSERT(haveSameValues(*stepperPL, *defaultPL, true))
    }
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK, ConstructingFromDefaults)
{
  double dt = 0.025;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_DIRK_SinCos.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  //RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
  RCP<SinCosModel<double> > model =
    Teuchos::rcp(new SinCosModel<double>(scm_pl));

  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperDIRK<double> > stepper =
    Teuchos::rcp(new Tempus::StepperDIRK<double>(model));

  // Setup TimeStepControl ------------------------------------
  RCP<Tempus::TimeStepControl<double> > timeStepControl =
    Teuchos::rcp(new Tempus::TimeStepControl<double>());
  ParameterList tscPL = pl->sublist("Default Integrator")
                           .sublist("Time Step Control");
  timeStepControl->setStepType (tscPL.get<std::string>("Integrator Step Type"));
  timeStepControl->setInitIndex(tscPL.get<int>   ("Initial Time Index"));
  timeStepControl->setInitTime (tscPL.get<double>("Initial Time"));
  timeStepControl->setFinalTime(tscPL.get<double>("Final Time"));
  timeStepControl->setInitTimeStep(dt);

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
  icState->setOrder   (stepper->getOrder());
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
  double timeFinal =pl->sublist("Default Integrator")
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
  std::cout << "  Stepper = SDIRK 2 Stage 2nd order" << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Exact solution   : " << get_ele(*(x_exact), 0) << "   "
                                       << get_ele(*(x_exact), 1) << std::endl;
  std::cout << "  Computed solution: " << get_ele(*(x      ), 0) << "   "
                                       << get_ele(*(x      ), 1) << std::endl;
  std::cout << "  Difference       : " << get_ele(*(xdiff  ), 0) << "   "
                                       << get_ele(*(xdiff  ), 1) << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0), 0.841470, 1.0e-4 );
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1), 0.540304, 1.0e-4 );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK, SinCos)
{
  std::vector<std::string> RKMethods;
  RKMethods.push_back("RK Backward Euler");
  RKMethods.push_back("IRK 1 Stage Theta Method");
  RKMethods.push_back("SDIRK 1 Stage 1st order");
  RKMethods.push_back("SDIRK 2 Stage 2nd order");
  RKMethods.push_back("SDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage Theta Method");
  RKMethods.push_back("SDIRK 3 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 5th order");

  std::vector<double> RKMethodErrors;
  RKMethodErrors.push_back(0.0124201);
  RKMethodErrors.push_back(5.20785e-05);
  RKMethodErrors.push_back(0.0124201);
  RKMethodErrors.push_back(2.52738e-05);
  RKMethodErrors.push_back(1.40223e-06);
  RKMethodErrors.push_back(2.17004e-07);
  RKMethodErrors.push_back(5.20785e-05);
  RKMethodErrors.push_back(6.41463e-08);
  RKMethodErrors.push_back(3.30631e-10);
  RKMethodErrors.push_back(1.35728e-11);

  for(std::vector<std::string>::size_type m = 0; m != RKMethods.size(); m++) {

    std::string RKMethod_ = RKMethods[m];
    std::replace(RKMethod_.begin(), RKMethod_.end(), ' ', '_');
    std::replace(RKMethod_.begin(), RKMethod_.end(), '/', '.');
    std::vector<double> StepSize;
    std::vector<double> ErrorNorm;
    const int nTimeStepSizes = 3; // 7 for error plots
    double dt = 0.05;
    double order = 0.0;
    for (int n=0; n<nTimeStepSizes; n++) {

      // Read params from .xml file
      RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_DIRK_SinCos.xml");

      //std::ofstream ftmp("PL.txt");
      //pList->print(ftmp);
      //ftmp.close();

      // Setup the SinCosModel
      RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
      RCP<SinCosModel<double> > model =
        Teuchos::rcp(new SinCosModel<double>(scm_pl));

      // Set the Stepper
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);
      pl->sublist("Default Stepper").set("Stepper Type", RKMethods[m]);
      if (RKMethods[m] == "IRK 1 Stage Theta Method" ||
          RKMethods[m] == "EDIRK 2 Stage Theta Method") {
        pl->sublist("Default Stepper").set<double>("theta", 0.5);
      } else if (RKMethods[m] == "SDIRK 2 Stage 2nd order") {
        pl->sublist("Default Stepper").set("gamma", 0.2928932188134524);
      } else if (RKMethods[m] == "SDIRK 2 Stage 3rd order") {
        pl->sublist("Default Stepper").set("3rd Order A-stable", true);
        pl->sublist("Default Stepper").set("2nd Order L-stable", false);
        pl->sublist("Default Stepper").set("gamma", 0.7886751345948128);
      }

      dt /= 2;

      // Setup the Integrator and reset initial time step
      pl->sublist("Default Integrator")
         .sublist("Time Step Control").set("Initial Time Step", dt);
      RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::integratorBasic<double>(pl, model);
      order = integrator->getStepper()->getOrder();

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

      // Test if at 'Final Time'
      double time = integrator->getTime();
      double timeFinal = pl->sublist("Default Integrator")
        .sublist("Time Step Control").get<double>("Final Time");
      double tol = 100.0 * std::numeric_limits<double>::epsilon();
      TEST_FLOATING_EQUALITY(time, timeFinal, tol);

      // Time-integrated solution and the exact solution
      RCP<Thyra::VectorBase<double> > x = integrator->getX();
      RCP<const Thyra::VectorBase<double> > x_exact =
        model->getExactSolution(time).get_x();

      // Plot sample solution and exact solution
      if (n == 0) {
        std::ofstream ftmp("Tempus_"+RKMethod_+"_SinCos.dat");
        RCP<const SolutionHistory<double> > solutionHistory =
          integrator->getSolutionHistory();
        int nStates = solutionHistory->getNumStates();
        RCP<const Thyra::VectorBase<double> > x_exact_plot;
        for (int i=0; i<nStates; i++) {
          RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
          double time = solutionState->getTime();
          RCP<const Thyra::VectorBase<double> > x_plot = solutionState->getX();
          x_exact_plot = model->getExactSolution(time).get_x();
          ftmp << time << "   "
               << Thyra::get_ele(*(x_plot), 0) << "   "
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

    std::ofstream ftmp("Tempus_"+RKMethod_+"_SinCos-Error.dat");
    double error0 = 0.8*ErrorNorm[0];
    for (int n=0; n<(int)StepSize.size(); n++) {
      ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
           << error0*(pow(StepSize[n]/StepSize[0],order)) << std::endl;
    }
    ftmp.close();

    //if (RKMethods[m] == "SDIRK 5 Stage 4th order") {
    //  StepSize.pop_back();  StepSize.pop_back();
    //  ErrorNorm.pop_back(); ErrorNorm.pop_back();
    //} else if (RKMethods[m] == "SDIRK 5 Stage 5th order") {
    //  StepSize.pop_back();  StepSize.pop_back();  StepSize.pop_back();
    //  ErrorNorm.pop_back(); ErrorNorm.pop_back(); ErrorNorm.pop_back();
    //}

    // Check the order and intercept
    double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
    std::cout << "  Stepper = " << RKMethods[m] << std::endl;
    std::cout << "  =========================" << std::endl;
    std::cout << "  Expected order: " << order << std::endl;
    std::cout << "  Observed order: " << slope << std::endl;
    std::cout << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY( slope, order, 0.03 );
    TEST_FLOATING_EQUALITY( ErrorNorm[0], RKMethodErrors[m], 5.0e-4 );

  }
  Teuchos::TimeMonitor::summarize();
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK, VanDerPol)
{
  std::string RKMethod_ = "SDIRK 2 Stage 2nd order";
  std::vector<RCP<Thyra::VectorBase<double>>> solutions;
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 3;  // 8 for error plot
  double dt = 0.05;
  double order = 0.0;
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_DIRK_VanDerPol.xml");

    // Setup the VanDerPolModel
    RCP<ParameterList> vdpm_pl = sublist(pList, "VanDerPolModel", true);
    RCP<VanDerPolModel<double> > model =
      Teuchos::rcp(new VanDerPolModel<double>(vdpm_pl));

    // Set the Stepper
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Default Stepper").set("Stepper Type", RKMethod_);
    pl->sublist("Default Stepper").set("gamma", 0.2928932188134524);

    // Set the step size
    dt /= 2;
    if (n == nTimeStepSizes-1) dt /= 10.0;

    // Setup the Integrator and reset initial time step
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::integratorBasic<double>(pl, model);
    order = integrator->getStepper()->getOrder();

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
      .sublist("Time Step Control").get<double>("Final Time");
    double tol = 100.0 * std::numeric_limits<double>::epsilon();
    TEST_FLOATING_EQUALITY(time, timeFinal, tol);

    // Store off the final solution and step size
    auto solution = Thyra::createMember(model->get_x_space());
    Thyra::copy(*(integrator->getX()),solution.ptr());
    solutions.push_back(solution);
    StepSize.push_back(dt);

    // Output finest temporal solution for plotting
    // This only works for ONE MPI process
    if ((n == 0) or (n == nTimeStepSizes-1)) {
      std::string fname = "Tempus_DIRK_VanDerPol-Ref.dat";
      if (n == 0) fname = "Tempus_DIRK_VanDerPol.dat";
      std::ofstream ftmp(fname);
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      int nStates = solutionHistory->getNumStates();
      for (int i=0; i<nStates; i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        RCP<const Thyra::VectorBase<double> > x = solutionState->getX();
        double ttime = solutionState->getTime();
        ftmp << ttime << "   " << get_ele(*x, 0) << "   " << get_ele(*x, 1)
             << std::endl;
      }
      ftmp.close();
    }
  }

  // Calculate the error - use the most temporally refined mesh for
  // the reference solution.
  auto ref_solution = solutions[solutions.size()-1];
  std::vector<double> StepSizeCheck;
  for (std::size_t i=0; i < (solutions.size()-1); ++i) {
    auto tmp = solutions[i];
    Thyra::Vp_StV(tmp.ptr(), -1.0, *ref_solution);
    const double L2norm = Thyra::norm_2(*tmp);
    StepSizeCheck.push_back(StepSize[i]);
    ErrorNorm.push_back(L2norm);
  }

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSizeCheck,ErrorNorm);
  std::cout << "  Stepper = " << RKMethod_ << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Expected order: " << order << std::endl;
  std::cout << "  Observed order: " << slope << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope, order, 0.06 );
  out << "\n\n ** Slope on " << RKMethod_ << " = " << slope
      << "\n" << std::endl;

  // Write error data
  {
    std::ofstream ftmp("Tempus_DIRK_VanDerPol-Error.dat");
    double error0 = 0.8*ErrorNorm[0];
    for (std::size_t n = 0; n < StepSizeCheck.size(); n++) {
      ftmp << StepSizeCheck[n]  << "   " << ErrorNorm[n] << "   "
           << error0*(pow(StepSize[n]/StepSize[0],order)) << std::endl;
    }
    ftmp.close();
  }
  Teuchos::TimeMonitor::summarize();
}


} // namespace Tempus_Test
