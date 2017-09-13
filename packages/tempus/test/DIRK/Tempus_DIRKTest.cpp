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

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

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
  RKMethods.push_back("SDIRK 1 Stage 1st order");
  RKMethods.push_back("SDIRK 2 Stage 2nd order");
  RKMethods.push_back("SDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage 3rd order");
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

    if (RKMethods[m] == "SDIRK 2 Stage 2nd order") {
      // Construct in the same order as default.
      RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
      RCP<ParameterList> solverPL = Teuchos::parameterList();
      *solverPL  = *(sublist(stepperPL, "Default Solver", true));
      tempusPL->sublist("Default Stepper").remove("Default Solver");
      tempusPL->sublist("Default Stepper")
           .set<double>("gamma", 0.2928932188134524);
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

      TEST_ASSERT(haveSameValues(*stepperPL,*defaultPL))
    }

    // Test constructor IntegratorBasic(model, stepperType)
    {
      RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::integratorBasic<double>(model, RKMethods[m]);

      RCP<ParameterList> stepperPL = sublist(tempusPL, "Default Stepper", true);
      RCP<ParameterList> defaultPL =
        integrator->getStepper()->getDefaultParameters();
      defaultPL->remove("Description");

      TEST_ASSERT(haveSameValues(*stepperPL,*defaultPL))
    }
  }
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK, SinCos)
{
  std::vector<std::string> RKMethods;
  RKMethods.push_back("SDIRK 1 Stage 1st order");
  RKMethods.push_back("SDIRK 2 Stage 2nd order");
  RKMethods.push_back("SDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage 3rd order");
  RKMethods.push_back("SDIRK 3 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 5th order");

  std::vector<double> RKMethodErrors;
  RKMethodErrors.push_back(0.0486418);
  RKMethodErrors.push_back(0.000404087);
  RKMethodErrors.push_back(8.91985e-05);
  RKMethodErrors.push_back(1.38785e-05);
  RKMethodErrors.push_back(1.61576e-05);
  RKMethodErrors.push_back(8.46096e-08);
  RKMethodErrors.push_back(1.38629e-08);

  for(std::vector<std::string>::size_type m = 0; m != RKMethods.size(); m++) {

    std::string RKMethod_ = RKMethods[m];
    std::replace(RKMethod_.begin(), RKMethod_.end(), ' ', '_');
    std::replace(RKMethod_.begin(), RKMethod_.end(), '/', '.');
    std::vector<double> StepSize;
    std::vector<double> ErrorNorm;
    const int nTimeStepSizes = 7;
    double dt = 0.2;
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
      if (RKMethods[m] == "SDIRK 2 Stage 2nd order") {
        pl->sublist("Default Stepper").set("gamma", 0.2928932190);
      } else if (RKMethods[m] == "SDIRK 2 Stage 3rd order") {
        pl->sublist("Default Stepper").set("3rd Order A-stable", true);
        pl->sublist("Default Stepper").set("2nd Order L-stable", false);
        pl->sublist("Default Stepper").set("gamma", 0.7886751347);
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

    if (RKMethods[m] == "SDIRK 5 Stage 4th order") {
      StepSize.pop_back();  StepSize.pop_back();
      ErrorNorm.pop_back(); ErrorNorm.pop_back();
    } else if (RKMethods[m] == "SDIRK 5 Stage 5th order") {
      StepSize.pop_back();  StepSize.pop_back();  StepSize.pop_back();
      ErrorNorm.pop_back(); ErrorNorm.pop_back(); ErrorNorm.pop_back();
    }

    // Check the order and intercept
    double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
    std::cout << "  Stepper = " << RKMethods[m] << std::endl;
    std::cout << "  =========================" << std::endl;
    std::cout << "  Expected order: " << order << std::endl;
    std::cout << "  Observed order: " << slope << std::endl;
    std::cout << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY( slope, order, 0.01 );
    TEST_FLOATING_EQUALITY( ErrorNorm[0], RKMethodErrors[m], 1.0e-4 );

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
  const int nTimeStepSizes = 5;  // 8 for error plot
  double dt = 0.20;
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
  TEST_FLOATING_EQUALITY( slope, order, 0.05 );
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
