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
#include "Tempus_WrapperModelEvaluatorPairIMEX_Basic.hpp"
#include "Tempus_StepperIMEX_RK.hpp"

#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ImplicitModel.hpp"
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
#define TEST_CONSTRUCTING_FROM_DEFAULTS
#define TEST_VANDERPOL


#ifdef TEST_CONSTRUCTING_FROM_DEFAULTS
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK, ConstructingFromDefaults)
{
  double dt = 0.025;

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_IMEX_RK_VanDerPol.xml");
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);

  // Setup the explicit VanDerPol ModelEvaluator
  RCP<ParameterList> vdpmPL = sublist(pList, "VanDerPolModel", true);
  RCP<VanDerPol_IMEX_ExplicitModel<double> > explicitModel =
    Teuchos::rcp(new VanDerPol_IMEX_ExplicitModel<double>(vdpmPL));

  // Setup the implicit VanDerPol ModelEvaluator (reuse vdpmPL)
  RCP<VanDerPol_IMEX_ImplicitModel<double> > implicitModel =
    Teuchos::rcp(new VanDerPol_IMEX_ImplicitModel<double>(vdpmPL));

  // Setup the IMEX Pair ModelEvaluator
  RCP<Tempus::WrapperModelEvaluatorPairIMEX_Basic<double> > model =
      Teuchos::rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
                                             explicitModel, implicitModel));


  // Setup Stepper for field solve ----------------------------
  RCP<Tempus::StepperIMEX_RK<double> > stepper =
    Teuchos::rcp(new Tempus::StepperIMEX_RK<double>(model));

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

  // Check the order and intercept
  std::cout << "  Stepper = " << stepper->description() << std::endl;
  std::cout << "  =========================" << std::endl;
  std::cout << "  Computed solution: " << get_ele(*(x      ), 0) << "   "
                                       << get_ele(*(x      ), 1) << std::endl;
  std::cout << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY(get_ele(*(x), 0),  1.810210, 1.0e-4 );
  TEST_FLOATING_EQUALITY(get_ele(*(x), 1), -0.754602, 1.0e-4 );
}
#endif // TEST_CONSTRUCTING_FROM_DEFAULTS


#ifdef TEST_VANDERPOL
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(IMEX_RK, VanDerPol)
{
  std::vector<std::string> stepperTypes;
  stepperTypes.push_back("IMEX RK 1st order");
  stepperTypes.push_back("IMEX RK SSP2"     );
  stepperTypes.push_back("IMEX RK ARS 233"  );
  stepperTypes.push_back("General IMEX RK"  );

  std::vector<double> stepperOrders;
  stepperOrders.push_back(1.07964);
  stepperOrders.push_back(2.00408);
  stepperOrders.push_back(2.70655);
  stepperOrders.push_back(2.00211);

  std::vector<double> stepperErrors;
  stepperErrors.push_back(0.0046423);
  stepperErrors.push_back(0.0154534);
  stepperErrors.push_back(0.000298908);
  stepperErrors.push_back(0.0071546);

  std::vector<double> stepperInitDt;
  stepperInitDt.push_back(0.0125);
  stepperInitDt.push_back(0.05);
  stepperInitDt.push_back(0.05);
  stepperInitDt.push_back(0.05);

  std::vector<std::string>::size_type m;
  for(m = 0; m != stepperTypes.size(); m++) {

    std::string stepperType = stepperTypes[m];
    std::string stepperName = stepperTypes[m];
    std::replace(stepperName.begin(), stepperName.end(), ' ', '_');
    std::replace(stepperName.begin(), stepperName.end(), '/', '.');

    RCP<Tempus::IntegratorBasic<double> > integrator;
    std::vector<RCP<Thyra::VectorBase<double>>> solutions;
    std::vector<RCP<Thyra::VectorBase<double>>> solutionsDot;
    std::vector<double> StepSize;
    std::vector<double> xErrorNorm;
    std::vector<double> xDotErrorNorm;

    const int nTimeStepSizes = 3;  // 6 for error plot
    double dt = stepperInitDt[m];
    double time = 0.0;
    for (int n=0; n<nTimeStepSizes; n++) {

      // Read params from .xml file
      RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_IMEX_RK_VanDerPol.xml");

      // Setup the explicit VanDerPol ModelEvaluator
      RCP<ParameterList> vdpmPL = sublist(pList, "VanDerPolModel", true);
      RCP<VanDerPol_IMEX_ExplicitModel<double> > explicitModel =
        Teuchos::rcp(new VanDerPol_IMEX_ExplicitModel<double>(vdpmPL));

      // Setup the implicit VanDerPol ModelEvaluator (reuse vdpmPL)
      RCP<VanDerPol_IMEX_ImplicitModel<double> > implicitModel =
        Teuchos::rcp(new VanDerPol_IMEX_ImplicitModel<double>(vdpmPL));

      // Setup the IMEX Pair ModelEvaluator
      RCP<Tempus::WrapperModelEvaluatorPairIMEX_Basic<double> > model =
          Teuchos::rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
                                                 explicitModel, implicitModel));

      // Set the Stepper
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);
      if (stepperType == "General IMEX RK"){
          // use the appropriate stepper sublist
          pl->sublist("Default Integrator").set("Stepper Name", "General IMEX RK");
      }  else {
          pl->sublist("Default Stepper").set("Stepper Type", stepperType);
      }

      // Set the step size
      if (n == nTimeStepSizes-1) dt /= 10.0;
      else dt /= 2;

      // Setup the Integrator and reset initial time step
      pl->sublist("Default Integrator")
         .sublist("Time Step Control").set("Initial Time Step", dt);
      integrator = Tempus::integratorBasic<double>(pl, model);

      // Integrate to timeMax
      bool integratorStatus = integrator->advanceTime();
      TEST_ASSERT(integratorStatus)

      // Test if at 'Final Time'
      time = integrator->getTime();
      double timeFinal =pl->sublist("Default Integrator")
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
        std::string fname = "Tempus_"+stepperName+"_VanDerPol-Ref.dat";
        if (n == 0) fname = "Tempus_"+stepperName+"_VanDerPol.dat";
        RCP<const SolutionHistory<double> > solutionHistory =
          integrator->getSolutionHistory();
        writeSolution(fname, solutionHistory);
      }
    }

    // Check the order and intercept
    double xSlope = 0.0;
    double xDotSlope = 0.0;
    RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
    //double order = stepper->getOrder();
    writeOrderError("Tempus_"+stepperName+"_VanDerPol-Error.dat",
                    stepper, StepSize,
                    solutions,    xErrorNorm,    xSlope,
                    solutionsDot, xDotErrorNorm, xDotSlope);

    TEST_FLOATING_EQUALITY( xSlope,        stepperOrders[m],   0.02 );
    TEST_FLOATING_EQUALITY( xErrorNorm[0], stepperErrors[m], 1.0e-4 );
    // xDot not yet available for IMEX_RK.
    //TEST_FLOATING_EQUALITY( xDotSlope,        1.74898, 0.10 );
    //TEST_FLOATING_EQUALITY( xDotErrorNorm[0],  1.0038, 1.0e-4 );

  }
  //Teuchos::TimeMonitor::summarize();
}
#endif // TEST_VANDERPOL


} // namespace Tempus_Test
