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

#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ImplicitModel.hpp"
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
TEUCHOS_UNIT_TEST(IMEX_RK, VanDerPol)
{
  std::vector<std::string> stepperTypes;
  stepperTypes.push_back("IMEX RK 1st order");
  stepperTypes.push_back("IMEX RK SSP2"     );
  stepperTypes.push_back("IMEX RK ARS 233"  );
  //stepperTypes.push_back("General IMEX RK"  );

  std::vector<double> stepperOrders;
  stepperOrders.push_back(1.21571);
  stepperOrders.push_back(1.94113);
  stepperOrders.push_back(3.14676);
  //stepperOrders.push_back(1.0);

  std::vector<double> stepperErrors;
  stepperErrors.push_back(0.136124);
  stepperErrors.push_back(0.0269125);
  stepperErrors.push_back(0.0309342);
  //stepperErrors.push_back(1.38785e-05);

  std::vector<double> stepperInitDt;
  stepperInitDt.push_back(0.0125);
  stepperInitDt.push_back(0.05);
  stepperInitDt.push_back(0.05);
  //stepperInitDt.push_back(0.025);

  std::vector<std::string>::size_type m;
  for(m = 0; m != stepperTypes.size(); m++) {

    std::string stepperType = stepperTypes[m];
    std::string stepperName = stepperTypes[m];
    std::replace(stepperName.begin(), stepperName.end(), ' ', '_');
    std::replace(stepperName.begin(), stepperName.end(), '/', '.');

    std::vector<RCP<Thyra::VectorBase<double>>> solutions;
    std::vector<double> StepSize;
    std::vector<double> ErrorNorm;
    const int nTimeStepSizes = 3;
    double dt = stepperInitDt[m];
    double order = 0.0;
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
      pl->sublist("Default Stepper").set("Stepper Type", stepperType);

      // Set the step size
      if (n == nTimeStepSizes-1) dt /= 10.0;
      else dt /= 2;

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
        std::string fname = "Tempus_"+stepperName+"_VanDerPol-Ref.dat";
        if (n == 0) fname = "Tempus_"+stepperName+"_VanDerPol.dat";
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
    std::cout << "  Stepper = " << stepperType << std::endl;
    std::cout << "  =========================" << std::endl;
    std::cout << "  Expected order: " << order << std::endl;
    std::cout << "  Observed order: " << slope << std::endl;
    std::cout << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY( slope, stepperOrders[m], 0.02 );
    TEST_FLOATING_EQUALITY( ErrorNorm[0], stepperErrors[m], 1.0e-4 );

    // Write error data
    {
      std::ofstream ftmp("Tempus_"+stepperName+"_VanDerPol-Error.dat");
      double error0 = 0.8*ErrorNorm[0];
      for (std::size_t n = 0; n < StepSizeCheck.size(); n++) {
        ftmp << StepSizeCheck[n]  << "   " << ErrorNorm[n] << "   "
             << error0*(pow(StepSize[n]/StepSize[0],order)) << std::endl;
      }
      ftmp.close();
    }
  }
  Teuchos::TimeMonitor::summarize();
}


} // namespace Tempus_Test
