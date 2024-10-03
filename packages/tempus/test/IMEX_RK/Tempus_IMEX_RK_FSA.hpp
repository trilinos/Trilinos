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

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultProductVector.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorForwardSensitivity.hpp"
#include "Tempus_WrapperModelEvaluatorPairIMEX_Basic.hpp"

#include "../TestModels/VanDerPolModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ExplicitModel.hpp"
#include "../TestModels/VanDerPol_IMEX_ImplicitModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Test {

using Teuchos::getParametersFromXmlFile;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::sublist;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// ************************************************************
// ************************************************************
void test_vdp_fsa(const bool use_combined_method,
                  const bool use_dfdp_as_tangent, Teuchos::FancyOStream& out,
                  bool& success)
{
  std::vector<std::string> stepperTypes;
  stepperTypes.push_back("IMEX RK 1st order");
  stepperTypes.push_back("IMEX RK SSP2");
  stepperTypes.push_back("IMEX RK ARS 233");
  stepperTypes.push_back("General IMEX RK");

  std::vector<double> stepperOrders;
  std::vector<double> stepperErrors;
  if (use_combined_method) {
    stepperOrders.push_back(1.1198);
    stepperOrders.push_back(1.98931);
    stepperOrders.push_back(2.60509);
    stepperOrders.push_back(1.992);

    stepperErrors.push_back(0.00619674);
    stepperErrors.push_back(0.294989);
    stepperErrors.push_back(0.0062125);
    stepperErrors.push_back(0.142489);
  }
  else {
    stepperOrders.push_back(1.1198);
    stepperOrders.push_back(2.05232);
    stepperOrders.push_back(2.43013);
    stepperOrders.push_back(2.07975);

    stepperErrors.push_back(0.00619674);
    stepperErrors.push_back(0.0534172);
    stepperErrors.push_back(0.00224533);
    stepperErrors.push_back(0.032632);
  }
  std::vector<double> stepperInitDt;
  stepperInitDt.push_back(0.0125);
  stepperInitDt.push_back(0.05);
  stepperInitDt.push_back(0.05);
  stepperInitDt.push_back(0.05);

  Teuchos::RCP<const Teuchos::Comm<int>> comm =
      Teuchos::DefaultComm<int>::getComm();

  std::vector<std::string>::size_type m;
  for (m = 0; m != stepperTypes.size(); m++) {
    std::string stepperType = stepperTypes[m];
    std::string stepperName = stepperTypes[m];
    std::replace(stepperName.begin(), stepperName.end(), ' ', '_');
    std::replace(stepperName.begin(), stepperName.end(), '/', '.');

    std::vector<RCP<Thyra::VectorBase<double>>> solutions;
    std::vector<RCP<Thyra::VectorBase<double>>> sensitivities;
    std::vector<double> StepSize;
    std::vector<double> ErrorNorm;
    const int nTimeStepSizes = 3;  // 6 for error plot
    double dt                = stepperInitDt[m];
    double order             = 0.0;
    for (int n = 0; n < nTimeStepSizes; n++) {
      // Read params from .xml file
      RCP<ParameterList> pList =
          getParametersFromXmlFile("Tempus_IMEX_RK_VanDerPol.xml");

      // Setup the explicit VanDerPol ModelEvaluator
      RCP<ParameterList> vdpmPL = sublist(pList, "VanDerPolModel", true);
      vdpmPL->set("Use DfDp as Tangent", use_dfdp_as_tangent);
      RCP<VanDerPol_IMEX_ExplicitModel<double>> explicitModel =
          Teuchos::rcp(new VanDerPol_IMEX_ExplicitModel<double>(vdpmPL));

      // Setup the implicit VanDerPol ModelEvaluator (reuse vdpmPL)
      RCP<VanDerPol_IMEX_ImplicitModel<double>> implicitModel =
          Teuchos::rcp(new VanDerPol_IMEX_ImplicitModel<double>(vdpmPL));

      // Setup the IMEX Pair ModelEvaluator
      RCP<Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>> model =
          Teuchos::rcp(new Tempus::WrapperModelEvaluatorPairIMEX_Basic<double>(
              explicitModel, implicitModel));

      // Setup sensitivities
      RCP<ParameterList> pl  = sublist(pList, "Tempus", true);
      ParameterList& sens_pl = pl->sublist("Sensitivities");
      if (use_combined_method)
        sens_pl.set("Sensitivity Method", "Combined");
      else {
        sens_pl.set("Sensitivity Method", "Staggered");
        sens_pl.set("Reuse State Linear Solver", true);
      }
      sens_pl.set("Use DfDp as Tangent", use_dfdp_as_tangent);
      ParameterList& interp_pl = pl->sublist("Default Integrator")
                                     .sublist("Solution History")
                                     .sublist("Interpolator");
      interp_pl.set("Interpolator Type", "Lagrange");
      interp_pl.set("Order", 2);  // All RK methods here are at most 3rd order

      // Set the Stepper
      if (stepperType == "General IMEX RK") {
        // use the appropriate stepper sublist
        pl->sublist("Default Integrator")
            .set("Stepper Name", "General IMEX RK");
      }
      else {
        pl->sublist("Default Stepper").set("Stepper Type", stepperType);
      }

      // Set the step size
      if (n == nTimeStepSizes - 1)
        dt /= 10.0;
      else
        dt /= 2;

      // Setup the Integrator and reset initial time step
      pl->sublist("Default Integrator")
          .sublist("Time Step Control")
          .set("Initial Time Step", dt);
      pl->sublist("Default Integrator")
          .sublist("Time Step Control")
          .remove("Time Step Control Strategy");
      RCP<Tempus::IntegratorForwardSensitivity<double>> integrator =
          Tempus::createIntegratorForwardSensitivity<double>(pl, model);
      order = integrator->getStepper()->getOrder();

      // Integrate to timeMax
      bool integratorStatus = integrator->advanceTime();
      TEST_ASSERT(integratorStatus)

      // Test if at 'Final Time'
      double time      = integrator->getTime();
      double timeFinal = pl->sublist("Default Integrator")
                             .sublist("Time Step Control")
                             .get<double>("Final Time");
      double tol = 100.0 * std::numeric_limits<double>::epsilon();
      TEST_FLOATING_EQUALITY(time, timeFinal, tol);

      // Store off the final solution and step size
      auto solution    = Thyra::createMember(model->get_x_space());
      auto sensitivity = Thyra::createMember(model->get_x_space());
      Thyra::copy(*(integrator->getX()), solution.ptr());
      Thyra::copy(*(integrator->getDxDp()->col(0)), sensitivity.ptr());
      solutions.push_back(solution);
      sensitivities.push_back(sensitivity);
      StepSize.push_back(dt);

      // Output finest temporal solution for plotting
      if (comm->getRank() == 0 && ((n == 0) || (n == nTimeStepSizes - 1))) {
        typedef Thyra::DefaultMultiVectorProductVector<double> DMVPV;

        std::string fname = "Tempus_" + stepperName + "_VanDerPol_Sens-Ref.dat";
        if (n == 0) fname = "Tempus_" + stepperName + "_VanDerPol_Sens.dat";
        std::ofstream ftmp(fname);
        RCP<const SolutionHistory<double>> solutionHistory =
            integrator->getSolutionHistory();
        int nStates = solutionHistory->getNumStates();
        for (int i = 0; i < nStates; i++) {
          RCP<const SolutionState<double>> solutionState =
              (*solutionHistory)[i];
          RCP<const DMVPV> x_prod =
              Teuchos::rcp_dynamic_cast<const DMVPV>(solutionState->getX());
          RCP<const Thyra::VectorBase<double>> x =
              x_prod->getMultiVector()->col(0);
          RCP<const Thyra::VectorBase<double>> dxdp =
              x_prod->getMultiVector()->col(1);
          double ttime = solutionState->getTime();
          ftmp << std::fixed << std::setprecision(7) << ttime << "   "
               << std::setw(11) << get_ele(*x, 0) << "   " << std::setw(11)
               << get_ele(*x, 1) << "   " << std::setw(11) << get_ele(*dxdp, 0)
               << "   " << std::setw(11) << get_ele(*dxdp, 1) << std::endl;
        }
        ftmp.close();
      }
    }

    // Calculate the error - use the most temporally refined mesh for
    // the reference solution.
    auto ref_solution    = solutions[solutions.size() - 1];
    auto ref_sensitivity = sensitivities[solutions.size() - 1];
    std::vector<double> StepSizeCheck;
    for (std::size_t i = 0; i < (solutions.size() - 1); ++i) {
      auto sol = solutions[i];
      auto sen = sensitivities[i];
      Thyra::Vp_StV(sol.ptr(), -1.0, *ref_solution);
      Thyra::Vp_StV(sen.ptr(), -1.0, *ref_sensitivity);
      const double L2norm_sol = Thyra::norm_2(*sol);
      const double L2norm_sen = Thyra::norm_2(*sen);
      const double L2norm =
          std::sqrt(L2norm_sol * L2norm_sol + L2norm_sen * L2norm_sen);
      StepSizeCheck.push_back(StepSize[i]);
      ErrorNorm.push_back(L2norm);

      out << " n = " << i << " dt = " << StepSize[i] << " error = " << L2norm
          << std::endl;
    }

    // Check the order and intercept
    double slope =
        computeLinearRegressionLogLog<double>(StepSizeCheck, ErrorNorm);
    out << "  Stepper = " << stepperType << std::endl;
    out << "  =========================" << std::endl;
    out << "  Expected order: " << order << std::endl;
    out << "  Observed order: " << slope << std::endl;
    out << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY(slope, stepperOrders[m], 0.02);
    TEST_FLOATING_EQUALITY(ErrorNorm[0], stepperErrors[m], 1.0e-4);

    // Write error data
    {
      std::ofstream ftmp("Tempus_" + stepperName + "_VanDerPol_Sens-Error.dat");
      double error0 = 0.8 * ErrorNorm[0];
      for (std::size_t n = 0; n < StepSizeCheck.size(); n++) {
        ftmp << StepSizeCheck[n] << "   " << ErrorNorm[n] << "   "
             << error0 * (pow(StepSize[n] / StepSize[0], order)) << std::endl;
      }
      ftmp.close();
    }
  }
  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
