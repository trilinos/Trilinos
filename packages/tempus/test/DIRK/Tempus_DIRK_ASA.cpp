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

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorAdjointSensitivity.hpp"

#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultProductVector.hpp"

#include "../TestModels/SinCosModel.hpp"
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
TEUCHOS_UNIT_TEST(DIRK, SinCos_ASA)
{
  std::vector<std::string> RKMethods;
  RKMethods.push_back("General DIRK");
  RKMethods.push_back("RK Backward Euler");
  RKMethods.push_back("DIRK 1 Stage Theta Method");
  RKMethods.push_back("RK Implicit 1 Stage 1st order Radau IA");
  RKMethods.push_back("SDIRK 2 Stage 2nd order");
  RKMethods.push_back("RK Implicit 2 Stage 2nd order Lobatto IIIB");
  RKMethods.push_back("SDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage 3rd order");
  RKMethods.push_back("EDIRK 2 Stage Theta Method");
  RKMethods.push_back("SDIRK 3 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 4th order");
  RKMethods.push_back("SDIRK 5 Stage 5th order");

  std::vector<double> RKMethodErrors;
  RKMethodErrors.push_back(8.48235e-05);
  RKMethodErrors.push_back(0.0383339);
  RKMethodErrors.push_back(0.000221028);
  RKMethodErrors.push_back(0.0428449);
  RKMethodErrors.push_back(8.48235e-05);
  RKMethodErrors.push_back(0.000297933);
  RKMethodErrors.push_back(4.87848e-06);
  RKMethodErrors.push_back(7.30827e-07);
  RKMethodErrors.push_back(0.000272997);
  RKMethodErrors.push_back(3.10132e-07);
  RKMethodErrors.push_back(7.56838e-10);
  RKMethodErrors.push_back(1.32374e-10);

  Teuchos::RCP<const Teuchos::Comm<int> > comm =
      Teuchos::DefaultComm<int>::getComm();

  for (std::vector<std::string>::size_type m = 0; m != RKMethods.size(); m++) {
    std::string RKMethod_ = RKMethods[m];
    std::replace(RKMethod_.begin(), RKMethod_.end(), ' ', '_');
    std::replace(RKMethod_.begin(), RKMethod_.end(), '/', '.');
    std::vector<double> StepSize;
    std::vector<double> ErrorNorm;
    const int nTimeStepSizes = 2;  // 7 for error plots
    double dt                = 0.05;
    double order             = 0.0;
    for (int n = 0; n < nTimeStepSizes; n++) {
      // Read params from .xml file
      RCP<ParameterList> pList =
          getParametersFromXmlFile("Tempus_DIRK_SinCos.xml");

      // Setup the SinCosModel
      RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
      RCP<SinCosModel<double> > model =
          Teuchos::rcp(new SinCosModel<double>(scm_pl));

      // Set the Stepper
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);
      pl->sublist("Default Stepper").set("Stepper Type", RKMethods[m]);
      if (RKMethods[m] == "DIRK 1 Stage Theta Method" ||
          RKMethods[m] == "EDIRK 2 Stage Theta Method") {
        pl->sublist("Default Stepper").set<double>("theta", 0.5);
      }
      else if (RKMethods[m] == "SDIRK 2 Stage 2nd order") {
        pl->sublist("Default Stepper").set("gamma", 0.2928932188134524);
      }
      else if (RKMethods[m] == "SDIRK 2 Stage 3rd order") {
        pl->sublist("Default Stepper")
            .set<std::string>("Gamma Type", "3rd Order A-stable");
      }

      dt /= 2;

      // Setup sensitivities
      ParameterList& sens_pl = pl->sublist("Sensitivities");
      sens_pl.set("Mass Matrix Is Identity", true);  // Necessary for explicit
      ParameterList& interp_pl = pl->sublist("Default Integrator")
                                     .sublist("Solution History")
                                     .sublist("Interpolator");
      interp_pl.set("Interpolator Type", "Lagrange");
      interp_pl.set("Order", 4);  // All RK methods here are at most 5th order

      // Setup the Integrator and reset initial time step
      pl->sublist("Default Integrator")
          .sublist("Time Step Control")
          .set("Initial Time Step", dt);
      RCP<Tempus::IntegratorAdjointSensitivity<double> > integrator =
          Tempus::createIntegratorAdjointSensitivity<double>(pl, model);
      order = integrator->getStepper()->getOrder();

      // Initial Conditions
      // During the Integrator construction, the initial SolutionState
      // is set by default to model->getNominalVales().get_x().  However,
      // the application can set it also by
      // integrator->initializeSolutionHistory.
      RCP<Thyra::VectorBase<double> > x0 =
          model->getNominalValues().get_x()->clone_v();
      const int num_param = model->get_p_space(0)->dim();
      RCP<Thyra::MultiVectorBase<double> > DxDp0 =
          Thyra::createMembers(model->get_x_space(), num_param);
      for (int i = 0; i < num_param; ++i)
        Thyra::assign(DxDp0->col(i).ptr(),
                      *(model->getExactSensSolution(i, 0.0).get_x()));
      integrator->initializeSolutionHistory(0.0, x0, Teuchos::null,
                                            Teuchos::null, DxDp0, Teuchos::null,
                                            Teuchos::null);

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

      // Time-integrated solution and the exact solution along with
      // sensitivities (relying on response g(x) = x).  Note we must transpose
      // dg/dp since the integrator returns it in gradient form.
      RCP<const Thyra::VectorBase<double> > x         = integrator->getX();
      RCP<const Thyra::MultiVectorBase<double> > DgDp = integrator->getDgDp();
      RCP<Thyra::MultiVectorBase<double> > DxDp =
          Thyra::createMembers(model->get_x_space(), num_param);
      {
        Thyra::ConstDetachedMultiVectorView<double> dgdp_view(*DgDp);
        Thyra::DetachedMultiVectorView<double> dxdp_view(*DxDp);
        const int num_g = DgDp->domain()->dim();
        for (int i = 0; i < num_g; ++i)
          for (int j = 0; j < num_param; ++j) dxdp_view(i, j) = dgdp_view(j, i);
      }
      RCP<const Thyra::VectorBase<double> > x_exact =
          model->getExactSolution(time).get_x();
      RCP<Thyra::MultiVectorBase<double> > DxDp_exact =
          Thyra::createMembers(model->get_x_space(), num_param);
      for (int i = 0; i < num_param; ++i)
        Thyra::assign(DxDp_exact->col(i).ptr(),
                      *(model->getExactSensSolution(i, time).get_x()));

      // Plot sample solution and exact solution
      if (comm->getRank() == 0 && n == nTimeStepSizes - 1) {
        typedef Thyra::DefaultProductVector<double> DPV;
        typedef Thyra::DefaultMultiVectorProductVector<double> DMVPV;

        std::ofstream ftmp("Tempus_" + RKMethod_ + "_SinCos_AdjSens.dat");
        RCP<const SolutionHistory<double> > solutionHistory =
            integrator->getSolutionHistory();
        for (int i = 0; i < solutionHistory->getNumStates(); i++) {
          RCP<const SolutionState<double> > solutionState =
              (*solutionHistory)[i];
          const double time_i = solutionState->getTime();
          RCP<const DPV> x_prod_plot =
              Teuchos::rcp_dynamic_cast<const DPV>(solutionState->getX());
          RCP<const Thyra::VectorBase<double> > x_plot =
              x_prod_plot->getVectorBlock(0);
          RCP<const DMVPV> adjoint_prod_plot =
              Teuchos::rcp_dynamic_cast<const DMVPV>(
                  x_prod_plot->getVectorBlock(1));
          RCP<const Thyra::MultiVectorBase<double> > adjoint_plot =
              adjoint_prod_plot->getMultiVector();
          RCP<const Thyra::VectorBase<double> > x_exact_plot =
              model->getExactSolution(time_i).get_x();
          ftmp << std::fixed << std::setprecision(7) << time_i << std::setw(11)
               << get_ele(*(x_plot), 0) << std::setw(11)
               << get_ele(*(x_plot), 1) << std::setw(11)
               << get_ele(*(adjoint_plot->col(0)), 0) << std::setw(11)
               << get_ele(*(adjoint_plot->col(0)), 1) << std::setw(11)
               << get_ele(*(adjoint_plot->col(1)), 0) << std::setw(11)
               << get_ele(*(adjoint_plot->col(1)), 1) << std::setw(11)
               << get_ele(*(x_exact_plot), 0) << std::setw(11)
               << get_ele(*(x_exact_plot), 1) << std::endl;
        }
        ftmp.close();
      }

      // Calculate the error
      RCP<Thyra::VectorBase<double> > xdiff         = x->clone_v();
      RCP<Thyra::MultiVectorBase<double> > DxDpdiff = DxDp->clone_mv();
      Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
      Thyra::V_VmV(DxDpdiff.ptr(), *DxDp_exact, *DxDp);
      StepSize.push_back(dt);
      double L2norm = Thyra::norm_2(*xdiff);
      L2norm *= L2norm;
      Teuchos::Array<double> L2norm_DxDp(num_param);
      Thyra::norms_2(*DxDpdiff, L2norm_DxDp());
      for (int i = 0; i < num_param; ++i)
        L2norm += L2norm_DxDp[i] * L2norm_DxDp[i];
      L2norm = std::sqrt(L2norm);
      ErrorNorm.push_back(L2norm);

      // out << " n = " << n << " dt = " << dt << " error = " << L2norm
      //     << std::endl;
    }

    if (comm->getRank() == 0) {
      std::ofstream ftmp("Tempus_" + RKMethod_ + "_SinCos_AdjSens-Error.dat");
      double error0 = 0.8 * ErrorNorm[0];
      for (int n = 0; n < (int)StepSize.size(); n++) {
        ftmp << StepSize[n] << "   " << ErrorNorm[n] << "   "
             << error0 * (pow(StepSize[n] / StepSize[0], order)) << std::endl;
      }
      ftmp.close();
    }

    // if (RKMethods[m] == "SDIRK 5 Stage 4th order") {
    //   StepSize.pop_back();  StepSize.pop_back();
    //   ErrorNorm.pop_back(); ErrorNorm.pop_back();
    // } else if (RKMethods[m] == "SDIRK 5 Stage 5th order") {
    //   StepSize.pop_back();  StepSize.pop_back();  StepSize.pop_back();
    //   ErrorNorm.pop_back(); ErrorNorm.pop_back(); ErrorNorm.pop_back();
    // }

    // Check the order and intercept
    double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
    out << "  Stepper = " << RKMethods[m] << std::endl;
    out << "  =========================" << std::endl;
    out << "  Expected order: " << order << std::endl;
    out << "  Observed order: " << slope << std::endl;
    out << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY(slope, order, 0.03);
    TEST_FLOATING_EQUALITY(ErrorNorm[0], RKMethodErrors[m], 5.0e-4);
  }
  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
