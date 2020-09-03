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

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorAdjointSensitivity.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

#include <fstream>
#include <limits>
#include <sstream>
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
TEUCHOS_UNIT_TEST(BDF2, SinCos_ASA)
{
  std::vector<double> StepSize;
  std::vector<double> ErrorNorm;
  const int nTimeStepSizes = 7;
  double dt = 0.2;
  double order = 0.0;
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();
  Teuchos::RCP<Teuchos::FancyOStream> my_out =
    Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  my_out->setProcRankAndSize(comm->getRank(), comm->getSize());
  my_out->setOutputToRootOnly(0);
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BDF2_SinCos_SA.xml");

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    RCP<SinCosModel<double> > model =
      Teuchos::rcp(new SinCosModel<double>(scm_pl));

    dt /= 2;

    // Setup sensitivities
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    //ParameterList& sens_pl = pl->sublist("Sensitivities");
    ParameterList& interp_pl =
      pl->sublist("Default Integrator").sublist("Solution History")
                                       .sublist("Interpolator");
    interp_pl.set("Interpolator Type", "Lagrange");
    interp_pl.set("Order", 1);

    // Setup the Integrator and reset initial time step
    pl->sublist("Default Integrator")
       .sublist("Time Step Control").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorAdjointSensitivity<double> > integrator =
      Tempus::integratorAdjointSensitivity<double>(pl, model);
    order = integrator->getStepper()->getOrder();

    // Initial Conditions
    double t0 = pl->sublist("Default Integrator")
      .sublist("Time Step Control").get<double>("Initial Time");
    RCP<const Thyra::VectorBase<double> > x0 =
      model->getExactSolution(t0).get_x();
    const int num_param = model->get_p_space(0)->dim();
    RCP<Thyra::MultiVectorBase<double> > DxDp0 =
      Thyra::createMembers(model->get_x_space(), num_param);
    for (int i=0; i<num_param; ++i)
      Thyra::assign(DxDp0->col(i).ptr(),
                    *(model->getExactSensSolution(i, t0).get_x()));
    integrator->initializeSolutionHistory(t0, x0, Teuchos::null, Teuchos::null,
                                DxDp0, Teuchos::null, Teuchos::null);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Default Integrator")
       .sublist("Time Step Control").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution along with
    // sensitivities (relying on response g(x) = x).  Note we must transpose
    // dg/dp since the integrator returns it in gradient form.
    RCP<const Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::MultiVectorBase<double> > DgDp = integrator->getDgDp();
    RCP<Thyra::MultiVectorBase<double> > DxDp =
      Thyra::createMembers(model->get_x_space(), num_param);
    {
      Thyra::ConstDetachedMultiVectorView<double> dgdp_view(*DgDp);
      Thyra::DetachedMultiVectorView<double> dxdp_view(*DxDp);
      const int num_g = DgDp->domain()->dim();
      for (int i=0; i<num_g; ++i)
        for (int j=0; j<num_param; ++j)
          dxdp_view(i,j) = dgdp_view(j,i);
    }
    RCP<const Thyra::VectorBase<double> > x_exact =
      model->getExactSolution(time).get_x();
    RCP<Thyra::MultiVectorBase<double> > DxDp_exact =
      Thyra::createMembers(model->get_x_space(), num_param);
    for (int i=0; i<num_param; ++i)
      Thyra::assign(DxDp_exact->col(i).ptr(),
                    *(model->getExactSensSolution(i, time).get_x()));

    // Plot sample solution and exact solution
    if (comm->getRank() == 0 && n == nTimeStepSizes-1) {
      typedef Thyra::DefaultProductVector<double> DPV;
      typedef Thyra::DefaultMultiVectorProductVector<double> DMVPV;

      std::ofstream ftmp("Tempus_BDF2_SinCos_AdjSens.dat");
      RCP<const SolutionHistory<double> > solutionHistory =
        integrator->getSolutionHistory();
      for (int i=0; i<solutionHistory->getNumStates(); i++) {
        RCP<const SolutionState<double> > solutionState = (*solutionHistory)[i];
        const double time_i = solutionState->getTime();
        RCP<const DPV> x_prod_plot =
          Teuchos::rcp_dynamic_cast<const DPV>(solutionState->getX());
        RCP<const Thyra::VectorBase<double> > x_plot =
          x_prod_plot->getVectorBlock(0);
        RCP<const DMVPV > adjoint_prod_plot =
          Teuchos::rcp_dynamic_cast<const DMVPV>(x_prod_plot->getVectorBlock(1));
        RCP<const Thyra::MultiVectorBase<double> > adjoint_plot =
          adjoint_prod_plot->getMultiVector();
        RCP<const Thyra::VectorBase<double> > x_exact_plot =
          model->getExactSolution(time_i).get_x();
        ftmp << std::fixed << std::setprecision(7)
             << time_i
             << std::setw(11) << get_ele(*(x_plot), 0)
             << std::setw(11) << get_ele(*(x_plot), 1)
             << std::setw(11) << get_ele(*(adjoint_plot->col(0)), 0)
             << std::setw(11) << get_ele(*(adjoint_plot->col(0)), 1)
             << std::setw(11) << get_ele(*(adjoint_plot->col(1)), 0)
             << std::setw(11) << get_ele(*(adjoint_plot->col(1)), 1)
             << std::setw(11) << get_ele(*(x_exact_plot), 0)
             << std::setw(11) << get_ele(*(x_exact_plot), 1)
             << std::endl;
      }
      ftmp.close();
    }

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    RCP<Thyra::MultiVectorBase<double> > DxDpdiff = DxDp->clone_mv();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
    Thyra::V_VmV(DxDpdiff.ptr(), *DxDp_exact, *DxDp);
    StepSize.push_back(dt);
    double L2norm = Thyra::norm_2(*xdiff);
    L2norm *= L2norm;
    Teuchos::Array<double> L2norm_DxDp(num_param);
    Thyra::norms_2(*DxDpdiff, L2norm_DxDp());
    for (int i=0; i<num_param; ++i)
      L2norm += L2norm_DxDp[i]*L2norm_DxDp[i];
    L2norm = std::sqrt(L2norm);
    ErrorNorm.push_back(L2norm);

    //*my_out << " n = " << n << " dt = " << dt << " error = " << L2norm
    //        << std::endl;
  }

  // Check the order and intercept
  double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
  *my_out << "  Stepper = BDF2" << std::endl;
  *my_out << "  =========================" << std::endl;
  *my_out << "  Expected order: " << order << std::endl;
  *my_out << "  Observed order: " << slope << std::endl;
  *my_out << "  =========================" << std::endl;
  TEST_FLOATING_EQUALITY( slope,          1.95006, 0.015 );
  TEST_FLOATING_EQUALITY( ErrorNorm[0], 0.0137394, 1.0e-4 );

  if (comm->getRank() == 0) {
    std::ofstream ftmp("Tempus_BDF2_SinCos_AdjSens-Error.dat");
    double error0 = 0.8*ErrorNorm[0];
    for (int n=0; n<nTimeStepSizes; n++) {
      ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
           << error0*(StepSize[n]/StepSize[0]) << std::endl;
    }
    ftmp.close();
  }

}

} // namespace Tempus_Test
