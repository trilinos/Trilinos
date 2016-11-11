#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <vector>

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(DynamicERK, SinCos)
{

    const int nTimeStepSizes = 7;
    double dt = 0.2;
    double order = 0.0;
    std::string RKMethod_ = "dynamic";
    std::vector<double> StepSize;
    std::vector<double> ErrorNorm;

    for (int n = 0; n < nTimeStepSizes; n++){
      // Read params from .xml file
      RCP<ParameterList> pList =
        getParametersFromXmlFile("Tempus_DynamicRK_SinCos.xml");

      //std::ofstream ftmp("PL.txt");
      //pList->print(ftmp);
      //ftmp.close();

      // Setup the SinCosModel
      RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
      //RCP<SinCosModel<double> > model = sineCosineModel(scm_pl);
      RCP<SinCosModel<double> > model =
        Teuchos::rcp(new SinCosModel<double>(scm_pl));

      // Set the Stepper
      RCP<ParameterList> pl = sublist(pList, "Tempus", true);
      dt /= 2;

      // Setup the Integrator and reset initial time step
      pl->sublist("Demo Integrator").set("Initial Time Step", dt);

      RCP<Tempus::IntegratorBasic<double> > integrator =
        Tempus::integratorBasic<double>(pl, model);
      order = integrator->getStepper()->getOrder();

      // Integrate to timeMax
      bool integratorStatus = integrator->advanceTime();
      TEST_ASSERT(integratorStatus)

      // Test if at 'Final Time'
      double time = integrator->getTime();
      double timeFinal = pl->sublist("Demo Integrator").get<double>("Final Time");
      TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

/*
      // Time-integrated solution and the exact solution
      RCP<Thyra::VectorBase<double> > x = integrator->getX();
      RCP<const Thyra::VectorBase<double> > x_exact =
        model->getExactSolution(time).get_x();

      // Plot sample solution and exact solution
      if (n == 0) {
        std::ofstream ftmp("Tempus_"+RKMethod_+"_SinCos.dat");
        RCP<SolutionHistory<double> > solutionHistory =
          integrator->getSolutionHistory();
        int nStates = solutionHistory->getNumStates();
        RCP<const Thyra::VectorBase<double> > x_exact_plot;
        for (int i=0; i<nStates; i++) {
          RCP<SolutionState<double> > solutionState = (*solutionHistory)[i];
          double time = solutionState->getTime();
          RCP<Thyra::VectorBase<double> > x_plot = solutionState->getX();
          x_exact_plot = model->getExactSolution(time).get_x();
          ftmp << time << "   "
               << get_ele(*(x_plot), 0) << "   "
               << get_ele(*(x_plot), 1) << "   "
               << get_ele(*(x_exact_plot), 0) << "   "
               << get_ele(*(x_exact_plot), 1) << std::endl;
        }
        ftmp.close();
      }

      // Calculate the error
      RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
      Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
      StepSize.push_back(dt);
      const double L2norm = Thyra::norm_2(*xdiff);
      ErrorNorm.push_back(L2norm);

    // Check the order and intercept
    double slope = computeLinearRegressionLogLog<double>(StepSize, ErrorNorm);
    std::cout << "  Stepper = " << RKMethod_ << std::endl;
    std::cout << "  =========================" << std::endl;
    std::cout << "  Expected order: " << order << std::endl;
    std::cout << "  Observed order: " << slope << std::endl;
    std::cout << "  =========================" << std::endl;
    TEST_FLOATING_EQUALITY( slope, order, 0.01 );
    TEST_FLOATING_EQUALITY( ErrorNorm[0], RKMethodErrors[m], 1.0e-4 );

    std::ofstream ftmp("Tempus_"+RKMethod_+"_SinCos-Error.dat");
    double error0 = 0.8*ErrorNorm[0];
    for (int n=0; n<nTimeStepSizes; n++) {
      ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
           << error0*(pow(StepSize[n]/StepSize[0],order)) << std::endl;
    }
    ftmp.close();
    */

    }
  }

} // namespace Tempus_Test
