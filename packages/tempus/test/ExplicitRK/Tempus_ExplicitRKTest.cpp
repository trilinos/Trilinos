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

TEUCHOS_UNIT_TEST(ExplicitRK, SinCos)
{
  std::vector<std::string> RKMethods;
  RKMethods.push_back("RK Forward Euler");
  RKMethods.push_back("RK Explicit 4 Stage");
  RKMethods.push_back("RK Explicit 3/8 Rule");
  RKMethods.push_back("RK Explicit 4 Stage 3rd order by Runge");
  RKMethods.push_back("RK Explicit 5 Stage 3rd order by Kinnmark and Gray");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order TVD");
  RKMethods.push_back("RK Explicit 3 Stage 3rd order by Heun");
  RKMethods.push_back("RK Explicit 2 Stage 2nd order by Runge");
  RKMethods.push_back("RK Explicit Trapezoidal");
  RKMethods.push_back("General ERK");
  std::vector<double> RKMethodErrors;
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
  RKMethodErrors.push_back(8.33251e-07);

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
        getParametersFromXmlFile("Tempus_ExplicitRK_SinCos.xml");

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
      pl->sublist("Demo Stepper").set("Stepper Type", RKMethods[m]);

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

    std::ofstream ftmp("Tempus_"+RKMethod_+"_SinCos-Error.dat");
    double error0 = 0.8*ErrorNorm[0];
    for (int n=0; n<nTimeStepSizes; n++) {
      ftmp << StepSize[n]  << "   " << ErrorNorm[n] << "   "
           << error0*(pow(StepSize[n]/StepSize[0],order)) << std::endl;
    }
    ftmp.close();

  }
}

} // namespace Tempus_Test
