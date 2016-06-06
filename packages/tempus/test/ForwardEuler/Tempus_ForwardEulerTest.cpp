#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;
using Teuchos::Array;

using Tempus::IntegratorBasic;

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(ForwardEuler, SinCos)
{
  Array<double> logStepSize;
  Array<double> logErrorNorm;
  const int nTimeStepSizes = 7;
  double dt = 0.2;
  for (int n=0; n<nTimeStepSizes; n++) {

    // Read params from .xml file
    RCP<ParameterList> pList =
      getParametersFromXmlFile("ForwardEuler_SinCos.xml");

    // Setup the SinCosModel
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    //RCP<SinCosModel> model    = sineCosineModel(scm_pl);
    RCP<SinCosModel> model    = Teuchos::rcp(new SinCosModel(scm_pl));

    dt /= 2;

    // Setup the Integrator and reset initial time step
    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    pl->sublist("Integrator Basic").set("Initial Time Step", dt);
    RCP<Tempus::IntegratorBasic<double> > integrator =
      integratorBasic<double>(pl, model);

    // Integrate to timeMax
    bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus)

    // Test if at 'Final Time'
    double time = integrator->getTime();
    double timeFinal =pl->sublist("Integrator Basic").get<double>("Final Time");
    TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

    // Time-integrated solution and the exact solution
    Teuchos::RCP<Thyra::VectorBase<double> > x = integrator->getX();
    RCP<const Thyra::VectorBase<double> > x_exact =
      model->getExactSolution(time).get_x();

    // Calculate the error
    RCP<Thyra::VectorBase<double> > xdiff = x->clone_v();
    Thyra::V_StVpStV(xdiff.ptr(), 1.0, *x_exact, -1.0, *(x));
    logStepSize.push_back(log(dt));
    const double L2norm = Thyra::norm_2(*xdiff);
    logErrorNorm.push_back(log(L2norm));
  }

  // Check the order and intercept
  double slope = 0.0;
  double yIntercept = 0.0;
  computeLinearRegression<double>(logStepSize,logErrorNorm,slope,yIntercept);
  TEST_FLOATING_EQUALITY( slope, 1.0, 0.01 );
  TEST_FLOATING_EQUALITY( exp(logErrorNorm[0]), 0.051123, 1.0e-4 );
}

} // namespace Tempus_Test
