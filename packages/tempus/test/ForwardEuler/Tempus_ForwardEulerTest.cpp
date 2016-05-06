#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include "../TestModels/SinCosModel.hpp"

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(ForwardEuler, SinCos)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("ForwardEuler_SinCos.xml");

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  //RCP<SinCosModel> model    = sineCosineModel(scm_pl);
  RCP<SinCosModel> model    = Teuchos::rcp(new SinCosModel(scm_pl));

  // Setup the Integrator
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  RCP<Tempus::IntegratorBasic<double> > integrator =
    integratorBasic<double>(pl, model);
}

} // namespace Tempus_Test
