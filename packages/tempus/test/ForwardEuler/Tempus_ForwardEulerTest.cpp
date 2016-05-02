#include <Teuchos_UnitTestHarness.hpp>
#include "Tempus_IntegratorBasic.hpp"

using Teuchos::RCP;
using Teuchos::ParameterList;
using Teuchos::getParametersFromXmlFile;

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(ForwardEuler, SinCos)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("ForwardEuler_SinCos.xml");

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = Teuchos::sublist(pList, "SinCosModel");
  RCP<SinCosModel> model    = sinCosModel(scm_pl);

  // Setup the Integrator
  RCP<ParameterList> i_pl = Teuchos::sublist(pList, "Tempus");
  RCP<IntegratorBasic<double> > integrator = IntegratorBasic(i_pl, model)
}

}
