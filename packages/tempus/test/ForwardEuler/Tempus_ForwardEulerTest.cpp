#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tempus_IntegratorBasic.hpp"
//#include "Tempus_IntegratorBasic_impl.hpp"

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
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel");
  //RCP<SinCosModel> model    = sineCosineModel(scm_pl);
  RCP<SinCosModel> model    = Teuchos::rcp(new SinCosModel(scm_pl));

  // Setup the Integrator
  RCP<ParameterList> pl = sublist(pList, "Tempus");
  const Teuchos::RCP<Thyra::NonlinearSolverBase<double> >& solver=Teuchos::null;
  //RCP<Tempus::IntegratorBasic<double> > integrator = intBasic<double>(pl, model, solver);
  RCP<Tempus::IntegratorBasic<double> > integrator = Teuchos::rcp(new Tempus::IntegratorBasic<double>(pl, model,solver));

}

}
