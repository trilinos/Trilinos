#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Tuple.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_STK_PeriodicBC_Matcher.hpp"

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

namespace panzer {

  TEUCHOS_UNIT_TEST(periodic_mesh, add_get_vector)
  {
    using Teuchos::RCP;
    using Teuchos::Tuple;

    panzer_stk::SquareQuadMeshFactory mesh_factory;

    // setup mesh
    /////////////////////////////////////////////
    RCP<panzer_stk::STK_Interface> mesh;
    {
       RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",1);
       pl->set("X Elements",6);
       pl->set("Y Elements",4);
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    CoordMatcher x_matcher(0);
    CoordMatcher y_matcher(1);
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("top","bottom",x_matcher));
    mesh->addPeriodicBC(panzer_stk::buildPeriodicBC_Matcher("left","right",y_matcher));

    std::vector<RCP<const panzer_stk::PeriodicBC_MatcherBase> > & mod_vec = mesh->getPeriodicBCVector();
    TEST_EQUALITY(mod_vec.size(),2);

    const std::vector<RCP<const panzer_stk::PeriodicBC_MatcherBase> > & const_vec = mesh.getConst()->getPeriodicBCVector();
    TEST_EQUALITY(const_vec.size(),2);
  }

}
