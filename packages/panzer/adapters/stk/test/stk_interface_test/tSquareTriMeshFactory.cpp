#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tSquareTriMeshFactory, defaults)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;
   
   SquareTriMeshFactory factory; 
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
 
   if(mesh->isWritable());
      mesh->writeToExodus("square-tri.exo");

   // minimal requirements
   TEST_ASSERT(mesh!=Teuchos::null);
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),4);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),2*25);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),25+60);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),36);

   int numprocs = stk::parallel_machine_size(MPI_COMM_WORLD);
   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
}

}
