#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, defaults)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;
   
   CubeHexMeshFactory factory; 
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   if(mesh->isWritable());
      mesh->writeToExodus("Cube.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),125);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),3*25*6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),3*30*6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),6*6*6);
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, element_counts)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",4);
   pl->set("Z Elements",5);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   if(mesh->isWritable());
      mesh->writeToExodus("Cube_oddelmt.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, allblock)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int xe = 4, ye = 5, ze = 2;
   int bx = 4, by = 2, bz = 3;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("Z Blocks",bz);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Z Elements",ze);

   xe *= bx; ye *= by; ze *= bz;
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   if(mesh->isWritable());
      mesh->writeToExodus("Cube_allblock.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),(std::size_t) bx*by*bz);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),(std::size_t) xe*ye*ze);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),(std::size_t) xe*ye*(ze+1)+xe*(ye+1)*ze+(xe+1)*ye*ze);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),(std::size_t) xe*(ye+1)*(ze+1)+(xe+1)*(ye+1)*ze+(xe+1)*ye*(ze+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(std::size_t) (xe+1)*(ye+1)*(ze+1));
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, two_block)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",5);
   pl->set("Y Elements",10);
   pl->set("Z Elements",5);
   
   CubeHexMeshFactory factory; 
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);
 
   if(mesh->isWritable());
      mesh->writeToExodus("Cube_2block.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),2);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),2*5*10*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),10*10*6+10*5*11+10*5*11);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),10*11*6+10*6*11+11*5*11);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),11*11*6);
}

}
