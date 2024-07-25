// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_LineMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

namespace panzer_stk {

TEUCHOS_UNIT_TEST(tLineMeshFactory, defaults)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   LineMeshFactory factory;
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   TEST_EQUALITY(mesh->getPeriodicBCVector().size(),0);

   if(mesh->isWritable())
      mesh->writeToExodus("Line.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),1);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),6);
}

TEUCHOS_UNIT_TEST(tLineMeshFactory, element_counts)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("X Elements",2);

   LineMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   if(mesh->isWritable())
      mesh->writeToExodus("Line_oddelmt.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),1);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),3);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),3);
}

TEUCHOS_UNIT_TEST(tLineMeshFactory, allblock)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int xe = 4;
   int bx = 4;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("X Elements",xe);

   xe *= bx;

   LineMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   if(mesh->isWritable())
      mesh->writeToExodus("Line_allblock.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),1);
   TEST_EQUALITY(mesh->getNumElementBlocks(),(std::size_t) bx);
   TEST_EQUALITY(mesh->getNumSidesets(),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),(std::size_t) xe);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),(std::size_t) xe+1);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(std::size_t) xe+1);
}

TEUCHOS_UNIT_TEST(tLineMeshFactory, two_block)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("X Elements",5);

   LineMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   if(mesh->isWritable())
      mesh->writeToExodus("Line_2block.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),1);
   TEST_EQUALITY(mesh->getNumElementBlocks(),2);
   TEST_EQUALITY(mesh->getNumSidesets(),2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),10);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),11);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),11);
}

}
