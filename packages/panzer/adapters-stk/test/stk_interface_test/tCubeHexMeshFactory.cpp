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
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

#include "Shards_BasicTopologies.hpp"

#include "Ioss_DatabaseIO.h"
#include "Ioss_IOFactory.h"
#include "Ioss_Region.h"
#include "Ioss_EdgeBlock.h"
#include "Ioss_FaceBlock.h"

namespace panzer_stk {

void edge_face_block_test_helper(Teuchos::FancyOStream &out,
                                 bool &success,
                                 Teuchos::RCP<Teuchos::ParameterList> pl,
                                 std::string exodus_filename,
                                 uint32_t expected_edge_block_count,
                                 uint32_t expected_face_block_count)
{
   CubeHexMeshFactory factory;
   factory.setParameterList(pl);
   Teuchos::RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   if(mesh->isWritable())
      mesh->writeToExodus(exodus_filename.c_str());

   {
   Ioss::DatabaseIO *db_io = Ioss::IOFactory::create("exodus",
                                                     exodus_filename.c_str(),
                                                     Ioss::READ_MODEL);
   TEST_ASSERT(db_io);

   Ioss::Region region(db_io);
   TEST_ASSERT(db_io->ok() == true);

   auto all_edge_blocks = region.get_edge_blocks();
   TEST_ASSERT(all_edge_blocks.size() == expected_edge_block_count);
   auto all_face_blocks = region.get_face_blocks();
   TEST_ASSERT(all_face_blocks.size() == expected_face_block_count);

   if (expected_edge_block_count == 1) {
      std::vector<stk::mesh::Entity> edges;
      mesh->getMyEdges(edges);
      TEST_ASSERT(all_edge_blocks[0]->entity_count() == (int64_t)edges.size());
   }
   if (expected_face_block_count == 1) {
      std::vector<stk::mesh::Entity> faces;
      mesh->getMyFaces(faces);
      TEST_ASSERT(all_face_blocks[0]->entity_count() == (int64_t)faces.size());
   }
   }
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, defaults)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   CubeHexMeshFactory factory;
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   TEST_EQUALITY(mesh->getPeriodicBCVector().size(),0);

   if(mesh->isWritable())
      mesh->writeToExodus("CubeHex.exo");

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

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, periodic_input)
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
   Teuchos::ParameterList & pbcs = pl->sublist("Periodic BCs");
   pbcs.set<int>("Count",1);
   pbcs.set("Periodic Condition 1","yz-coord left;right");

   CubeHexMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   TEST_EQUALITY(mesh->getPeriodicBCVector().size(),1);
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

   if(mesh->isWritable())
      mesh->writeToExodus("CubeHex_oddelmt.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),2*(4+1)*(5+1)+4*(2+1)*(5+1)+5*(2+1)*(4+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getFaceRank()),2*4*(5+1)+2*5*(4+1)+4*5*(2+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, disable_subcells)
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
   pl->set("Build Subcells",false);

   CubeHexMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   if(mesh->isWritable())
      mesh->writeToExodus("CubeHex_disable_subcells.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),4*2*5);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),2*4*2 + 2*5*2 + 4*5*2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),0);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getFaceRank()),2*4*2 + 2*5*2 + 4*5*2);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(4+1)*(2+1)*(5+1));
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, default_edge_face_blocks)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int xe = 2, ye = 2, ze = 2;
   int bx = 1, by = 1, bz = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("Z Blocks",bz);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Z Elements",ze);

   std::size_t expected_edge_block_count = 0;
   std::size_t expected_face_block_count = 0;

   edge_face_block_test_helper(out, success, pl,
                               "CubeHex_default_edge_face_blocks.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, create_edge_blocks_pl)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int xe = 2, ye = 2, ze = 2;
   int bx = 1, by = 1, bz = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("Z Blocks",bz);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Z Elements",ze);
   pl->set("Create Edge Blocks",true);

   std::size_t expected_edge_block_count = 1;
   std::size_t expected_face_block_count = 0;

   edge_face_block_test_helper(out, success, pl,
                               "CubeHex_create_edge_blocks_pl.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, create_face_blocks_pl)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int xe = 2, ye = 2, ze = 2;
   int bx = 1, by = 1, bz = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("Z Blocks",bz);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Z Elements",ze);
   pl->set("Create Face Blocks",true);

   std::size_t expected_edge_block_count = 0;
   std::size_t expected_face_block_count = 1;

   edge_face_block_test_helper(out, success, pl,
                               "CubeHex_create_face_blocks_pl.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}


TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, create_edge_face_blocks_pl)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int xe = 2, ye = 2, ze = 2;
   int bx = 1, by = 1, bz = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("Z Blocks",bz);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Z Elements",ze);
   pl->set("Create Edge Blocks",true);
   pl->set("Create Face Blocks",true);

   std::size_t expected_edge_block_count = 1;
   std::size_t expected_face_block_count = 1;

   edge_face_block_test_helper(out, success, pl,
                               "CubeHex_create_edge_face_blocks_pl.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, multiblock_create_edge_face_blocks_pl)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int xe = 2, ye = 2, ze = 2;
   int bx = 2, by = 1, bz = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("Z Blocks",bz);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Z Elements",ze);
   pl->set("Create Edge Blocks",true);
   pl->set("Create Face Blocks",true);

   std::size_t expected_edge_block_count = 1;
   std::size_t expected_face_block_count = 1;

   edge_face_block_test_helper(out, success, pl,
                               "CubeHex_multiblock_create_edge_face_blocks_pl.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, allblock)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   int rank = stk::parallel_machine_rank(MPI_COMM_WORLD);

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

   if(mesh->isWritable())
      mesh->writeToExodus("CubeHex_allblock.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),(std::size_t) bx*by*bz);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),(std::size_t) xe*ye*ze);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getSideRank()),(std::size_t) xe*ye*(ze+1)+xe*(ye+1)*ze+(xe+1)*ye*ze);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getEdgeRank()),(std::size_t) xe*(ye+1)*(ze+1)+(xe+1)*(ye+1)*ze+(xe+1)*ye*(ze+1));
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),(std::size_t) (xe+1)*(ye+1)*(ze+1));

   std::vector<std::string> sidesets, nodesets;
   mesh->getNodesetNames(nodesets);
   TEST_EQUALITY(nodesets.size(),1);

   std::vector<stk::mesh::Entity> nodes;
   mesh->getMyNodes("origin","eblock-0_0_0",nodes);
   if(rank==0) {
      std::vector<std::size_t> localNodeIds;
      std::vector<stk::mesh::Entity> elements;

      TEST_EQUALITY(nodes.size(),1);
      workset_utils::getNodeElements(*mesh,"eblock-0_0_0",nodes,localNodeIds,elements);

      TEST_EQUALITY(localNodeIds.size(),1);
      TEST_EQUALITY(elements.size(),1);
      TEST_EQUALITY(mesh->elementGlobalId(elements[0]),1);
      TEST_EQUALITY(localNodeIds[0],0);
   }
   else {
      TEST_EQUALITY(nodes.size(),0);
   }
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

   if(mesh->isWritable())
      mesh->writeToExodus("CubeHex_2block.exo");

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

TEUCHOS_UNIT_TEST(tCubeHexMeshFactory, sub_two_block)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;
   int size; MPI_Comm_size(MPI_COMM_WORLD, &size);

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",2);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Procs",2);
   pl->set("Y Procs",2);
   pl->set("Z Procs",2);
   pl->set("X Elements",5);
   pl->set("Y Elements",10);
   pl->set("Z Elements",5);

   CubeHexMeshFactory factory;
   factory.setParameterList(pl);
   RCP<STK_Interface> mesh;
   if(size!=8) {
      TEST_THROW(factory.buildMesh(MPI_COMM_WORLD),std::logic_error);
      return;
   }
   else {
      mesh = factory.buildMesh(MPI_COMM_WORLD);
   }
   TEST_ASSERT(mesh!=Teuchos::null);

   if(mesh->isWritable())
      mesh->writeToExodus("CubeHex_sub_2block.exo");

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
