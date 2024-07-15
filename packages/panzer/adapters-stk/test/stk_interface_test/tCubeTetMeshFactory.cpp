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
#include "Panzer_STK_CubeTetMeshFactory.hpp"

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
   CubeTetMeshFactory factory;
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

TEUCHOS_UNIT_TEST(tCubeTetMeshFactory, defaults)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

/*
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",1);
   pl->set("Y Blocks",1);
   pl->set("Z Blocks",1);
   pl->set("X Elements",2);
   pl->set("Y Elements",1);
   pl->set("Z Elements",1);
*/

   CubeTetMeshFactory factory;
   // factory.setParameterList(pl);
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
   TEST_ASSERT(mesh!=Teuchos::null);

   TEST_EQUALITY(mesh->getPeriodicBCVector().size(),0);

   if(mesh->isWritable())
      mesh->writeToExodus("CubeTet.exo");

   // minimal requirements
   TEST_ASSERT(not mesh->isModifiable());

   TEST_EQUALITY(mesh->getDimension(),3);
   TEST_EQUALITY(mesh->getNumElementBlocks(),1);
   TEST_EQUALITY(mesh->getNumSidesets(),6);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getElementRank()),12*125);
   TEST_EQUALITY(mesh->getEntityCounts(mesh->getNodeRank()),6*6*6+5*5*5);

   // check for nodeset
   std::vector<std::string> nodesets;
   mesh->getNodesetNames(nodesets);

   TEST_EQUALITY(nodesets.size(),1);
   TEST_EQUALITY(nodesets[0],"origin");
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
                               "CubeTet_default_edge_face_blocks.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}

TEUCHOS_UNIT_TEST(tCubeTetMeshFactory, create_edge_blocks_pl)
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
                               "CubeTet_create_edge_blocks_pl.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}

TEUCHOS_UNIT_TEST(tCubeTetMeshFactory, create_face_blocks_pl)
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
                               "CubeTet_create_face_blocks_pl.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}


TEUCHOS_UNIT_TEST(tCubeTetMeshFactory, create_edge_face_blocks_pl)
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
                               "CubeTet_create_edge_face_blocks_pl.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}


TEUCHOS_UNIT_TEST(tCubeTetMeshFactory, multiblock_create_edge_face_blocks_pl)
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
                               "CubeTet_multiblock_create_edge_face_blocks_pl.exo",
                               expected_edge_block_count,
                               expected_face_block_count);
}

}
