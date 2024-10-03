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
#include "Panzer_STK_SquareTriMeshFactory.hpp"

#include "Shards_BasicTopologies.hpp"

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"

#include "Ioss_DatabaseIO.h"
#include "Ioss_IOFactory.h"
#include "Ioss_Region.h"
#include "Ioss_EdgeBlock.h"

namespace panzer_stk {

void edge_block_test_helper(Teuchos::FancyOStream &out,
                            bool &success,
                            Teuchos::RCP<Teuchos::ParameterList> pl,
                            std::string exodus_filename,
                            uint32_t expected_edge_block_count)
{
   SquareTriMeshFactory factory;
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

   if (expected_edge_block_count == 1) {
      std::vector<stk::mesh::Entity> edges;
      mesh->getMyEdges(edges);
      TEST_ASSERT(all_edge_blocks[0]->entity_count() == (int64_t)edges.size());
   }
   }
}

TEUCHOS_UNIT_TEST(tSquareTriMeshFactory, defaults)
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcpFromRef;

   SquareTriMeshFactory factory;
   RCP<STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

   if(mesh->isWritable())
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

   int mpi_numprocs = -1;
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);
   int mpi_rank = -1;
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
   TEST_EQUALITY(numprocs,mpi_numprocs);
   TEST_EQUALITY(rank,mpi_rank);

   // check for nodeset
   std::vector<std::string> nodesets;
   mesh->getNodesetNames(nodesets);

   TEST_EQUALITY(nodesets.size(),1);
   TEST_EQUALITY(nodesets[0],"origin");
}

TEUCHOS_UNIT_TEST(tSquareTriMeshFactory, default_edge_face_blocks)
{
   using Teuchos::RCP;

   int xe = 2, ye = 2;
   int bx = 1, by = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);

   std::size_t expected_edge_block_count = 0;

   edge_block_test_helper(out, success, pl,
                          "SquareTri_default_edge_blocks.exo",
                          expected_edge_block_count);
}

TEUCHOS_UNIT_TEST(tSquareTriMeshFactory, create_edge_blocks_pl)
{
   using Teuchos::RCP;

   int xe = 2, ye = 2;
   int bx = 1, by = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Create Edge Blocks",true);

   std::size_t expected_edge_block_count = 1;

   edge_block_test_helper(out, success, pl,
                          "SquareTri_create_edge_blocks_pl.exo",
                          expected_edge_block_count);
}

TEUCHOS_UNIT_TEST(tSquareTriMeshFactory, multiblock_create_edge_blocks_pl)
{
   using Teuchos::RCP;

   int xe = 2, ye = 2;
   int bx = 2, by = 1;

   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
   pl->set("X Blocks",bx);
   pl->set("Y Blocks",by);
   pl->set("X Elements",xe);
   pl->set("Y Elements",ye);
   pl->set("Create Edge Blocks",true);

   std::size_t expected_edge_block_count = 1;

   edge_block_test_helper(out, success, pl,
                          "SquareTri_create_edge_blocks_pl.exo",
                          expected_edge_block_count);
}
}
