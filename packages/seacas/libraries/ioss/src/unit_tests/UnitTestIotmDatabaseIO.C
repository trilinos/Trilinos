// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "Ionit_Initializer.h"
#include "Ioss_DBUsage.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_Hex8.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_Shell4.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "text_mesh/Iotm_DatabaseIO.h"
#include "text_mesh/Iotm_TextMeshTopologyMapping.h"
#include <fmt/core.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "Ioss_CodeTypes.h"
#include "Ioss_Field.h"
#include "Ioss_ParallelUtils.h"

namespace {

  Iotm::DatabaseIO *create_input_db_io(const std::string &meshDesc)
  {
    Ioss::Init::Initializer init_db;

    Ioss::DatabaseUsage   db_usage = Ioss::READ_MODEL;
    Ioss::PropertyManager properties;

    properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));

    auto *db_io = new Iotm::DatabaseIO(nullptr, meshDesc, db_usage,
                                       Ioss::ParallelUtils::comm_world(), properties);
    return db_io;
  }

  Iotm::DatabaseIO *create_output_db_io(const std::string &filename)
  {
    Ioss::Init::Initializer init_db;
    Ioss::DatabaseUsage     db_usage = Ioss::WRITE_RESULTS;
    Ioss::PropertyManager   properties;

    properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));

    auto *db_io = new Iotm::DatabaseIO(nullptr, filename, db_usage,
                                       Ioss::ParallelUtils::comm_world(), properties);
    return db_io;
  }

  int get_parallel_size()
  {
    return Ioss::ParallelUtils(Ioss::ParallelUtils::comm_world()).parallel_size();
  }

  int get_parallel_rank()
  {
    return Ioss::ParallelUtils(Ioss::ParallelUtils::comm_world()).parallel_rank();
  }

  bool include_entity(const Ioss::GroupingEntity *entity)
  {
    assert(entity);

    // Check whether entity has "omitted" property...
    bool omitted =
        (entity->property_exists("omitted")) && (entity->get_property("omitted").get_int() == 1);

    return !omitted;
  }

  TEST(TextMesh, twoHexesSerial)
  {
    if (get_parallel_size() != 1) {
      GTEST_SKIP();
    }

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";

    Iotm::DatabaseIO *db_io = create_input_db_io(meshDesc);
    db_io->set_surface_split_type(Ioss::SPLIT_BY_ELEMENT_BLOCK);

    Ioss::Region region(db_io);

    EXPECT_TRUE(nullptr != db_io);
    EXPECT_TRUE(db_io->ok());
    EXPECT_EQ("TextMesh", db_io->get_format());

    const std::vector<Ioss::ElementBlock *> &element_blocks = region.get_element_blocks();
    EXPECT_EQ(2u, element_blocks.size());

    EXPECT_EQ(1u, element_blocks[0]->entity_count());
    EXPECT_EQ(1u, element_blocks[1]->entity_count());

    EXPECT_EQ(Ioss::Hex8::name, element_blocks[0]->topology()->name());
    EXPECT_EQ(Ioss::Hex8::name, element_blocks[1]->topology()->name());
  }

  TEST(TextMesh, twoHexesParallel)
  {
    if (get_parallel_size() != 2) {
      GTEST_SKIP();
    }

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";

    Iotm::DatabaseIO *db_io = create_input_db_io(meshDesc);
    db_io->set_surface_split_type(Ioss::SPLIT_BY_ELEMENT_BLOCK);

    Ioss::Region region(db_io);

    EXPECT_TRUE(nullptr != db_io);
    EXPECT_TRUE(db_io->ok());
    EXPECT_EQ("TextMesh", db_io->get_format());

    const std::vector<Ioss::ElementBlock *> &element_blocks = region.get_element_blocks();
    EXPECT_EQ(2u, element_blocks.size());

    if (get_parallel_rank() == 0) {
      EXPECT_EQ(1u, element_blocks[0]->entity_count());
      EXPECT_EQ(0u, element_blocks[1]->entity_count());
    }
    else {
      EXPECT_EQ(0u, element_blocks[0]->entity_count());
      EXPECT_EQ(1u, element_blocks[1]->entity_count());
    }

    EXPECT_EQ(Ioss::Hex8::name, element_blocks[0]->topology()->name());
    EXPECT_EQ(Ioss::Hex8::name, element_blocks[1]->topology()->name());
  }

  TEST(TextMesh, twoHexesParallel_skipBlock1)
  {
    if (get_parallel_size() != 2) {
      GTEST_SKIP();
    }

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";

    const std::vector<std::string> omittedBlocks{"block_1"};

    Iotm::DatabaseIO *db_io = create_input_db_io(meshDesc);
    db_io->set_surface_split_type(Ioss::SPLIT_BY_ELEMENT_BLOCK);
    db_io->set_block_omissions(omittedBlocks);

    Ioss::Region region(db_io);

    EXPECT_TRUE(nullptr != db_io);
    EXPECT_TRUE(db_io->ok());
    EXPECT_EQ("TextMesh", db_io->get_format());

    const std::vector<Ioss::ElementBlock *> &element_blocks = region.get_element_blocks();
    EXPECT_EQ(2u, element_blocks.size());

    EXPECT_FALSE(include_entity(element_blocks[0]));
    EXPECT_TRUE(include_entity(element_blocks[1]));
  }

  TEST(TextMesh, surfaceToBlockMapping_noSplit)
  {
    if (get_parallel_size() != 2) {
      GTEST_SKIP();
    }

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "1,2,HEX_8,2,9,10,3,6,11,12,7,block_2\n"
                           "|sideset:name=left_surf;data=1,4";

    Iotm::DatabaseIO *db_io = create_input_db_io(meshDesc);

    Ioss::Region region(db_io);

    EXPECT_TRUE(nullptr != db_io);
    EXPECT_TRUE(db_io->ok());
    EXPECT_EQ("TextMesh", db_io->get_format());

    const std::vector<Ioss::ElementBlock *> &element_blocks = region.get_element_blocks();
    EXPECT_EQ(2u, element_blocks.size());

    const std::vector<Ioss::SideSet *> &sidesets = region.get_sidesets();
    EXPECT_EQ(1u, sidesets.size());

    {
      std::string      sideblockName("LEFT_SURF");
      Ioss::SideBlock *sideblock = sidesets[0]->get_side_block(sideblockName);
      EXPECT_TRUE(nullptr != sideblock);

      std::vector<std::string> touchingBlocks;
      db_io->compute_block_membership(sideblock, touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_1"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
  }

  TEST(TextMesh, surfaceToBlockMapping_splitByBlock)
  {
    if (get_parallel_size() != 1) {
      GTEST_SKIP();
    }

    std::string meshDesc = "0,1,TRI_3_2D,3,1,4,block_1\n"
                           "0,2,TRI_3_2D,1,2,4,block_1\n"
                           "0,3,TRI_3_2D,2,5,4,block_1\n"
                           "0,4,TRI_3_2D,5,7,4,block_2\n"
                           "0,5,TRI_3_2D,7,6,4,block_2\n"
                           "0,6,TRI_3_2D,6,3,4,block_2\n"
                           "|coordinates: 0,0,0.1,0,0,0.1,0.05,0.1,0.1,0.1,0,0.2,0.1,0.2"
                           "|dimension:2"
                           "|sideset:name=skinned_surf; skin=all; split=block"
                           "|sideset:name=shared_surf; data=1,1,6,1; split=block"
                           "|sideset:name=owned_surf; data=5,1; split=block";

    Iotm::DatabaseIO *db_io = create_input_db_io(meshDesc);

    Ioss::Region region(db_io);

    EXPECT_TRUE(nullptr != db_io);
    EXPECT_TRUE(db_io->ok());
    EXPECT_EQ("TextMesh", db_io->get_format());

    const std::vector<Ioss::ElementBlock *> &element_blocks = region.get_element_blocks();
    EXPECT_EQ(2u, element_blocks.size());

    const std::vector<Ioss::SideSet *> &sidesets = region.get_sidesets();
    EXPECT_EQ(3u, sidesets.size());

    Iotm::IossTopologyMapping topologyMapping;
    topologyMapping.initialize_topology_map();

    auto get_topology_name = [&topologyMapping](const std::string &textMeshTopologyName) {
      return topologyMapping.topology(textMeshTopologyName).name();
    };

    Ioss::SideSet *skinned_surf = region.get_sideset("SKINNED_SURF");
    EXPECT_TRUE(nullptr != skinned_surf);
    {
      std::vector<std::string> touchingBlocks;
      skinned_surf->block_membership(touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_1", "BLOCK_2"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
    {
      std::string sideblockName("SKINNED_SURF_BLOCK_1_" + get_topology_name("LINE_2"));
      sideblockName = Ioss::Utils::uppercase(sideblockName);

      Ioss::SideBlock *sideblock = skinned_surf->get_side_block(sideblockName);
      EXPECT_TRUE(nullptr != sideblock);

      std::vector<std::string> touchingBlocks;
      db_io->compute_block_membership(sideblock, touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_1"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
    {
      std::string sideblockName("SKINNED_SURF_BLOCK_2_" + get_topology_name("LINE_2"));
      sideblockName = Ioss::Utils::uppercase(sideblockName);

      Ioss::SideBlock *sideblock = skinned_surf->get_side_block(sideblockName);
      EXPECT_TRUE(nullptr != sideblock);

      std::vector<std::string> touchingBlocks;
      db_io->compute_block_membership(sideblock, touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_2"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }

    Ioss::SideSet *shared_surf = region.get_sideset("SHARED_SURF");
    EXPECT_TRUE(nullptr != shared_surf);
    {
      std::vector<std::string> touchingBlocks;
      shared_surf->block_membership(touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_1", "BLOCK_2"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
    {
      std::string sideblockName("SHARED_SURF_BLOCK_1_" + get_topology_name("LINE_2"));
      sideblockName = Ioss::Utils::uppercase(sideblockName);

      Ioss::SideBlock *sideblock = shared_surf->get_side_block(sideblockName);
      EXPECT_TRUE(nullptr != sideblock);

      std::vector<std::string> touchingBlocks;
      db_io->compute_block_membership(sideblock, touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_1"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
    {
      std::string sideblockName("SHARED_SURF_BLOCK_2_" + get_topology_name("LINE_2"));
      sideblockName = Ioss::Utils::uppercase(sideblockName);

      Ioss::SideBlock *sideblock = shared_surf->get_side_block(sideblockName);
      EXPECT_TRUE(nullptr != sideblock);

      std::vector<std::string> touchingBlocks;
      db_io->compute_block_membership(sideblock, touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_2"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }

    Ioss::SideSet *owned_surf = region.get_sideset("OWNED_SURF");
    EXPECT_TRUE(nullptr != owned_surf);
    {
      std::vector<std::string> touchingBlocks;
      owned_surf->block_membership(touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_2"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
    {
      std::string sideblockName("OWNED_SURF_BLOCK_2_" + get_topology_name("LINE_2"));
      sideblockName = Ioss::Utils::uppercase(sideblockName);

      Ioss::SideBlock *sideblock = owned_surf->get_side_block(sideblockName);
      EXPECT_TRUE(nullptr != sideblock);

      std::vector<std::string> touchingBlocks;
      db_io->compute_block_membership(sideblock, touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_2"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
  }

  TEST(TextMesh, surfaceToBlockMapping_splitByTopology)
  {
    if (get_parallel_size() != 1) {
      GTEST_SKIP();
    }

    std::string meshDesc = "0,1,PYRAMID_5,1,2,3,4,5,block_1\n"
                           "0,2,TET_4,2,3,5,6,block_2"
                           "|sideset:name=surface_1; skin=all; split=topology";

    Iotm::DatabaseIO *db_io = create_input_db_io(meshDesc);

    Ioss::Region region(db_io);

    EXPECT_TRUE(nullptr != db_io);
    EXPECT_TRUE(db_io->ok());
    EXPECT_EQ("TextMesh", db_io->get_format());

    const std::vector<Ioss::ElementBlock *> &element_blocks = region.get_element_blocks();
    EXPECT_EQ(2u, element_blocks.size());

    const std::vector<Ioss::SideSet *> &sidesets = region.get_sidesets();
    EXPECT_EQ(1u, sidesets.size());

    Iotm::IossTopologyMapping topologyMapping;
    topologyMapping.initialize_topology_map();

    auto get_topology_name = [&topologyMapping](const std::string &textMeshTopologyName) {
      return topologyMapping.topology(textMeshTopologyName).name();
    };

    Ioss::SideSet *surf_1 = region.get_sideset("SURFACE_1");
    EXPECT_TRUE(nullptr != surf_1);
    {
      std::vector<std::string> touchingBlocks;
      surf_1->block_membership(touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_1", "BLOCK_2"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
    {
      std::string sideblockName("SURFACE_" + get_topology_name("PYRAMID_5") + "_" +
                                get_topology_name("QUAD_4") + "_1");
      sideblockName = Ioss::Utils::uppercase(sideblockName);

      Ioss::SideBlock *sideblock = surf_1->get_side_block(sideblockName);
      EXPECT_TRUE(nullptr != sideblock);

      std::vector<std::string> touchingBlocks;
      db_io->compute_block_membership(sideblock, touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_1"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
    {
      std::string sideblockName("SURFACE_" + get_topology_name("PYRAMID_5") + "_" +
                                get_topology_name("TRI_3") + "_1");
      sideblockName = Ioss::Utils::uppercase(sideblockName);

      Ioss::SideBlock *sideblock = surf_1->get_side_block(sideblockName);
      EXPECT_TRUE(nullptr != sideblock);

      std::vector<std::string> touchingBlocks;
      db_io->compute_block_membership(sideblock, touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_1"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
    {
      std::string sideblockName("SURFACE_" + get_topology_name("TET_4") + "_" +
                                get_topology_name("TRI_3") + "_1");
      sideblockName = Ioss::Utils::uppercase(sideblockName);

      Ioss::SideBlock *sideblock = surf_1->get_side_block(sideblockName);
      EXPECT_TRUE(nullptr != sideblock);

      std::vector<std::string> touchingBlocks;
      db_io->compute_block_membership(sideblock, touchingBlocks);

      std::vector<std::string> goldTouchingBlocks{"BLOCK_2"};
      EXPECT_EQ(goldTouchingBlocks, touchingBlocks);
    }
  }

} // namespace
