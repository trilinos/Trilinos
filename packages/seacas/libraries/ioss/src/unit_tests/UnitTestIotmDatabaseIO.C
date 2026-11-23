// Copyright(C) 1999-2020, 2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"
#include <chrono>
#include <functional>
#include <future>

#include "Ionit_Initializer.h"
#include "Ioss_CopyDatabase.h"
#include "Ioss_DBUsage.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_Hex8.h"
#include "Ioss_MeshCopyOptions.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_Property.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_Shell4.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "text_mesh/Iotm_DatabaseIO.h"
#include "text_mesh/Iotm_TextMeshTopologyMapping.h"
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "Ioss_CodeTypes.h"
#include "Ioss_Field.h"
#include "Ioss_ParallelUtils.h"

namespace {

  int db_api_int_size(Ioss::DatabaseIO *db)
  {
    assert(db != nullptr);
    return db->int_byte_size_api();
  }

  template <typename INT>
  std::vector<int64_t> get_element_ids_from_block_impl(const Ioss::ElementBlock *block)
  {
    std::vector<int64_t> elemIds;
    std::vector<INT>     ids;

    block->get_field_data("ids", ids);

    for (INT id : ids) {
      elemIds.push_back(static_cast<int64_t>(id));
    }

    return elemIds;
  }

  std::vector<int64_t> get_element_ids_from_block(const Ioss::ElementBlock *block)
  {
    if (db_api_int_size(block->get_database()) == 4) {
      return get_element_ids_from_block_impl<int>(block);
    }
    else {
      return get_element_ids_from_block_impl<int64_t>(block);
    }
  }

  template <typename INT>
  bool get_element_conn_from_block_impl(int64_t elemId, const Ioss::ElementBlock *block,
                                        std::vector<int64_t> &elemConn)
  {
    const Ioss::ElementTopology *topo = nullptr;

    std::vector<INT> connectivity;
    std::vector<INT> elemIds;

    block->get_field_data("ids", elemIds);
    block->get_field_data("connectivity", connectivity);

    topo = block->topology();

    size_t elementCount = elemIds.size();
    int    nodesPerElem = topo->number_nodes();

    for (size_t i = 0; i < elementCount; ++i) {
      INT *conn = &connectivity[i * nodesPerElem];
      auto id   = static_cast<int64_t>(elemIds[i]);

      if (id == elemId) {
        for (int j = 0; j < nodesPerElem; j++) {
          elemConn.push_back(conn[j]);
        }

        return true;
      }
    }

    return false;
  }

  template <typename INT>
  std::vector<int64_t> get_element_conn_impl(Ioss::Region &region, int64_t elemId)
  {
    std::vector<int64_t> elemConn;

    const Ioss::ElementBlockContainer &elemBlocks = region.get_element_blocks();

    for (const Ioss::ElementBlock *block : elemBlocks) {
      if (get_element_conn_from_block_impl<INT>(elemId, block, elemConn)) {
        return elemConn;
      }
    }
    return std::vector<int64_t>{};
  }

  std::vector<int64_t> get_element_conn(Ioss::Region &region, int64_t elemId)
  {
    std::vector<int64_t> elemConn;

    if (db_api_int_size(region.get_database()) == 4) {
      return get_element_conn_impl<int>(region, elemId);
    }
    else {
      return get_element_conn_impl<int64_t>(region, elemId);
    }
  }

  void define_element_transient(Ioss::Region &o_region, const std::string &elemFieldName)
  {
    o_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

    for (Ioss::ElementBlock *o_eb : o_region.get_element_blocks()) {
      size_t      num_elem = o_eb->entity_count();
      std::string storage  = "scalar";

      Ioss::Field field(elemFieldName, Ioss::Field::REAL, storage, 1, Ioss::Field::Field::TRANSIENT,
                        num_elem);
      o_eb->field_add(field);
    }
    o_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  }

  void write_element_transient(Ioss::Region &o_region, const std::string &elemFieldName)
  {
    int numTimeSteps = o_region.get_implicit_property("state_count").get_int();

    o_region.begin_mode(Ioss::STATE_TRANSIENT);
    int step = o_region.add_state((double)numTimeSteps);
    o_region.begin_state(step);

    for (Ioss::ElementBlock *o_eb : o_region.get_element_blocks()) {
      size_t num_elem = o_eb->entity_count();

      std::vector<double>  field_data(num_elem);
      std::vector<int64_t> elem_ids = get_element_ids_from_block(o_eb);

      for (size_t i = 0; i < elem_ids.size(); i++) {
        field_data[i] = (double)elem_ids[i];
      }

      o_eb->put_field_data(elemFieldName, field_data);
    }

    o_region.end_state(step);
    o_region.end_mode(Ioss::STATE_TRANSIENT);
  }

  Iotm::DatabaseIO *create_input_db_io(const std::string &meshDesc,
                                       Ioss_MPI_Comm      comm = Ioss::ParallelUtils::comm_world())
  {
    Ioss::Init::Initializer init_db;

    Ioss::DatabaseUsage   db_usage = Ioss::READ_MODEL;
    Ioss::PropertyManager properties;

    properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));

    auto *db_io = new Iotm::DatabaseIO(nullptr, meshDesc, db_usage, comm, properties);
    return db_io;
  }

  Ioss::DatabaseIO *create_output_db_io(const std::string &outputFile,
                                        Ioss_MPI_Comm      comm = Ioss::ParallelUtils::comm_world())
  {
    Ioss::DatabaseUsage   db_usage = Ioss::WRITE_RESTART;
    Ioss::PropertyManager properties;

    properties.add(Ioss::Property("FLUSH_INTERVAL", 1));
    properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));

    Ioss::DatabaseIO *db_io =
        Ioss::IOFactory::create("exodusII", outputFile, db_usage, comm, properties);
    return db_io;
  }

  int get_parallel_size(Ioss_MPI_Comm comm = Ioss::ParallelUtils::comm_world())
  {
    return Ioss::ParallelUtils(comm).parallel_size();
  }

  int get_parallel_rank(Ioss_MPI_Comm comm = Ioss::ParallelUtils::comm_world())
  {
    return Ioss::ParallelUtils(comm).parallel_rank();
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
    EXPECT_TRUE(nullptr != db_io);
    db_io->set_surface_split_type(Ioss::SPLIT_BY_ELEMENT_BLOCK);

    Ioss::Region region(db_io);

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

  TEST(TextMesh, inputDuplicateMeshCommSelf)
  {
    if (get_parallel_size() != 2) {
      GTEST_SKIP();
    }

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1";

    Iotm::DatabaseIO *db_io = create_input_db_io(meshDesc, Ioss::ParallelUtils::comm_self());
    EXPECT_TRUE(nullptr != db_io);

    Ioss::Region region(db_io);

    EXPECT_TRUE(db_io->ok());
    EXPECT_EQ("TextMesh", db_io->get_format());

    // With COMM_SELF, the mesh should be duplicated on each MPI rank
    const std::vector<Ioss::ElementBlock *> &element_blocks = region.get_element_blocks();
    EXPECT_EQ(1u, element_blocks.size());
    EXPECT_EQ(1u, element_blocks[0]->entity_count());
    EXPECT_EQ(Ioss::Hex8::name, element_blocks[0]->topology()->name());

    const Ioss::NodeBlockContainer &node_blocks = region.get_node_blocks();
    EXPECT_EQ(1u, node_blocks.size());
    EXPECT_EQ(8u, node_blocks[0]->entity_count());

    std::vector<int64_t> elemConn = get_element_conn(region, 1);
    std::vector<int64_t> goldNodeIds{1, 2, 3, 4, 5, 6, 7, 8};
    EXPECT_EQ(goldNodeIds, elemConn);
  }

  TEST(TextMesh, outputMeshCommSplit)
  {
    if (get_parallel_size() != 4) {
      GTEST_SKIP();
    }

    Ioss_MPI_Comm comm = Ioss::ParallelUtils::comm_self();
#ifdef SEACAS_HAVE_MPI
    int color = get_parallel_rank() % 2;
    int key   = get_parallel_rank();
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
#endif

    std::string       meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1"
                                 "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1";
    Iotm::DatabaseIO *db_i     = create_input_db_io(meshDesc, comm);
    ASSERT_FALSE(db_i == nullptr || !db_i->ok(true));

    Ioss::Region region_i(db_i, "region_i");

    std::string outputFile;
    {
      std::ostringstream os;
      os << "output_file.e.";
      os << get_parallel_rank();

      outputFile = os.str();
    }

    Ioss::DatabaseIO *db_o = create_output_db_io(outputFile, comm);
    ASSERT_FALSE(db_o == nullptr || !db_o->ok(true));

    // NOTE: 'region_o' owns 'db_o' pointer at this time
    Ioss::Region region_o(db_o, "region_o");

    Ioss::MeshCopyOptions options{};
    options.verbose           = true;
    options.output_summary    = true;
    options.debug             = false;
    options.ints_64_bit       = false;
    options.data_storage_type = 1;
    options.add_proc_id       = false;
    Ioss::copy_database(region_i, region_o, options);

    const std::string elemFieldName = "elem_id_data";
    define_element_transient(region_o, elemFieldName);
    write_element_transient(
        region_o, elemFieldName); // Remove this line and enable the pre-processor macro below

    unlink(db_o->decoded_filename().c_str());
  }
} // namespace
