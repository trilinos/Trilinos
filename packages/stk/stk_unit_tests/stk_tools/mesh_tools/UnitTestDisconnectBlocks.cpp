#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocksImpl.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/environment/WallTime.hpp>
#include "stk_unit_test_utils/getOption.h"
#include <string>

namespace {

void create_sides_between_blocks(stk::mesh::BulkData& bulk,
                                 const std::string& block1Name,
                                 const std::string& block2Name,
                                 const std::string& sidePartName)
{
    stk::mesh::Part& block1 = *bulk.mesh_meta_data().get_part(block1Name);
    stk::mesh::Part& block2 = *bulk.mesh_meta_data().get_part(block2Name);
    stk::mesh::Part& sidePart = *bulk.mesh_meta_data().get_part(sidePartName);

    stk::mesh::Selector blockSelector = block1 | block2;
    stk::mesh::create_interior_block_boundary_sides(bulk, blockSelector, stk::mesh::PartVector{&sidePart});
}

void create_sideset(stk::mesh::BulkData& bulk,
                    const std::string& surfacePartName,
                    const std::string& blockPartName)
{
    stk::mesh::Part& blockPart = *bulk.mesh_meta_data().get_part(blockPartName);
    stk::mesh::Part& surfacePart = *bulk.mesh_meta_data().get_part(surfacePartName);

    bulk.mesh_meta_data().set_surface_to_block_mapping(&surfacePart,
                                                       stk::mesh::ConstPartVector{&blockPart});
    bulk.create_sideset(surfacePart);
}

void move_elems_from_block_to_block(stk::mesh::BulkData& bulk,
                                    const std::vector<stk::mesh::EntityId>& elemIDs,
                                    const std::string& fromBlockName,
                                    const std::string& toBlockName)
{
    stk::mesh::Part& fromBlock = *bulk.mesh_meta_data().get_part(fromBlockName);
    stk::mesh::Part& toBlock = *bulk.mesh_meta_data().get_part(toBlockName);

    stk::mesh::EntityVector elems;
    for(stk::mesh::EntityId elemID : elemIDs) {
        stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemID);
        ThrowRequireMsg(bulk.is_valid(elem), "Failed to find element with ID="<<elemID);
        elems.push_back(elem);
    }

    bulk.batch_change_entity_parts(elems, stk::mesh::PartVector{&toBlock}, stk::mesh::PartVector{&fromBlock});
}

unsigned get_num_common_nodes(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks)
{
  stk::mesh::Selector selector;

  if (blocks.size() == 1u) {
    selector = *blocks[0];
  }
  else {
    for (size_t i = 0; i < blocks.size()-1; ++i) {
      for (size_t j = i+1; j < blocks.size(); ++j) {
        selector |= (*blocks[i] & *blocks[j]);
      }
    }
  }

  for (stk::mesh::Bucket * bucket : bulk.get_buckets(stk::topology::NODE_RANK, selector)) {
    for (stk::mesh::Entity node : *bucket) {
      std::cout << "Common node: " << bulk.entity_key(node) << std::endl;
    }
  }

  return stk::mesh::count_selected_entities(selector, bulk.buckets(stk::topology::NODE_RANK));
}

unsigned get_num_total_nodes(const stk::mesh::BulkData & bulk)
{
  return stk::mesh::count_selected_entities(bulk.mesh_meta_data().universal_part(), bulk.buckets(stk::topology::NODE_RANK));
}

void check_disconnected_nodes(stk::mesh::BulkData & bulk)
{
  for (stk::mesh::Bucket * bucket : bulk.buckets(stk::topology::NODE_RANK)) {
    for (stk::mesh::Entity node : *bucket) {
      const unsigned numElems = bulk.num_elements(node);
      if (numElems == 0u) {
        std::cout << "Found disconnected node: " << bulk.entity_key(node) << std::endl;
      }
    }
  }
}

void output_mesh(const std::string & fileName, stk::mesh::BulkData & bulk)
{
  std::string writeOutput = stk::unit_test_util::get_option("--output", "off");
  if (writeOutput == "on") {
    stk::io::write_mesh("disconnect_1block_1quad.g", bulk);
  }
}

}//namespace


class TestDisconnectBlocks2D : public stk::unit_test_util::MeshFixture
{
protected:
    TestDisconnectBlocks2D() : stk::unit_test_util::MeshFixture(2)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    }

    stk::mesh::Part & create_part(const stk::topology topology, const std::string & blockName, int64_t blockId)
    {
      stk::mesh::MetaData& meta = get_meta();
      stk::mesh::Part& part = meta.declare_part_with_topology(blockName, topology);
      stk::io::put_io_part_attribute(part);
      meta.set_part_id(part, blockId);
      return part;
    }

    stk::mesh::PartVector setup_mesh_1block_1quad()
    {
      stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1";
      std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector block = block1;

      unsigned expectedNumCommonNodes = 4;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(block, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1};
    }

    stk::mesh::PartVector setup_mesh_2block_1quad()
    {
      stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

      std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1";
      std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector bothBlocks = block1 & block2;
      stk::mesh::Selector eitherBlocks = block1 | block2;

      unsigned expectedNumCommonNodes = 0;
      unsigned expectedNumTotalNodes  = 4;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));
      EXPECT_EQ(expectedNumTotalNodes,  stk::mesh::count_selected_entities(eitherBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector setup_mesh_2block_2quad()
    {
      stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

      std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                             "0,2,QUAD_4_2D,2,3,6,5,block_2";
      std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector bothBlocks = block1 & block2;

      unsigned expectedNumCommonNodes = 2;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector setup_mesh_2block_4quad_corner()
    {
      stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

      std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                             "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                             "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                             "0,4,QUAD_4_2D,5,6,9,8,block_1";
      std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector bothBlocks = block1 & block2;

      unsigned expectedNumCommonNodes = 3;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector setup_mesh_2block_4quad_swappedCorner()
    {
      stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

      std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                             "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                             "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                             "0,4,QUAD_4_2D,5,6,9,8,block_2";
      std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector bothBlocks = block1 & block2;

      unsigned expectedNumCommonNodes = 3;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector create_3_blocks_order1()
    {
      stk::mesh::PartVector blocks(3);
      blocks[0] = &create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      blocks[1] = &create_part(stk::topology::QUAD_4_2D, "block_2", 2);
      blocks[2] = &create_part(stk::topology::QUAD_4_2D, "block_3", 3);
      return blocks;
    }

    stk::mesh::PartVector create_3_blocks_order2()
    {
      stk::mesh::PartVector blocks(3);
      blocks[0] = &create_part(stk::topology::QUAD_4_2D, "block_2", 2);
      blocks[1] = &create_part(stk::topology::QUAD_4_2D, "block_3", 3);
      blocks[2] = &create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      return blocks;
    }

    stk::mesh::PartVector create_3_blocks_order3()
    {
      stk::mesh::PartVector blocks(3);
      blocks[0] = &create_part(stk::topology::QUAD_4_2D, "block_3", 3);
      blocks[1] = &create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      blocks[2] = &create_part(stk::topology::QUAD_4_2D, "block_2", 2);
      return blocks;
    }

    stk::mesh::PartVector create_3_blocks_order4()
    {
      stk::mesh::PartVector blocks(3);
      blocks[0] = &create_part(stk::topology::QUAD_4_2D, "block_3", 3);
      blocks[1] = &create_part(stk::topology::QUAD_4_2D, "block_2", 2);
      blocks[2] = &create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      return blocks;
    }

    stk::mesh::PartVector create_3_blocks_order5()
    {
      stk::mesh::PartVector blocks(3);
      blocks[0] = &create_part(stk::topology::QUAD_4_2D, "block_2", 2);
      blocks[1] = &create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      blocks[2] = &create_part(stk::topology::QUAD_4_2D, "block_3", 3);
      return blocks;
    }

    stk::mesh::PartVector create_3_blocks_order6()
    {
      stk::mesh::PartVector blocks(3);
      blocks[0] = &create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      blocks[1] = &create_part(stk::topology::QUAD_4_2D, "block_3", 3);
      blocks[2] = &create_part(stk::topology::QUAD_4_2D, "block_2", 2);
      return blocks;
    }

    void setup_mesh_3block_4quad_base(const stk::mesh::PartVector & blocks)
    {
      std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                             "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                             "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                             "0,4,QUAD_4_2D,5,6,9,8,block_3";
      std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector sharedByBlocks = (*blocks[0] & *blocks[1]) | (*blocks[0] & *blocks[2]) | (*blocks[1] & *blocks[2]);

      unsigned expectedNumCommonNodes = 4;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));
    }

    stk::mesh::PartVector setup_mesh_3block_4quad(const unsigned order)
    {
      stk::mesh::PartVector blocks;
      if (order == 1) {
        blocks = create_3_blocks_order1();
      } else if (order == 2) {
        blocks = create_3_blocks_order2();
      } else if (order == 3) {
        blocks = create_3_blocks_order3();
      } else if (order == 4) {
        blocks = create_3_blocks_order4();
      } else if (order == 5) {
        blocks = create_3_blocks_order5();
      } else if (order == 6) {
        blocks = create_3_blocks_order6();
      } else {
        std::cerr << "ERROR: Unexpected part ordinal ordering!!!" << std::endl;
        exit(1);
      }

      setup_mesh_3block_4quad_base(blocks);
      return blocks;
    }

    void test_mesh_3block_4quad(const unsigned order) {
      if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

      stk::mesh::PartVector blocks = setup_mesh_3block_4quad(order);

      stk::io::write_mesh("disconnect_3block_4quad_order" + std::to_string(order) + "_init.g", get_bulk());

      stk::tools::disconnect_all_blocks(get_bulk());

      EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
      EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));

      stk::io::write_mesh("disconnect_3block_4quad_order" + std::to_string(order) + ".g", get_bulk());
    }

    stk::mesh::PartVector setup_mesh_3block_4quad_keepLowerRight()
    {
      stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);
      stk::mesh::Part & block3 = create_part(stk::topology::QUAD_4_2D, "block_3", 3);

      std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                             "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                             "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                             "0,4,QUAD_4_2D,5,6,9,8,block_2";
      std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector sharedByBlocks = (block1 & block2) | (block1 & block3) | (block2 & block3);

      unsigned expectedNumCommonNodes = 4;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2, &block3};
    }

    stk::mesh::PartVector setup_mesh_2block_4quad_checkerboard()
    {
      stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

      std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                             "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                             "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                             "0,4,QUAD_4_2D,5,6,9,8,block_1";
      std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector sharedByBlocks = (block1 & block2);

      unsigned expectedNumCommonNodes = 5;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector setup_mesh_3block_4quad_checkerboard()
    {
      stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);
      stk::mesh::Part & block3 = create_part(stk::topology::QUAD_4_2D, "block_3", 3);

      std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                             "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
                             "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                             "0,4,QUAD_4_2D,5,6,9,8,block_1";
      std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector sharedByBlocks = (block1 & block2) | (block1 & block3) | (block2 & block3);

      unsigned expectedNumCommonNodes = 5;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2, &block3};
    }

    stk::mesh::PartVector setup_mesh_2block_2quad_diagonal()
    {
      stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
      stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

      std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                             "0,2,QUAD_4_2D,4,5,7,6,block_2";
      std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 1,2, 2,2 };
      stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

      stk::mesh::Selector sharedByBlocks = (block1 & block2);

      unsigned expectedNumCommonNodes = 1;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

};

TEST_F(TestDisconnectBlocks2D, disconnect_1block_1quad)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_1block_1quad();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(4u, get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(4u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_1block_1quad.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_1quad)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_1quad();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u, get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(4u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_1quad.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_2quad();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u, get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_2quad.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner();

        output_mesh("disconnect_2block_4quad_corner_init.g", get_bulk());

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_4quad_corner.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner();

        output_mesh("disconnect_2block_4quad_swappedCorner_init.g", get_bulk());

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_4quad_swappedCorner.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_order1)
{
    test_mesh_3block_4quad(1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_order2)
{
    test_mesh_3block_4quad(2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_order3)
{
    test_mesh_3block_4quad(3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_order4)
{
    test_mesh_3block_4quad(4);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_order5)
{
    test_mesh_3block_4quad(5);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_order6)
{
    test_mesh_3block_4quad(6);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight();

        output_mesh("disconnect_3block_4quad_keepLowerRight_init.g", get_bulk());

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_3block_4quad_keepLowerRight.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard();

        output_mesh("disconnect_2block_4quad_checkerboard_init.g", get_bulk());

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_4quad_checkerboard.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard();

        output_mesh("disconnect_3block_4quad_checkerboard_init.g", get_bulk());

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_3block_4quad_checkerboard.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad_diagonal)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_2quad_diagonal();

        output_mesh("disconnect_2block_2quad_diagonal_init.g", get_bulk());

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_2quad_diagonal.g", get_bulk());
    }
}


class TestDisconnectBlocks : public stk::unit_test_util::MeshFixture
{
protected:
    stk::mesh::Part & create_part(const stk::topology topology, const std::string & blockName, int64_t blockId)
    {
      stk::mesh::MetaData& meta = get_meta();
      stk::mesh::Part& part = meta.declare_part_with_topology(blockName, topology);
      stk::io::put_io_part_attribute(part);
      meta.set_part_id(part, blockId);
      return part;
    }

    stk::mesh::PartVector setup_mesh_1block_1hex()
    {
      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x1x1", get_bulk());

      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      unsigned expectedNumCommonNodes = 8;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(block1, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1};
    }

    stk::mesh::PartVector setup_mesh_2block_1hex()
    {
      stk::mesh::Part & block2 = create_part(stk::topology::HEX_8, "block_2", 2);
      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x1x1", get_bulk());
      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      stk::mesh::Selector bothBlocks = block1 & block2;
      stk::mesh::Selector eitherBlocks = block1 | block2;

      unsigned expectedNumCommonNodes = 0;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      unsigned expectedNumTotalNodes = 8;
      EXPECT_EQ(expectedNumTotalNodes, stk::mesh::count_selected_entities(eitherBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector setup_mesh_2block_2hex()
    {
      stk::mesh::Part & block2 = create_part(stk::topology::HEX_8, "block_2", 2);
      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x1x2", get_bulk());
      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

      stk::mesh::Selector bothBlocks = block1 & block2;

      unsigned expectedNumCommonNodes = 4;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector setup_mesh_2block_2hex_withSides()
    {
      stk::mesh::Part & block2 = create_part(stk::topology::HEX_8, "block_2", 2);
      create_part(stk::topology::QUAD_4, "surface_1", 1);

      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x1x2", get_bulk());
      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

      create_sideset(get_bulk(), "surface_1", "block_1");
      create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");

      stk::mesh::Selector bothBlocks = block1 & block2;

      unsigned expectedNumCommonSides = 1;
      unsigned expectedNumCommonNodes = 4;
      EXPECT_EQ(expectedNumCommonSides, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::FACE_RANK)));
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector setup_mesh_3block_4hex()
    {
      stk::mesh::Part & block2 = create_part(stk::topology::HEX_8, "block_2", 2);
      stk::mesh::Part & block3 = create_part(stk::topology::HEX_8, "block_3", 3);

      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x2x2", get_bulk());
      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{4}, "block_1", "block_3");

      stk::mesh::Selector blocks_1_and_2 = block1 & block2;
      stk::mesh::Selector blocks_1_and_3 = block1 & block3;
      stk::mesh::Selector blocks_2_and_3 = block2 & block3;

      stk::mesh::Selector sharedByBlocks = blocks_1_and_2 | blocks_1_and_3 | blocks_2_and_3;

      unsigned expectedNumCommonNodes = 8;
      unsigned nodesSharedByMultipleBlocks = stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::NODE_RANK));

      EXPECT_EQ(expectedNumCommonNodes, nodesSharedByMultipleBlocks);

      return {&block1, &block2, &block3};
    }

    stk::mesh::PartVector setup_mesh_3block_4hex_withSides()
    {
      stk::mesh::Part & block2 = create_part(stk::topology::HEX_8, "block_2", 2);
      stk::mesh::Part & block3 = create_part(stk::topology::HEX_8, "block_3", 3);
      create_part(stk::topology::QUAD_4, "surface_1", 1);
      create_part(stk::topology::QUAD_4, "surface_2", 2);

      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x2x2", get_bulk());
      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{4}, "block_1", "block_3");

      create_sideset(get_bulk(), "surface_1", "block_1");
      create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");
      create_sides_between_blocks(get_bulk(), "block_1", "block_3", "surface_1");
      create_sideset(get_bulk(), "surface_2", "block_2");
      create_sides_between_blocks(get_bulk(), "block_2", "block_3", "surface_2");

      stk::mesh::Selector blocks_1_and_2 = block1 & block2;
      stk::mesh::Selector blocks_1_and_3 = block1 & block3;
      stk::mesh::Selector blocks_2_and_3 = block2 & block3;
      stk::mesh::Selector sharedByBlocks = blocks_1_and_2 | blocks_1_and_3 | blocks_2_and_3;

      unsigned expectedNumCommonSides = 3;
      unsigned expectedNumCommonNodes = 8;
      unsigned sidesSharedByMultipleBlocks = stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::FACE_RANK));
      unsigned nodesSharedByMultipleBlocks = stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::NODE_RANK));

      EXPECT_EQ(expectedNumCommonSides, sidesSharedByMultipleBlocks);
      EXPECT_EQ(expectedNumCommonNodes, nodesSharedByMultipleBlocks);

      return {&block1, &block2, &block3};
    }

    stk::mesh::PartVector setup_mesh_2block_2cubeOfTet()
    {
      stk::mesh::Part & block2 = create_part(stk::topology::TET_4, "block_2", 2);

      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x1x2|tets", get_bulk());
      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12},
                                     "block_1", "block_2");

      stk::mesh::Selector bothBlocks = block1 & block2;

      unsigned expectedNumCommonNodes = 4;
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector setup_mesh_2block_2cubeOfTet_withSides()
    {
      stk::mesh::Part & block2 = create_part(stk::topology::TET_4, "block_2", 2);
      create_part(stk::topology::TRI_3, "surface_1", 1);

      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x1x2|tets", get_bulk());
      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12},
                                     "block_1", "block_2");

      create_sideset(get_bulk(), "surface_1", "block_1");
      create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");

      stk::mesh::Selector bothBlocks = block1 & block2;

      unsigned expectedNumCommonSides = 2;
      unsigned expectedNumCommonNodes = 4;
      EXPECT_EQ(expectedNumCommonSides, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::FACE_RANK)));
      EXPECT_EQ(expectedNumCommonNodes, stk::mesh::count_selected_entities(bothBlocks, get_bulk().buckets(stk::topology::NODE_RANK)));

      return {&block1, &block2};
    }

    stk::mesh::PartVector setup_mesh_3block_4cubeOfTet()
    {
      stk::mesh::Part & block2 = create_part(stk::topology::TET_4, "block_2", 2);
      stk::mesh::Part & block3 = create_part(stk::topology::TET_4, "block_3", 3);

      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x2x2|tets", get_bulk());
      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12}, "block_1", "block_2");
      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{19, 20, 21, 22, 23, 24}, "block_1", "block_3");

      stk::mesh::Selector blocks_1_and_2 = block1 & block2;
      stk::mesh::Selector blocks_1_and_3 = block1 & block3;
      stk::mesh::Selector blocks_2_and_3 = block2 & block3;

      stk::mesh::Selector sharedByBlocks = blocks_1_and_2 | blocks_1_and_3 | blocks_2_and_3;

      unsigned expectedNumCommonNodes = 8;
      unsigned nodesSharedByMultipleBlocks = stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::NODE_RANK));

      EXPECT_EQ(expectedNumCommonNodes, nodesSharedByMultipleBlocks);

      return {&block1, &block2, &block3};
    }

    stk::mesh::PartVector setup_mesh_3block_4cubeOfTet_withSides()
    {
      stk::mesh::Part & block2 = create_part(stk::topology::TET_4, "block_2", 2);
      stk::mesh::Part & block3 = create_part(stk::topology::TET_4, "block_3", 3);
      create_part(stk::topology::TRI_3, "surface_1", 1);
      create_part(stk::topology::TRI_3, "surface_2", 2);

      allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::io::fill_mesh("generated:1x2x2|tets", get_bulk());
      stk::mesh::Part & block1 = *get_meta().get_part("block_1");

      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12}, "block_1", "block_2");
      move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{19, 20, 21, 22, 23, 24}, "block_1", "block_3");

      create_sideset(get_bulk(), "surface_1", "block_1");
      create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");
      create_sides_between_blocks(get_bulk(), "block_1", "block_3", "surface_1");
      create_sideset(get_bulk(), "surface_2", "block_2");
      create_sides_between_blocks(get_bulk(), "block_2", "block_3", "surface_2");

      stk::mesh::Selector blocks_1_and_2 = block1 & block2;
      stk::mesh::Selector blocks_1_and_3 = block1 & block3;
      stk::mesh::Selector blocks_2_and_3 = block2 & block3;

      stk::mesh::Selector sharedByBlocks = blocks_1_and_2 | blocks_1_and_3 | blocks_2_and_3;

      unsigned expectedNumCommonSides = 6;
      unsigned expectedNumCommonNodes = 8;
      unsigned sidesSharedByMultipleBlocks = stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::FACE_RANK));
      unsigned nodesSharedByMultipleBlocks = stk::mesh::count_selected_entities(sharedByBlocks, get_bulk().buckets(stk::topology::NODE_RANK));

      EXPECT_EQ(expectedNumCommonSides, sidesSharedByMultipleBlocks);
      EXPECT_EQ(expectedNumCommonNodes, nodesSharedByMultipleBlocks);

      return {&block1, &block2, &block3};
    }

};

TEST_F(TestDisconnectBlocks, disconnect_1block_1hex)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_1block_1hex();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(8u, get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_1block_1hex.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_1hex)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_1hex();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u, get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_1hex.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2hex)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_2hex();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_2hex.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_2block_2hex_withSides)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_2hex_withSides();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_2hex_withSides.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2cubeOfTet)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_2cubeOfTet.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_2block_2cubeOfTet_withSides)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet_withSides();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_2block_2cubeOfTet_withSides.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks, disconnect_3block_4hex)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_3block_4hex();

        output_mesh("disconnect_3block_4hex_init.g", get_bulk());

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
        check_disconnected_nodes(get_bulk());

        output_mesh("disconnect_3block_4hex.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_3block_4hex_withSides)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_3block_4hex_withSides();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_3block_4hex_withSides.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks, disconnect_3block_4cubeOfTet)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
        check_disconnected_nodes(get_bulk());

        output_mesh("disconnect_3block_4cubeOfTet.g", get_bulk());
    }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_3block_4cubeOfTet_withSides)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet_withSides();

        stk::tools::disconnect_all_blocks(get_bulk());

        EXPECT_EQ(0u,  get_num_common_nodes(get_bulk(), blocks));
        EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));

        output_mesh("disconnect_3block_4cubeOfTet_withSides.g", get_bulk());
    }
}

TEST(DisconnectBlocks, input_mesh)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
        double startTime = stk::wall_time();

        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);

        std::string exodusFileName = stk::unit_test_util::get_option("-f", "");
        if (exodusFileName.empty()) return;

        stk::io::fill_mesh(exodusFileName, bulk);

        double meshReadTime = stk::wall_time();

        stk::tools::disconnect_all_blocks(bulk);

        double disconnectTime = stk::wall_time();

        stk::mesh::PartVector allBlocksInMesh;
        stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
        EXPECT_EQ(0u,  get_num_common_nodes(bulk, allBlocksInMesh));

        check_disconnected_nodes(bulk);

        stk::io::write_mesh("disconnected_" + exodusFileName, bulk);
        double meshWriteTime = stk::wall_time();

        std::cout << " Mesh read time = " << (meshReadTime - startTime) << " s" << std::endl;
        std::cout << "Disconnect time = " << (disconnectTime - meshReadTime) << " s" << std::endl;
        std::cout << "Mesh write time = " << (meshWriteTime - disconnectTime) << " s" << std::endl;

    }
}

