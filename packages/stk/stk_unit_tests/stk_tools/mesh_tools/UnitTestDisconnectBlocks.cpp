#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_unit_test_utils/getOption.h"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include <algorithm>
#include <fstream>
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <libgen.h>
#include <stk_io/FillMesh.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_tools/mesh_tools/DetectHingesImpl.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocksImpl.hpp>
#include <stk_tools/mesh_tools/DisconnectUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <string>

namespace {

struct BlockConnection {
  BlockConnection(const std::string& b1, const std::string& b2)
  : block1(b1), block2(b2), numExpectedIntersectingNodes(0) { }

  BlockConnection(const std::string& b1, const std::string& b2, unsigned expectedNumNodes)
  : block1(b1), block2(b2), numExpectedIntersectingNodes(expectedNumNodes) { }

  std::string block1;
  std::string block2 = 0;
  unsigned numExpectedIntersectingNodes = 0;
};

typedef std::vector<BlockConnection> BlockConnectionVector;

stk::tools::BlockPairVector convert_connection_vector_to_pair_vector(const stk::mesh::BulkData& bulk, const BlockConnectionVector& disconnectConnVector)
{
  stk::tools::BlockPairVector pairVector;
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::tools::impl::PartPairLess comparator;

  for(const BlockConnection& connection : disconnectConnVector) {
    stk::mesh::Part* part1 = meta.get_part(connection.block1);
    stk::mesh::Part* part2 = meta.get_part(connection.block2);

    stk::util::insert_keep_sorted_and_unique(stk::tools::impl::get_block_pair(part1, part2), pairVector, comparator);
  }

  return pairVector;
}

stk::tools::BlockPairVector get_local_reconnect_list(const stk::mesh::BulkData& bulk, const BlockConnectionVector& disconnectList)
{
  stk::tools::BlockPairVector convertedDisconnectList = convert_connection_vector_to_pair_vector(bulk, disconnectList);

  return stk::tools::impl::get_local_reconnect_list(bulk, convertedDisconnectList);
}

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

void create_all_boundary_sides(stk::mesh::BulkData& bulk, const std::string& sidePartName)
{
  stk::mesh::Part& sidePart = *bulk.mesh_meta_data().get_part(sidePartName);
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);

  stk::mesh::Selector blockSelector = stk::mesh::selectUnion(allBlocksInMesh);
  stk::mesh::create_exposed_block_boundary_sides(bulk, blockSelector, stk::mesh::PartVector{&sidePart});
}

unsigned get_num_surface_nodes(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blockPartNames)
{
  stk::mesh::PartVector blockParts;
  stk::mesh::EntityRank rank = bulk.mesh_meta_data().side_rank();
  for(const std::string& blockName : blockPartNames) {
    stk::mesh::Part* part = bulk.mesh_meta_data().get_part(blockName);
    ThrowRequire(part != nullptr);
    ThrowRequire(part->primary_entity_rank() == rank);
    blockParts.push_back(part);
  }
  stk::mesh::Selector blockSelector =  stk::mesh::selectUnion(blockParts) & bulk.mesh_meta_data().locally_owned_part();
  unsigned localCount = stk::mesh::count_selected_entities(blockSelector, bulk.buckets(stk::topology::NODE_RANK));
  return stk::get_global_sum(MPI_COMM_WORLD, localCount);
}

void create_sideset(stk::mesh::BulkData& bulk,
                    const std::string& surfacePartName,
                    const std::vector<std::string>& blockPartNames)
{
  stk::mesh::ConstPartVector blockParts;
  for(const std::string& blockName : blockPartNames) {
    stk::mesh::Part* part = bulk.mesh_meta_data().get_part(blockName);
    ThrowRequire(part != nullptr);
    blockParts.push_back(part);
  }
  stk::mesh::Part& surfacePart = *bulk.mesh_meta_data().get_part(surfacePartName);

  bulk.mesh_meta_data().set_surface_to_block_mapping(&surfacePart, blockParts);
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
    if (bulk.is_valid(elem) && bulk.bucket(elem).owned()) {
      elems.push_back(elem);
    }
  }

  bulk.batch_change_entity_parts(elems, stk::mesh::PartVector{&toBlock}, stk::mesh::PartVector{&fromBlock});
}

unsigned get_num_common_entities(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks, const stk::mesh::EntityRank rank)
{
  stk::mesh::Selector selector;
  if(!blocks.empty()) {
    for (size_t i = 0; i < blocks.size()-1; ++i) {
      for (size_t j = i+1; j < blocks.size(); ++j) {
        selector |= (*blocks[i] & *blocks[j]);
      }
    }
  }

  selector &= bulk.mesh_meta_data().locally_owned_part();

  unsigned localNumCommonEntities = stk::mesh::count_selected_entities(selector, bulk.buckets(rank));

  return stk::get_global_sum(MPI_COMM_WORLD, localNumCommonEntities);
}

unsigned get_num_intersecting_nodes(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks)
{
  return get_num_common_entities(bulk, blocks, stk::topology::NODE_RANK);
}

unsigned get_num_common_sides(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks)
{
  return get_num_common_entities(bulk, blocks, bulk.mesh_meta_data().side_rank());
}

unsigned get_num_total_entities(const stk::mesh::BulkData & bulk, const stk::mesh::EntityRank rank)
{
  unsigned localNumTotalEntities = stk::mesh::count_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(rank));

  return stk::get_global_sum(MPI_COMM_WORLD, localNumTotalEntities);
}

unsigned get_num_total_nodes(const stk::mesh::BulkData & bulk)
{
  return get_num_total_entities(bulk, stk::topology::NODE_RANK);
}

unsigned get_num_total_sides(const stk::mesh::BulkData & bulk)
{
  return get_num_total_entities(bulk, bulk.mesh_meta_data().side_rank());
}

bool check_orphaned_nodes(stk::mesh::BulkData & bulk)
{
  bool foundOrphanedNode = false;
  for (stk::mesh::Bucket * bucket : bulk.buckets(stk::topology::NODE_RANK)) {
    for (stk::mesh::Entity node : *bucket) {
      const unsigned numElems = bulk.num_elements(node);
      if (numElems == 0u) {
        foundOrphanedNode = true;
        std::cout << "[p" << bulk.parallel_rank() << "] Found orphaned node: " << bulk.entity_key(node) << std::endl;
      }
    }
  }
  return foundOrphanedNode;
}

void output_mesh(stk::mesh::BulkData & bulk, const std::string & fileName)
{
  std::string writeOutput = stk::unit_test_util::get_option("--output", "off");
  if (writeOutput == "on") {
    stk::io::write_mesh(fileName, bulk);
  }
}

void output_mesh(stk::mesh::BulkData & bulk)
{
  const std::string fileName = std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".g";
  output_mesh(bulk, fileName);
}

int get_debug_level()
{
  int level = stk::unit_test_util::get_command_line_option("--debug", 0);
  return std::max(level, 0);
}

// create mesh
stk::mesh::Part & create_part(stk::mesh::MetaData& meta, const stk::topology topology, const std::string & blockName, int64_t blockId)
{
  stk::mesh::Part& part = meta.declare_part_with_topology(blockName, topology);
  stk::io::put_io_part_attribute(part);
  meta.set_part_id(part, blockId);
  return part;
}

stk::mesh::PartVector setup_mesh_1block_1quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, {&block1}));

  return {&block1};
}

stk::mesh::PartVector setup_mesh_2block_1quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1";
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(4u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_2block_2quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2";
  }
  else {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_2";
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(2u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(6u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_2block_2quad_reversed(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_1";
  }
  else {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_1";
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(2u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(6u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_2block_4quad_corner(stk::mesh::BulkData& bulk, int decompPattern = 0)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else {
    if (decompPattern == 1) {
      // p0 for block_1 and p1 for block_2
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 2) {
      // p0 for bottom half and p1 for top half
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 3) {
      // p0 for non-face-adjacent block_1 element, p1 for block_2 and face-adjacent block_1 elements
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 4) {
      // p0 diagonal, p1 off-diagonal (checkerboard)
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(3u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_2block_4quad_swappedCorner(stk::mesh::BulkData& bulk, int decompPattern = 0)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (bulk.parallel_size() == 2) {
    if (decompPattern == 1) {
      // p0 for block_1 and p1 for block_2
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 2) {
      // p0 for bottom half and p1 for top half
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 3) {
      // p0 for non-face-adjacent block_2 element, p1 for block_1 and face-adjacent block_2 elements
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 4) {
      // p0 diagonal, p1 off-diagonal (checkerboard)
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 2) {
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 3) {
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 4) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "2,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(3u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector create_3_blocks_order1(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order2(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order3(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order4(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order5(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  return blocks;
}

stk::mesh::PartVector create_3_blocks_order6(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector blocks(3);
  blocks[0] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  blocks[1] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  blocks[2] = &create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  return blocks;
}

void setup_mesh_3block_4quad_base(stk::mesh::BulkData& bulk, stk::mesh::PartVector & blocks, unsigned decompPattern)
{
  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_3";
  }
  else if (bulk.parallel_size() == 2) {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else if (decompPattern == 2) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else if (decompPattern == 3) {
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "2,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else if (decompPattern == 2) {
      meshDesc = "2,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else if (decompPattern == 3) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_3";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(4u, get_num_intersecting_nodes(bulk, blocks));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));
}

stk::mesh::PartVector setup_mesh_3block_4quad(stk::mesh::BulkData& bulk, unsigned blockOrder, unsigned decompPattern)
{
  stk::mesh::PartVector blocks;
  if (blockOrder == 1) {
    blocks = create_3_blocks_order1(bulk);
  } else if (blockOrder == 2) {
    blocks = create_3_blocks_order2(bulk);
  } else if (blockOrder == 3) {
    blocks = create_3_blocks_order3(bulk);
  } else if (blockOrder == 4) {
    blocks = create_3_blocks_order4(bulk);
  } else if (blockOrder == 5) {
    blocks = create_3_blocks_order5(bulk);
  } else if (blockOrder == 6) {
    blocks = create_3_blocks_order6(bulk);
  } else {
    std::cerr << "ERROR: Unexpected part ordinal ordering!!!" << std::endl;
    exit(1);
  }

  setup_mesh_3block_4quad_base(bulk, blocks, decompPattern);
  return blocks;
}

void test_mesh_3block_4quad(stk::mesh::BulkData& bulk, unsigned blockOrder, unsigned decompPattern) {
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(bulk, blockOrder, decompPattern);

  output_mesh(bulk, "disconnect_3block_4quad_blockOrder" + std::to_string(blockOrder) + "_decomp" + std::to_string(decompPattern) + "_init.g");

  stk::tools::disconnect_all_blocks(bulk);

  EXPECT_EQ(0u,  get_num_intersecting_nodes(bulk, blocks));
  EXPECT_EQ(14u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  output_mesh(bulk, "disconnect_3block_4quad_blockOrder" + std::to_string(blockOrder) + "_decomp" + std::to_string(decompPattern) + ".g");
}

stk::mesh::PartVector setup_mesh_3block_4quad_keepLowerRight(stk::mesh::BulkData& bulk, unsigned decompPattern)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (bulk.parallel_size() == 2) {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 2) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 3) {
      meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "2,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 2) {
      meshDesc = "2,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else if (decompPattern == 3) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "2,4,QUAD_4_2D,5,6,9,8,block_2";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(4u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2, &block3};
}

stk::mesh::PartVector setup_mesh_2block_4quad_checkerboard(stk::mesh::BulkData& bulk, unsigned decompPattern)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 2) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 3) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(5u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_3block_4quad_checkerboard(stk::mesh::BulkData& bulk, unsigned decompPattern)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else {
    if (decompPattern == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 2) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else if (decompPattern == 3) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
                 "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "1,4,QUAD_4_2D,5,6,9,8,block_1";
    }
    else {
      std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
      exit(1);
    }
  }
  std::vector<double> coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(5u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2, &block3};
}

stk::mesh::PartVector setup_mesh_2block_2quad_diagonal(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);

  std::string meshDesc;
  if (bulk.parallel_size() == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
               "0,2,QUAD_4_2D,4,5,7,6,block_2";
  }
  else {
    meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
               "1,2,QUAD_4_2D,4,5,7,6,block_2";
  }
  std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 1,2, 2,2 };

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(1u, get_num_intersecting_nodes(bulk, {&block1, &block2}));
  EXPECT_EQ(7u, get_num_total_nodes(bulk));

  return {&block1, &block2};
}

stk::mesh::PartVector setup_mesh_3block_4quad_bowtie(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);

  std::string meshDesc = "0,1,QUAD_4_2D,1,2, 6, 5,block_1\n"
                         "0,2,QUAD_4_2D,3,4, 7, 6,block_2\n"
                         "0,3,QUAD_4_2D,8,6,11,10,block_3\n"
                         "0,4,QUAD_4_2D,6,9,13,12,block_1";
  std::vector<double> coordinates = { 0,0, 0.9,0, 1.1,0, 2,0, 0,0.9, 1,1, 2,0.9, 0,1.1, 2,1.1, 0,2, 0.9,2, 1.1,2, 2,2 };
  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ( 1u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3}));
  EXPECT_EQ(13u, get_num_total_nodes(bulk));

  return {&block1, &block2, &block3};
}

void fill_mesh_description_4block_4quad_np1(stk::mesh::BulkData& bulk, unsigned blockOrder,
    std::string& meshDesc, std::vector<double>& coordinates)
{
  ThrowRequire(bulk.parallel_size() == 1);

  if (blockOrder == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else if (blockOrder == 2) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (blockOrder == 3) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_4";
  }
  coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
}

void fill_mesh_description_4block_4quad_np2(stk::mesh::BulkData& bulk, unsigned blockOrder,
    std::string& meshDesc, std::vector<double>& coordinates)
{
  ThrowRequire(bulk.parallel_size() == 2);
  ThrowRequire(blockOrder <= 3 && blockOrder > 0);

  if (blockOrder == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "1,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else if (blockOrder == 2) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "1,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (blockOrder == 3) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
               "1,4,QUAD_4_2D,5,6,9,8,block_4";
  }
  coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
}

void fill_mesh_description_4block_4quad_np3(stk::mesh::BulkData& bulk, unsigned blockOrder,
    std::string& meshDesc, std::vector<double>& coordinates)
{
  ThrowRequire(bulk.parallel_size() == 3);
  ThrowRequire(blockOrder <= 3 && blockOrder > 0);

  if (blockOrder == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "2,4,QUAD_4_2D,5,6,9,8,block_1";
  }
  else if (blockOrder == 2) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
               "2,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (blockOrder == 3) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
               "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
               "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
               "2,4,QUAD_4_2D,5,6,9,8,block_4";
  }
  coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
}

void fill_mesh_description_4block_4quad_np4(stk::mesh::BulkData& bulk, unsigned blockOrder,
    std::string& meshDesc, std::vector<double>& coordinates)
{
  ThrowRequire(bulk.parallel_size() == 4);
  ThrowRequire(blockOrder <= 4 && blockOrder > 0);

  if (blockOrder == 1) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
        "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
        "2,3,QUAD_4_2D,4,5,8,7,block_3\n"
        "3,4,QUAD_4_2D,5,6,9,8,block_4";
  }
  else if (blockOrder == 2) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
        "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
        "2,3,QUAD_4_2D,4,5,8,7,block_2\n"
        "3,4,QUAD_4_2D,5,6,9,8,block_3";
  }
  else if (blockOrder == 3) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
        "1,2,QUAD_4_2D,2,3,6,5,block_4\n"
        "2,3,QUAD_4_2D,4,5,8,7,block_1\n"
        "3,4,QUAD_4_2D,5,6,9,8,block_2";
  }
  else if (blockOrder == 4) {
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
        "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
        "2,3,QUAD_4_2D,4,5,8,7,block_4\n"
        "3,4,QUAD_4_2D,5,6,9,8,block_1";
  }

  coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
}

stk::mesh::PartVector setup_mesh_4block_4quad(stk::mesh::BulkData& bulk, unsigned blockOrder)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);

  std::string meshDesc;
  std::vector<double> coordinates;
  switch(bulk.parallel_size()) {
  case 1:
    fill_mesh_description_4block_4quad_np1(bulk, blockOrder, meshDesc, coordinates);
    break;
  case 2:
    fill_mesh_description_4block_4quad_np2(bulk, blockOrder, meshDesc, coordinates);
    break;
  case 3:
    fill_mesh_description_4block_4quad_np3(bulk, blockOrder, meshDesc, coordinates);
    break;
  case 4:
    fill_mesh_description_4block_4quad_np4(bulk, blockOrder, meshDesc, coordinates);
    break;
  default:
    ThrowRequireMsg(false, "Unexpected proc count for this test\n");
    break;
  }

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(5u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3, &block4}));
  EXPECT_EQ(9u, get_num_total_nodes(bulk));

  return {&block1, &block2, &block3, &block4};
}

void fill_mesh_description_6block_6quad_np1(stk::mesh::BulkData& bulk, std::string& meshDesc, std::vector<double>& coordinates)
{
  ThrowRequire(bulk.parallel_size() == 1);

  meshDesc = "0,1,QUAD_4_2D,1,2,6,5,block_1\n"
             "0,2,QUAD_4_2D,2,3,7,6,block_2\n"
             "0,3,QUAD_4_2D,3,4,8,7,block_3\n"
             "0,4,QUAD_4_2D,5,6,10,9,block_4\n"
             "0,5,QUAD_4_2D,6,7,11,10,block_5\n"
             "0,6,QUAD_4_2D,7,8,12,11,block_6\n";

  coordinates = { 0,0, 1,0, 2,0, 3,0,
                  0,1, 1,1, 2,1, 3,1,
                  0,2, 1,2, 2,2, 3,2};
}

stk::mesh::PartVector setup_mesh_6block_6quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);
  stk::mesh::Part & block5 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_5", 5);
  stk::mesh::Part & block6 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_6", 6);

  std::string meshDesc;
  std::vector<double> coordinates;
  switch(bulk.parallel_size()) {
  case 1:
    fill_mesh_description_6block_6quad_np1(bulk, meshDesc, coordinates);
    break;
  default:
    ThrowRequireMsg(false, "Unexpected proc count for this test\n");
    break;
  }

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(8u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3, &block4, &block5, &block6}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));

  bulk.dump_mesh_per_proc("dump");
  output_mesh(bulk, "initial.g");

  return {&block1, &block2, &block3, &block4, &block5, &block6};
}

void fill_mesh_description_9block_9quad_np1(stk::mesh::BulkData& bulk, std::string& meshDesc, std::vector<double>& coordinates)
{
  ThrowRequire(bulk.parallel_size() == 1);

  meshDesc = "0,1,QUAD_4_2D,1,2,6,5,block_1\n"
             "0,2,QUAD_4_2D,2,3,7,6,block_2\n"
             "0,3,QUAD_4_2D,3,4,8,7,block_3\n"
             "0,4,QUAD_4_2D,5,6,10,9,block_4\n"
             "0,5,QUAD_4_2D,6,7,11,10,block_5\n"
             "0,6,QUAD_4_2D,7,8,12,11,block_6\n"
             "0,7,QUAD_4_2D,9,10,14,13,block_7\n"
             "0,8,QUAD_4_2D,10,11,15,14,block_8\n"
             "0,9,QUAD_4_2D,11,12,16,15,block_9\n";

  coordinates = { 0,0, 1,0, 2,0, 3,0,
                  0,1, 1,1, 2,1, 3,1,
                  0,2, 1,2, 2,2, 3,2,
                  0,3, 1,3, 2,3, 3,3};
}

stk::mesh::PartVector setup_mesh_9block_9quad(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part & block1 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_1", 1);
  stk::mesh::Part & block2 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_2", 2);
  stk::mesh::Part & block3 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_3", 3);
  stk::mesh::Part & block4 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_4", 4);
  stk::mesh::Part & block5 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_5", 5);
  stk::mesh::Part & block6 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_6", 6);
  stk::mesh::Part & block7 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_7", 7);
  stk::mesh::Part & block8 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_8", 8);
  stk::mesh::Part & block9 = create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4_2D, "block_9", 9);

  std::string meshDesc;
  std::vector<double> coordinates;
  switch(bulk.parallel_size()) {
  case 1:
    fill_mesh_description_9block_9quad_np1(bulk, meshDesc, coordinates);
    break;
  default:
    ThrowRequireMsg(false, "Unexpected proc count for this test\n");
    break;
  }

  stk::unit_test_util::setup_text_mesh(bulk, meshDesc, coordinates);

  EXPECT_EQ(12u, get_num_intersecting_nodes(bulk, {&block1, &block2, &block3, &block4, &block5, &block6, &block7, &block8, &block9}));
  EXPECT_EQ(16u, get_num_total_nodes(bulk));

  bulk.dump_mesh_per_proc("dump");
  output_mesh(bulk, "initial.g");

  return {&block1, &block2, &block3, &block4, &block5, &block6, &block7, &block8, &block9};
}

stk::mesh::PartVector setup_mesh_1block_1hex(stk::mesh::BulkData& bulk)
{
  stk::io::fill_mesh("generated:1x1x1", bulk);

  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, {block1}));
  EXPECT_EQ(8u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1};
}

stk::mesh::PartVector setup_mesh_2block_1hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::io::fill_mesh("generated:1x1x1", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ(8u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::io::fill_mesh("generated:1x1x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2hex_withInternalSides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);

  stk::io::fill_mesh("generated:1x1x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

  create_sideset(bulk, "surface_1", {"block_1"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ( 1u, get_num_common_sides(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_EQ( 1u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2hex_withExternalSides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);

  stk::io::fill_mesh("generated:1x1x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

  create_sideset(bulk, "surface_1", {"block_1", "block_2"});
  create_all_boundary_sides(bulk, "surface_1");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_EQ(10u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_3block_4hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);

  stk::io::fill_mesh("generated:1x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_3");

  EXPECT_EQ( 8u, get_num_intersecting_nodes(bulk, {block1, block2, block3}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3};
}

stk::mesh::PartVector setup_mesh_3block_4hex_withInternalSides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_2", 2);

  stk::io::fill_mesh("generated:1x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_3");

  create_sideset(bulk, "surface_1", {"block_1"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");
  create_sides_between_blocks(bulk, "block_1", "block_3", "surface_1");
  create_sideset(bulk, "surface_2", {"block_2"});
  create_sides_between_blocks(bulk, "block_2", "block_3", "surface_2");

  EXPECT_EQ( 8u, get_num_intersecting_nodes(bulk, {block1, block2, block3}));
  EXPECT_EQ( 3u, get_num_common_sides(bulk, {block1, block2, block3}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_EQ( 3u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3};
}

stk::mesh::PartVector setup_mesh_2block_2cubeOfTet(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_2", 2);

  stk::io::fill_mesh("generated:1x1x2|tets", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12},
      "block_1", "block_2");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_2block_2cubeOfTet_withInternalSides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_2", 2);
  create_part(bulk.mesh_meta_data(), stk::topology::TRI_3, "surface_1", 1);

  stk::io::fill_mesh("generated:1x1x2|tets", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12},
      "block_1", "block_2");

  create_sideset(bulk, "surface_1", {"block_1"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");

  EXPECT_EQ( 4u, get_num_intersecting_nodes(bulk, {block1, block2}));
  EXPECT_EQ( 2u, get_num_common_sides(bulk, {block1, block2}));
  EXPECT_EQ(12u, get_num_total_nodes(bulk));
  EXPECT_EQ( 2u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2};
}

stk::mesh::PartVector setup_mesh_3block_4cubeOfTet(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_3", 3);

  stk::io::fill_mesh("generated:1x2x2|tets", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{19, 20, 21, 22, 23, 24}, "block_1", "block_3");

  EXPECT_EQ( 8u, get_num_intersecting_nodes(bulk, {block1, block2, block3}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3};
}

stk::mesh::PartVector setup_mesh_3block_4cubeOfTet_withInternalSides(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::TET_4, "block_3", 3);
  create_part(bulk.mesh_meta_data(), stk::topology::TRI_3, "surface_1", 1);
  create_part(bulk.mesh_meta_data(), stk::topology::TRI_3, "surface_2", 2);

  stk::io::fill_mesh("generated:1x2x2|tets", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{19, 20, 21, 22, 23, 24}, "block_1", "block_3");

  create_sideset(bulk, "surface_1", {"block_1"});
  create_sides_between_blocks(bulk, "block_1", "block_2", "surface_1");
  create_sides_between_blocks(bulk, "block_1", "block_3", "surface_1");
  create_sideset(bulk, "surface_2", {"block_2"});
  create_sides_between_blocks(bulk, "block_2", "block_3", "surface_2");

  EXPECT_EQ( 8u, get_num_intersecting_nodes(bulk, {block1, block2, block3}));
  EXPECT_EQ( 6u, get_num_common_sides(bulk, {block1, block2, block3}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_EQ( 6u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3};
}

stk::mesh::PartVector setup_mesh_4block_4hex(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part * block4 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);

  stk::io::fill_mesh("generated:1x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{3}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_4");

  EXPECT_EQ(10u, get_num_intersecting_nodes(bulk, {block1, block2, block3, block4}));
  EXPECT_EQ(18u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3, block4};
}

stk::mesh::PartVector setup_mesh_4block_4hex_vertical_stack(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part * block4 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);

  stk::io::fill_mesh("generated:1x1x4", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{3}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_4");

  EXPECT_EQ(12u, get_num_intersecting_nodes(bulk, {block1, block2, block3, block4}));
  EXPECT_EQ(20u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3, block4};
}

stk::mesh::PartVector setup_mesh_4block_8hex_cube(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part * block4 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);

  stk::io::fill_mesh("generated:2x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{3}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{5}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{6}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7}, "block_1", "block_4");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{8}, "block_1", "block_4");

  EXPECT_EQ(15u, get_num_intersecting_nodes(bulk, {block1, block2, block3, block4}));
  EXPECT_EQ(27u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3, block4};
}

stk::mesh::PartVector setup_mesh_8block_8hex_cube(stk::mesh::BulkData& bulk)
{
  stk::mesh::Part * block2 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_2", 2);
  stk::mesh::Part * block3 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_3", 3);
  stk::mesh::Part * block4 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_4", 4);
  stk::mesh::Part * block5 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_5", 5);
  stk::mesh::Part * block6 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_6", 6);
  stk::mesh::Part * block7 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_7", 7);
  stk::mesh::Part * block8 = &create_part(bulk.mesh_meta_data(), stk::topology::HEX_8, "block_8", 8);

  stk::io::fill_mesh("generated:2x2x2", bulk);
  stk::mesh::Part * block1 = bulk.mesh_meta_data().get_part("block_1");

  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{3}, "block_1", "block_3");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{4}, "block_1", "block_4");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{5}, "block_1", "block_5");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{6}, "block_1", "block_6");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{7}, "block_1", "block_7");
  move_elems_from_block_to_block(bulk, std::vector<stk::mesh::EntityId>{8}, "block_1", "block_8");

  EXPECT_EQ(19u, get_num_intersecting_nodes(bulk, {block1, block2, block3, block4, block5, block6, block7, block8}));
  EXPECT_EQ(27u, get_num_total_nodes(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return {block1, block2, block3, block4, block5, block6, block7, block8};
}

std::vector<std::string> get_part_names(const stk::mesh::PartVector& parts)
{
  std::vector<std::string> partNames;

  for(stk::mesh::Part* part : parts) {
    partNames.push_back(part->name());
  }

  return partNames;
}

stk::mesh::PartVector setup_mesh_8block_8hex_withExternalSides(stk::mesh::BulkData& bulk)
{
  create_part(bulk.mesh_meta_data(), stk::topology::QUAD_4, "surface_1", 1);

  stk::mesh::PartVector createdParts = setup_mesh_8block_8hex_cube(bulk);

  create_sideset(bulk, "surface_1", get_part_names(createdParts));
  create_all_boundary_sides(bulk, "surface_1");

  EXPECT_EQ(27u, get_num_total_nodes(bulk));
  EXPECT_EQ(24u, get_num_total_sides(bulk));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  return createdParts;
}
}//namespace

class TestDisconnectBlocks2D : public stk::unit_test_util::MeshFixture
{
protected:
  TestDisconnectBlocks2D() : stk::unit_test_util::MeshFixture(2)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }

};

TEST_F(TestDisconnectBlocks2D, disconnect_1block_1quad)
{
  stk::mesh::PartVector blocks = setup_mesh_1block_1quad(get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(4u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_1block_1quad.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_1quad)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_1quad(get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(4u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_1quad.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad(get_bulk());

  output_mesh(get_bulk(), "disconnect_2block_2quad_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_2quad.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad_updateGraph)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad(get_bulk());

  get_bulk().initialize_face_adjacent_element_graph();
  stk::tools::disconnect_all_blocks(get_bulk());

  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
  if (get_bulk().is_valid(elem1)) {
    EXPECT_EQ(0u, get_bulk().get_face_adjacent_element_graph().get_num_connected_elems(elem1));
  }
  if (get_bulk().is_valid(elem2)) {
    EXPECT_EQ(0u, get_bulk().get_face_adjacent_element_graph().get_num_connected_elems(elem2));
  }

  create_all_sides(get_bulk(), get_meta().locally_owned_part(), {&get_meta().get_topology_root_part(stk::topology::LINE_2)}, false);

  unsigned localNumSides = stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::EDGE_RANK));
  unsigned numSides = stk::get_global_sum(get_bulk().parallel(), localNumSides);
  EXPECT_EQ(8u, numSides);
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad_reversed)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_reversed(get_bulk());

  output_mesh(get_bulk(), "disconnect_2block_2quad_reversed_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_2quad_reversed.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 4);

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp4_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_corner_decomp4.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 4);

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp4_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_swappedCorner_decomp4.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder1_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 1, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder1_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 1, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder1_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 1, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder2_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 2, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder2_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 2, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder2_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 2, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder3_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 3, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder3_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 3, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder3_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 3, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder4_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 4, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder4_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 4, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder4_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 4, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder5_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 5, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder5_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 5, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder5_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 5, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder6_decomp1)
{
  test_mesh_3block_4quad(get_bulk(), 6, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder6_decomp2)
{
  test_mesh_3block_4quad(get_bulk(), 6, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder6_decomp3)
{
  test_mesh_3block_4quad(get_bulk(), 6, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_keepLowerRight_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_4quad_checkerboard_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 1);

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp1_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp1.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 2);

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp2_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp2.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 3);

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp3_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_3block_4quad_checkerboard_decomp3.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad_diagonal)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_diagonal(get_bulk());

  output_mesh(get_bulk(), "disconnect_2block_2quad_diagonal_init.g");

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh(get_bulk(), "disconnect_2block_2quad_diagonal.g");
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_bowtie)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4quad_bowtie(get_bulk());

    output_mesh(get_bulk(), "disconnect_3block_4quad_bowtie_init.g");

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ( 0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4quad_bowtie.g");
  }
}

typedef TestDisconnectBlocks2D TestReconnectList2D;

void test_reconnect_list(const stk::mesh::BulkData& bulk, const BlockConnectionVector& disconnectList, const BlockConnectionVector& expectedReconnectList)
{
  stk::tools::impl::PartPairLess comparator;
  stk::tools::BlockPairVector reconnectList = get_local_reconnect_list(bulk, disconnectList);
  EXPECT_EQ(reconnectList.size(), expectedReconnectList.size());

  stk::tools::BlockPairVector convertedExpectedList = convert_connection_vector_to_pair_vector(bulk, expectedReconnectList);

  for(stk::tools::BlockPair& reconnectPair : reconnectList) {
    EXPECT_TRUE(std::binary_search(convertedExpectedList.begin(), convertedExpectedList.end(), reconnectPair, comparator))
        << " did not find: " << reconnectPair.first->name() << " and " << reconnectPair.second->name() << std::endl;
  }
}

TEST_F(TestReconnectList2D, reconnect_2block_2quad)
{
  stk::mesh::BulkData& bulk = get_bulk();
  setup_mesh_2block_2quad(bulk);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2") };
  BlockConnectionVector expectedReconnectList;

  test_reconnect_list(bulk, blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_3block_4quad_permutation1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  setup_mesh_3block_4quad(bulk, 1, 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_3") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_2"), BlockConnection("block_2", "block_3") };

  test_reconnect_list(bulk, blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_3block_4quad_permutation2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  setup_mesh_3block_4quad(bulk, 1, 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_3"), BlockConnection("block_2", "block_3") };

  test_reconnect_list(bulk, blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_3block_4quad_permutation3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  setup_mesh_3block_4quad(bulk, 1, 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_3") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_2"), BlockConnection("block_1", "block_3") };

  test_reconnect_list(bulk, blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_4block_4quad_permutation1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) { return; }
  setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_3"), BlockConnection("block_2", "block_3"), BlockConnection("block_2", "block_4"), BlockConnection("block_1", "block_4"), BlockConnection("block_3", "block_4") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_4block_4quad_permutation2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) { return; }
  setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_1","block_3") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_2", "block_3"), BlockConnection("block_2", "block_4"), BlockConnection("block_1", "block_4"), BlockConnection("block_3", "block_4") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_4block_4quad_permutation3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) { return; }
  setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_3", "block_4") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_2", "block_3"), BlockConnection("block_2", "block_4"), BlockConnection("block_1", "block_4"), BlockConnection("block_1", "block_3") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_4block_4quad_permutation4)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 4) { return; }
  setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_2", "block_4"), BlockConnection("block_3", "block_4"), BlockConnection("block_1", "block_3") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_2", "block_3"),  BlockConnection("block_1", "block_4") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_9block_9quad_permutation1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_1", "block_4") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_5"), BlockConnection("block_2", "block_4"), BlockConnection("block_2", "block_5"), BlockConnection("block_5", "block_4") };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}

TEST_F(TestReconnectList2D, reconnect_9block_9quad_permutation2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_5","block_2"), BlockConnection("block_5", "block_4"), BlockConnection("block_5", "block_6"), BlockConnection("block_5", "block_8") };
  BlockConnectionVector expectedReconnectList{ BlockConnection("block_1", "block_5"), BlockConnection("block_1", "block_2"), BlockConnection("block_1", "block_4"), BlockConnection("block_2", "block_4"),
                                               BlockConnection("block_3", "block_5"), BlockConnection("block_2", "block_6"), BlockConnection("block_2", "block_3"), BlockConnection("block_3", "block_6"),
                                               BlockConnection("block_5", "block_9"), BlockConnection("block_6", "block_8"), BlockConnection("block_8", "block_9"), BlockConnection("block_6", "block_9"),
                                               BlockConnection("block_5", "block_7"), BlockConnection("block_4", "block_8"), BlockConnection("block_4", "block_7"), BlockConnection("block_7", "block_8"),
                                             };

  test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
}


class TestDisconnectBlocks : public stk::unit_test_util::MeshFixture
{
protected:
  TestDisconnectBlocks() : stk::unit_test_util::MeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }
};

typedef TestDisconnectBlocks TestReconnectList;

TEST_F(TestReconnectList, reconnect_8block_8hex_permutation1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnect { BlockConnection("block_1","block_2"), BlockConnection("block_2","block_4"), BlockConnection("block_2","block_6"),
                                                   BlockConnection("block_5","block_7"), BlockConnection("block_5","block_6"), BlockConnection("block_5","block_1"),
                                                   BlockConnection("block_8","block_6"), BlockConnection("block_8","block_4"), BlockConnection("block_8","block_7"),
                                                   BlockConnection("block_3","block_7"), BlockConnection("block_3","block_4"), BlockConnection("block_3","block_1")};
    BlockConnectionVector expectedReconnectList { BlockConnection("block_2","block_3"), BlockConnection("block_2","block_5"), BlockConnection("block_2","block_8"), BlockConnection("block_2","block_7"),
                                                  BlockConnection("block_3","block_5"), BlockConnection("block_3","block_8"), BlockConnection("block_5","block_8"), BlockConnection("block_1","block_8"),
                                                  BlockConnection("block_1","block_4"), BlockConnection("block_1","block_6"), BlockConnection("block_1","block_7"), BlockConnection("block_3","block_6"),
                                                  BlockConnection("block_4","block_6"), BlockConnection("block_4","block_7"), BlockConnection("block_6","block_7"), BlockConnection("block_4","block_5")
                                                };
    test_reconnect_list(get_bulk(), blockPairsToDisconnect, expectedReconnectList);
  }
}

TEST_F(TestDisconnectBlocks, disconnect_1block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_1block_1hex(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_1block_1hex.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_2block_1hex(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_1hex.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2hex.g");
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_2block_2hex_withInternalSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex_withInternalSides(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2hex_withInternalSides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2hex_withExternalSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex_withExternalSides(get_bulk());

    EXPECT_EQ(12u, get_num_surface_nodes(get_bulk(), {"surface_1"}));

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(16u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2hex_withExternalSides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_8block_8hex_withExternalSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_withExternalSides(get_bulk());

    EXPECT_EQ(26u, get_num_surface_nodes(get_bulk(), {"surface_1"}));

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(56u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(64u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_8block_8hex_withExternalSides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2cubeOfTet.g");
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_2block_2cubeOfTet_withInternalSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet_withInternalSides(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2cubeOfTet_withInternalSides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_3block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex(get_bulk());

    output_mesh(get_bulk(), "disconnect_3block_4hex_init.g");

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4hex.g");
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_3block_4hex_withInternalSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex_withInternalSides(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4hex_withInternalSides.g");
  }
}

TEST_F(TestDisconnectBlocks, disconnect_3block_4cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4cubeOfTet.g");
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_3block_4cubeOfTet_withInternalSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet_withInternalSides(get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_3block_4cubeOfTet_withInternalSides.g");
  }
}

TEST(DisconnectBlocks, input_mesh)
{
  double startTime = stk::wall_time();

  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  std::string exodusFileName = stk::unit_test_util::get_option("-f", "");
  if (exodusFileName.empty()) return;

  stk::io::fill_mesh_with_auto_decomp(exodusFileName, bulk);

  double meshReadTime = stk::wall_time();

  stk::tools::disconnect_all_blocks(bulk);

  double disconnectTime = stk::wall_time();

  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  EXPECT_EQ(0u,  get_num_intersecting_nodes(bulk, allBlocksInMesh));
  EXPECT_FALSE(check_orphaned_nodes(bulk));

  output_mesh(bulk, "disconnected_" + exodusFileName);
  double meshWriteTime = stk::wall_time();

  if (bulk.parallel_rank() == 0) {
    std::cout << " Mesh read time = " << (meshReadTime - startTime) << " s" << std::endl;
    std::cout << "Disconnect time = " << (disconnectTime - meshReadTime) << " s" << std::endl;
    std::cout << "Mesh write time = " << (meshWriteTime - disconnectTime) << " s" << std::endl;
  }
}

typedef TestDisconnectBlocks2D TestReconnectBlocks2D;
typedef TestDisconnectBlocks   TestReconnectBlocks;
typedef TestDisconnectBlocks2D TestDisconnectFullAlgorithm2D;
typedef TestDisconnectBlocks2D TestDisconnectFullAlgorithmPartial2D;
typedef TestDisconnectBlocks   TestDisconnectFullAlgorithm;
typedef TestDisconnectBlocks   TestDisconnectFullAlgorithmPartial;

void test_one_block_reconnect(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  ThrowRequire(allBlocksInMesh.size() == 1u);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;
  info.debugLevel = get_debug_level();

  std::vector<stk::tools::BlockPair> blockPairs;

  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::reconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  output_mesh(bulk);
}

void test_two_block_reconnect(stk::mesh::BulkData& bulk, const stk::mesh::PartVector& blocks,
                              const unsigned expectedCommonNodeCount)
{
  ThrowRequire(blocks.size() >= 2u);
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;
  info.debugLevel = get_debug_level();

  std::vector<stk::tools::BlockPair> blockPairs;

  stk::tools::impl::insert_block_pair(blocks[0], blocks[1], blockPairs);

  EXPECT_EQ(expectedCommonNodeCount, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::reconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(expectedCommonNodeCount, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  output_mesh(bulk);
}

TEST_F(TestReconnectBlocks2D, reconnect_1block_1quad)
{
  stk::mesh::PartVector blocks = setup_mesh_1block_1quad(get_bulk());
  test_one_block_reconnect(get_bulk());
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_2quad)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad(get_bulk());
  test_two_block_reconnect(get_bulk(), blocks, 2u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_2quad_diagonal)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_diagonal(get_bulk());
  test_two_block_reconnect(get_bulk(), blocks, 1u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_2quad_reversed)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_reversed(get_bulk());
  test_two_block_reconnect(get_bulk(), blocks, 2u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 1);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 2);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 3);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(get_bulk(), 4);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 1);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 2);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 3);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(get_bulk(), 4);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 1);
  test_two_block_reconnect(get_bulk(), blocks, 5u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 2);
  test_two_block_reconnect(get_bulk(), blocks, 5u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(get_bulk(), 3);
  test_two_block_reconnect(get_bulk(), blocks, 5u);
}

void test_reconnect_block_pairs(stk::mesh::BulkData& bulk, stk::mesh::PartVector& blocks,
                                const std::vector<BlockConnection>& reconnectPairs,
                                unsigned expectedGlobalInitialCommonNodes,
                                unsigned expectedGlobalConnectCommonNodes)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = false;
  info.debugLevel = get_debug_level();

  std::vector<stk::tools::BlockPair> disconnectBlockPairs;
  for(unsigned i = 0; i < allBlocksInMesh.size() - 1; i++) {
    for(unsigned j = i+1; j < allBlocksInMesh.size(); j++) {
      stk::tools::impl::insert_block_pair(allBlocksInMesh[i], allBlocksInMesh[j], disconnectBlockPairs);
    }
  }

  std::vector<stk::tools::BlockPair> reconnectBlockPairs;

  for(const BlockConnection& connectPair : reconnectPairs) {
    ThrowRequire(connectPair.block1 != connectPair.block2);

    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    ThrowRequire(block1 != nullptr && block2 != nullptr);

    stk::tools::impl::insert_block_pair(block1, block2, reconnectBlockPairs);
  }

  EXPECT_EQ(expectedGlobalInitialCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, disconnectBlockPairs, info);
  bulk.dump_mesh_per_proc("dump");
  output_mesh(bulk, "test.g");
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::reconnect_block_pairs(bulk, reconnectBlockPairs, info);
  EXPECT_EQ(expectedGlobalConnectCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  stk::mesh::PartVector connectedBlockPair(2);
  for(const BlockConnection& connectPair : reconnectPairs) {
    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    connectedBlockPair[0] = block1;
    connectedBlockPair[1] = block2;

    EXPECT_EQ(connectPair.numExpectedIntersectingNodes, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " first block: " << block1->name() << " second block: " << block2->name();
  }
  output_mesh(bulk);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp1_v2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2), BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2), BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp2_v2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp2_v3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 2, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2), BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder3_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 3, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_3","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder3_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 3, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_3","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder3_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 3, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_3","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder4_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 4, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_2",2), BlockConnection("block_3","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder4_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 4, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_2",2), BlockConnection("block_3","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder4_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 4, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_2",2), BlockConnection("block_3","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder5_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 5, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder5_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 5, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder5_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 5, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder6_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 6, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_3",2), BlockConnection("block_1","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder6_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 6, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_3",2), BlockConnection("block_1","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder6_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 6, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_3",2), BlockConnection("block_1","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_keepLowerRight_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_keepLowerRight_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_keepLowerRight_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(get_bulk(), 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",3), BlockConnection("block_1","block_3",3)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 5);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",3), BlockConnection("block_1","block_3",3)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 5);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(get_bulk(), 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",3), BlockConnection("block_1","block_3",3)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 5);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_bowtie)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_bowtie(get_bulk());
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",1), BlockConnection("block_1","block_3",1)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 1, 1);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_4",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1_v2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_4",2), BlockConnection("block_2","block_3",1)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1_v3)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_4",2), BlockConnection("block_3","block_4",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1_v4)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_2","block_3",1), BlockConnection("block_3","block_4",2)/*, BlockConnection("block_1","block_4",1)*/};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1_v5)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_1","block_2",2), BlockConnection("block_2","block_4",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 4);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 2);
  output_mesh(get_bulk(), "initial.g");
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",1)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 2);
}

TEST_F(TestReconnectBlocks, reconnect_1block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_1block_1hex(get_bulk());
    test_one_block_reconnect(get_bulk());
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_2block_1hex(get_bulk());
    test_two_block_reconnect(get_bulk(), blocks, 0u);
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_2hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex(get_bulk());
    test_two_block_reconnect(get_bulk(), blocks, 4u);
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_2cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet(get_bulk());
    test_two_block_reconnect(get_bulk(), blocks, 4u);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",4), BlockConnection("block_2","block_1",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 6);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4hex_v2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4hex_v3)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet(get_bulk());
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",4), BlockConnection("block_2","block_1",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 6);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_1",4), BlockConnection("block_3","block_2",2) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 10, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex_v2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_3",2), BlockConnection("block_1","block_4",2) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 10, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex_v3)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_1","block_2",4), BlockConnection("block_2","block_4",4), BlockConnection("block_4","block_3",4),
                                          BlockConnection("block_2","block_3",2), BlockConnection("block_1","block_4",2)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 10, 8);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex_vertical_stack)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 4) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex_vertical_stack(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_1","block_2",4), BlockConnection("block_3","block_4",4) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 12, 8);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_4hex_vertical_stack_v2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 4) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4hex_vertical_stack(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_3",4) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 12, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_8hex_cube)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_4",6) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 15, 6);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_8hex_cube_v2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_3",3) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 15, 3);
  }
}

TEST_F(TestReconnectBlocks, reconnect_4block_8hex_cube_v3)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_2","block_3",3), BlockConnection("block_3","block_4",6) };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 15, 6);
  }
}

TEST_F(TestReconnectBlocks, DISABLED_reconnect_8block_8hex_cube)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector reconnectPairs{ BlockConnection("block_1","block_3",4),
                                          BlockConnection("block_3","block_4",4),
                                          BlockConnection("block_4","block_8",4),
                                          BlockConnection("block_8","block_7",4),
                                          BlockConnection("block_7","block_5",4),
                                          BlockConnection("block_5","block_6",4)
                                        };
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 19, 16);
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_2hex_withExternalSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex_withExternalSides(get_bulk());

    EXPECT_EQ(12u, get_num_surface_nodes(get_bulk(), {"surface_1"}));

    stk::tools::impl::LinkInfo info;
    info.preserveOrphans = true;
    info.debugLevel = get_debug_level();

    stk::tools::disconnect_all_blocks(get_bulk(), info);

    EXPECT_EQ(16u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_2block_2hex_withExternalSides.g");

    stk::tools::impl::reconnect_block_pairs(get_bulk(), {stk::tools::BlockPair(blocks[0], blocks[1])}, info);

    EXPECT_EQ(12u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(4u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(12u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk());
  }
}

TEST_F(TestReconnectBlocks, reconnect_8block_8hex_withExternalSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_withExternalSides(get_bulk());

    EXPECT_EQ(26u, get_num_surface_nodes(get_bulk(), {"surface_1"}));

    stk::tools::impl::LinkInfo info;
    info.preserveOrphans = true;
    info.debugLevel = get_debug_level();

    stk::tools::disconnect_all_blocks(get_bulk(), info);

    EXPECT_EQ(56u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(64u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk(), "disconnect_8block_8hex_withExternalSides.g");

    stk::tools::BlockPairVector blockPairs { stk::tools::impl::get_block_pair(blocks[0], blocks[1]),
                                             stk::tools::impl::get_block_pair(blocks[0], blocks[2]),
                                             stk::tools::impl::get_block_pair(blocks[0], blocks[3]),
                                             stk::tools::impl::get_block_pair(blocks[1], blocks[2]),
                                             stk::tools::impl::get_block_pair(blocks[1], blocks[3]),
                                             stk::tools::impl::get_block_pair(blocks[2], blocks[3])
                                           };

    stk::tools::impl::reconnect_block_pairs(get_bulk(), blockPairs, info);

    EXPECT_EQ(45u, get_num_surface_nodes(get_bulk(), {"surface_1"}));
    EXPECT_EQ(10u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(50u,get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh(get_bulk());
  }
}

stk::tools::BlockPairVector convert_connection_vector_to_pair_vector(const stk::mesh::BulkData& bulk, BlockConnectionVector& disconnectConnVector)
{
  stk::tools::BlockPairVector pairVector;
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();

  for(BlockConnection& connection : disconnectConnVector) {
    stk::mesh::Part* part1 = meta.get_part(connection.block1);
    stk::mesh::Part* part2 = meta.get_part(connection.block2);

    pairVector.push_back(stk::tools::impl::get_block_pair(part1, part2));
  }

  return pairVector;
}

void test_populate_blocks_to_reconnect(stk::mesh::BulkData& bulk, const stk::tools::BlockPairVector& expectedReconnectBlockPairs, BlockConnectionVector& disconnectPairs)
{
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::BlockPairVector orderedBlockPairsInMesh;
  stk::tools::BlockPairVector blockPairsToReconnect;
  stk::tools::BlockPairVector blockPairsToDisconnect = convert_connection_vector_to_pair_vector(bulk, disconnectPairs);

  stk::tools::impl::fill_ordered_block_pairs(allBlocksInMesh, orderedBlockPairsInMesh);

  stk::tools::impl::populate_blocks_to_reconnect(bulk, orderedBlockPairsInMesh, blockPairsToDisconnect, blockPairsToReconnect);
  EXPECT_EQ(expectedReconnectBlockPairs.size(), blockPairsToReconnect.size());

  for(const stk::tools::BlockPair& blockPair : expectedReconnectBlockPairs) {
    EXPECT_TRUE(std::binary_search(blockPairsToReconnect.begin(), blockPairsToReconnect.end(), blockPair, stk::tools::impl::PartPairLess()))
        << "first block name: " << blockPair.first->name() << " second block name: " << blockPair.second->name() << std::endl;
  }
}

TEST_F(TestDisconnectBlocks2D, test_populate_blocks_to_reconnect_2block_2quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad(bulk);
  BlockConnectionVector disconnectPairs{BlockConnection("block_1","block_2",2)};
  stk::tools::BlockPairVector expectedReconnectBlockPairs;

  test_populate_blocks_to_reconnect(bulk, expectedReconnectBlockPairs, disconnectPairs);
}

TEST_F(TestDisconnectBlocks2D, test_populate_blocks_to_reconnect_3block_4quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(bulk, 1, 1);

  stk::mesh::Part* block1 = meta.get_part("block_1");
  stk::mesh::Part* block2 = meta.get_part("block_2");
  stk::mesh::Part* block3 = meta.get_part("block_3");

  stk::tools::BlockPairVector expectedReconnectBlockPairs;
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block1, block3));
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block2, block3));

  BlockConnectionVector disconnectPairs{BlockConnection("block_1","block_2",2)};

  test_populate_blocks_to_reconnect(bulk, expectedReconnectBlockPairs, disconnectPairs);
}

TEST_F(TestDisconnectBlocks2D, test_populate_blocks_to_reconnect_3block_4quad_v2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(bulk, 1, 1);
  BlockConnectionVector disconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1", "block_3", 2)};

  stk::mesh::Part* block2 = meta.get_part("block_2");
  stk::mesh::Part* block3 = meta.get_part("block_3");

  stk::tools::BlockPairVector expectedReconnectBlockPairs;
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block2, block3));

  test_populate_blocks_to_reconnect(bulk, expectedReconnectBlockPairs, disconnectPairs);
}

TEST_F(TestDisconnectBlocks2D, test_populate_blocks_to_reconnect_4block_4quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }
  stk::mesh::BulkData& bulk = get_bulk();
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(bulk, 1);

  stk::mesh::Part* block1 = meta.get_part("block_1");
  stk::mesh::Part* block2 = meta.get_part("block_2");
  stk::mesh::Part* block3 = meta.get_part("block_3");
  stk::mesh::Part* block4 = meta.get_part("block_4");

  stk::tools::BlockPairVector expectedReconnectBlockPairs;
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block1, block3));
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block1, block4));
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block2, block3));
  expectedReconnectBlockPairs.push_back(stk::tools::impl::get_block_pair(block3, block4));

  BlockConnectionVector disconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_4","block_2",2)};

  test_populate_blocks_to_reconnect(bulk, expectedReconnectBlockPairs, disconnectPairs);
}

void print_block_pair_list(const stk::tools::BlockPairVector& blockPairs)
{
  return;
  std::cout << "Block pairs to reconnect" << std::endl;
  for(const stk::tools::BlockPair& blockPair : blockPairs) {
    std::cout << "\t block pair: {" << blockPair.first->name() << "," << blockPair.second->name() <<"}" << std::endl;
  }
  std::cout << std::endl;
}

void test_block_disconnect_full_algorithm_partial_disconnect(stk::mesh::BulkData& bulk,
                                                             BlockConnectionVector& blockPairConnectionsToDisconnect,
                                                             unsigned expectedGlobalInitialCommonNodes,
                                                             unsigned expectedGlobalConnectCommonNodesAfterHingeSnip,
                                                             unsigned expectedGlobalTotalNumberOfNodes)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;
  info.debugLevel = get_debug_level();

  stk::tools::BlockPairVector blockPairsToDisconnect = convert_connection_vector_to_pair_vector(bulk, blockPairConnectionsToDisconnect);
  stk::tools::BlockPairVector blockPairsToReconnect;
  blockPairsToReconnect = get_local_reconnect_list(bulk, blockPairConnectionsToDisconnect);
  print_block_pair_list(blockPairsToReconnect);

  info.print_debug_msg_p0(1) << "PRINTING RECONNECT PAIRS" << std::endl;
  for(stk::tools::BlockPair& blockPair : blockPairsToReconnect) {
    info.print_debug_msg_p0(1,false) << "\t{" << blockPair.first->name() << ", " << blockPair.second->name() << "}" << std::endl;
  }

  EXPECT_EQ(expectedGlobalInitialCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, blockPairsToDisconnect, info);
  output_mesh(bulk, "dis_seq_disconnected_mesh.g");

  stk::mesh::PartVector connectedBlockPair(2);
  for(const stk::tools::BlockPair& connectPair :  blockPairsToDisconnect) {
    connectedBlockPair[0] = connectPair.first;
    connectedBlockPair[1] = connectPair.second;

    EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " block1: " << connectedBlockPair[0]->name() << " block2: " << connectedBlockPair[1]->name();
  }

  stk::tools::impl::reconnect_block_pairs(bulk, blockPairsToReconnect, info);

  output_mesh(bulk, "dis_seq_reconnected_mesh.g");

  stk::tools::impl::snip_all_hinges_between_blocks(bulk, info.debugLevel > 0);

  EXPECT_EQ(expectedGlobalConnectCommonNodesAfterHingeSnip, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  EXPECT_EQ(expectedGlobalTotalNumberOfNodes, get_num_total_nodes(bulk));

  for(const BlockConnection& connectPair :  blockPairConnectionsToDisconnect) {
    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    connectedBlockPair[0] = block1;
    connectedBlockPair[1] = block2;

    EXPECT_EQ(connectPair.numExpectedIntersectingNodes, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " block1: " << block1->name() << " block2: " << block2->name();
  }
  output_mesh(bulk, "dis_seq_hinge_snipped_mesh.g");
}

void test_block_disconnect_full_algorithm(stk::mesh::BulkData& bulk,
                                          BlockConnectionVector& blockPairConnectionsToDisconnect,
                                          unsigned expectedGlobalInitialCommonNodes,
                                          unsigned expectedGlobalConnectCommonNodesAfterHingeSnip,
                                          unsigned expectedGlobalTotalNumberOfNodes)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;
  info.debugLevel = get_debug_level();

  stk::tools::BlockPairVector orderedBlockPairsInMesh;
  stk::tools::impl::fill_ordered_block_pairs(allBlocksInMesh, orderedBlockPairsInMesh);

  stk::tools::BlockPairVector blockPairsToDisconnect = convert_connection_vector_to_pair_vector(bulk, blockPairConnectionsToDisconnect);
  stk::tools::BlockPairVector blockPairsToReconnect;
  stk::tools::impl::populate_blocks_to_reconnect(bulk, orderedBlockPairsInMesh, blockPairsToDisconnect, blockPairsToReconnect);
  print_block_pair_list(blockPairsToReconnect);

  info.print_debug_msg_p0(1) << "PRINTING RECONNECT PAIRS" << std::endl;
  for(stk::tools::BlockPair& blockPair : blockPairsToReconnect) {
    info.print_debug_msg_p0(1,false) << "\t{" << blockPair.first->name() << ", " << blockPair.second->name() << "}" << std::endl;
  }

  EXPECT_EQ(expectedGlobalInitialCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::disconnect_block_pairs(bulk, orderedBlockPairsInMesh, info);
  output_mesh(bulk, "dis_seq_disconnected_mesh.g");
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::impl::reconnect_block_pairs(bulk, blockPairsToReconnect, info);
//  EXPECT_EQ(expectedGlobalConnectCommonNodesAfterReconnect, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  output_mesh(bulk, "dis_seq_reconnected_mesh.g");

  stk::tools::impl::snip_all_hinges_between_blocks(bulk, info.debugLevel > 0);

  EXPECT_EQ(expectedGlobalConnectCommonNodesAfterHingeSnip, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  EXPECT_EQ(expectedGlobalTotalNumberOfNodes, get_num_total_nodes(bulk));

  stk::mesh::PartVector connectedBlockPair(2);
  for(const BlockConnection& connectPair :  blockPairConnectionsToDisconnect) {
    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    connectedBlockPair[0] = block1;
    connectedBlockPair[1] = block2;

    EXPECT_EQ(connectPair.numExpectedIntersectingNodes, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " block1: " << block1->name() << " block2: " << block2->name();
  }
  output_mesh(bulk, "dis_seq_hinge_snipped_mesh.g");
}

void test_connection_pairs(stk::mesh::BulkData& bulk, stk::mesh::PartVector& allBlocksInMesh,
                           BlockConnectionVector& blockPairConnectionsToDisconnect,
                           unsigned expectedFinalCommonNodeCount)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector connectedBlockPair(2);
  for(const BlockConnection& connectPair :  blockPairConnectionsToDisconnect) {
    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    connectedBlockPair[0] = block1;
    connectedBlockPair[1] = block2;

    EXPECT_EQ(connectPair.numExpectedIntersectingNodes, get_num_intersecting_nodes(bulk, connectedBlockPair)) << " block1: " << block1->name() << " block2: " << block2->name();
  }

  output_mesh(bulk);

  EXPECT_EQ(expectedFinalCommonNodeCount, get_num_intersecting_nodes(bulk, allBlocksInMesh));
}

void test_user_block_disconnect(stk::mesh::BulkData& bulk,
                               BlockConnectionVector& blockPairConnectionsToDisconnect,
                               unsigned expectedFinalCommonNodeCount)
{
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::BlockPairVector blockPairsToDisconnect = convert_connection_vector_to_pair_vector(bulk, blockPairConnectionsToDisconnect);

  stk::tools::disconnect_user_blocks(bulk, blockPairsToDisconnect);

  test_connection_pairs(bulk, allBlocksInMesh, blockPairConnectionsToDisconnect, expectedFinalCommonNodeCount);
}

stk::tools::BlockNamePairVector convert_connection_vector_to_name_pair_vector(const stk::mesh::BulkData& bulk, const BlockConnectionVector& blockNamePairConnections)
{
  stk::tools::BlockNamePairVector namePairVector;

  for(const BlockConnection& connection : blockNamePairConnections) {
    namePairVector.push_back(std::make_pair(connection.block1, connection.block2));
  }

  return namePairVector;
}

void test_named_user_block_disconnect(stk::mesh::BulkData& bulk,
                                      BlockConnectionVector& blockPairConnectionsToDisconnect,
                                      unsigned expectedFinalCommonNodeCount)
{
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::BlockNamePairVector blockPairsToDisconnect = convert_connection_vector_to_name_pair_vector(bulk, blockPairConnectionsToDisconnect);

  stk::tools::disconnect_user_blocks(bulk, blockPairsToDisconnect);

  test_connection_pairs(bulk, allBlocksInMesh, blockPairConnectionsToDisconnect, expectedFinalCommonNodeCount);
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_2block_2quad)
{
  stk::mesh::PartVector initialMesh = setup_mesh_2block_2quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{BlockConnection("block_1","block_2",0)};
  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnect, 2u, 0u, 8u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_2block_2quad)
{
  stk::mesh::PartVector initialMesh = setup_mesh_2block_2quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{BlockConnection("block_1","block_2",0)};
  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnect, 2u, 0u, 8u);
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_3block_4quad_blockOrder1_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_1","block_3",0)};
  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 4u, 2u, 12u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_3block_4quad_blockOrder1_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_1","block_3",0)};
  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 4u, 2u, 12u);
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_3block_4quad_blockOrder1_decomp1_pacman)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_3",1)};
  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 4u, 3u, 10u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_3block_4quad_blockOrder1_decomp1_pacman)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(get_bulk(), 1, 1);
  BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_3",1)};
  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 4u, 3u, 10u);
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_4block_4quad_blockOrder1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) <= 3) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_4",0), BlockConnection("block_3","block_4",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 5u, 3u, 12u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_4block_4quad_blockOrder1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) <= 3) {
    stk::mesh::PartVector blocks = setup_mesh_4block_4quad(get_bulk(), 1);
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_4",0), BlockConnection("block_3","block_4",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 5u, 3u, 12u);
  }
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_6block_6quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_6block_6quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_5",1), BlockConnection("block_3","block_6",0),
                                                BlockConnection("block_2","block_6",0), BlockConnection("block_3","block_5",0) };

  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnect, 8u, 8u, 14u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_6block_6quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_6block_6quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_5",1), BlockConnection("block_3","block_6",0),
                                                BlockConnection("block_2","block_6",0), BlockConnection("block_3","block_5",0) };

  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnect, 8u, 8u, 14u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_9block_9quad)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_1","block_2"), BlockConnection("block_1", "block_4") };

  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnect, 12u, 10u, 19u);
}

TEST_F(TestDisconnectFullAlgorithmPartial2D, disconnect_full_algorithm_9block_9quad_v2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_3"), BlockConnection("block_5","block_6"),
                                                BlockConnection("block_3","block_5"), BlockConnection("block_2","block_6"), BlockConnection("block_9", "block_6")};

  test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnect, 12u, 11u, 20u);
}

TEST_F(TestDisconnectFullAlgorithm2D, disconnect_full_algorithm_9block_9quad_v2)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  setup_mesh_9block_9quad(get_bulk());
//  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_3"), BlockConnection("block_5", "block_6"), BlockConnection("block_9", "block_6")};
  BlockConnectionVector blockPairsToDisconnect{ BlockConnection("block_2","block_3"), BlockConnection("block_5","block_6"),
                                                BlockConnection("block_3","block_5"), BlockConnection("block_2","block_6"), BlockConnection("block_9", "block_6")};

  test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnect, 12u, 11u, 20u);
}

TEST_F(TestDisconnectFullAlgorithm, disconnect_full_algorithm_4block_8hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_4",0), BlockConnection("block_1","block_2",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 15u, 9u, 36u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial, disconnect_full_algorithm_4block_8hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_4block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_2","block_4",0), BlockConnection("block_1","block_2",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 15u, 9u, 36u);
  }
}

TEST_F(TestDisconnectFullAlgorithm, disconnect_full_algorithm_8block_8hex_lower_corner)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 19u, 16u, 34u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial, disconnect_full_algorithm_8block_8hex_lower_corner)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 19u, 16u, 34u);
  }
}

TEST_F(TestDisconnectFullAlgorithm, disconnect_full_algorithm_8block_8hex_opposite_corners)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_7","block_8",0), BlockConnection("block_5","block_7",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 19u, 13u, 41u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial, disconnect_full_algorithm_8block_8hex_opposite_corners)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_7","block_8",0), BlockConnection("block_5","block_7",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 19u, 13u, 41u);
  }
}

TEST_F(TestDisconnectFullAlgorithm, disconnect_full_algorithm_8block_8hex_opposite_corners2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_5","block_7",0), BlockConnection("block_5","block_6",0), BlockConnection("block_5","block_1",0),
                                                      BlockConnection("block_8","block_6",0), BlockConnection("block_8","block_4",0), BlockConnection("block_8","block_7",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_3","block_4",0), BlockConnection("block_3","block_1",0)};
    test_block_disconnect_full_algorithm(get_bulk(), blockPairsToDisconnectPairs, 19u, 0u, 64u);
  }
}

TEST_F(TestDisconnectFullAlgorithmPartial, disconnect_full_algorithm_8block_8hex_opposite_corners2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());
    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_5","block_7",0), BlockConnection("block_5","block_6",0), BlockConnection("block_5","block_1",0),
                                                      BlockConnection("block_8","block_6",0), BlockConnection("block_8","block_4",0), BlockConnection("block_8","block_7",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_3","block_4",0), BlockConnection("block_3","block_1",0)};
    test_block_disconnect_full_algorithm_partial_disconnect(get_bulk(), blockPairsToDisconnectPairs, 19u, 0u, 64u);
  }
}

std::string get_basename(const std::string& path)
{
  char tempPath[path.length()+1];
  tempPath[path.length()] = 0;
  path.copy(tempPath, path.length());
  return std::string(basename(tempPath));
}

TEST(TestDisconnectInputFile, input_mesh)
{
  double startTime = stk::wall_time();

  stk::mesh::MetaData meta(3);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);

  std::string exodusFileName = stk::unit_test_util::get_option("-exoFile", "");
  if (exodusFileName.empty()) return;

  std::string disconnectBlockFile = stk::unit_test_util::get_option("-blockFile", "");
  if (disconnectBlockFile.empty()) return;

  stk::io::fill_mesh_with_auto_decomp(exodusFileName, bulk);

  std::ifstream infile(disconnectBlockFile);

  std::string block1;
  std::string block2;
  stk::tools::BlockPairVector disconnectBlockVec;

  while(infile >> block1 >> block2) {
    stk::mesh::Part* block1Part = meta.get_part(block1);
    stk::mesh::Part* block2Part = meta.get_part(block2);

    ThrowRequire(block1Part != nullptr);
    ThrowRequire(block2Part != nullptr);
    ThrowRequire(block1Part != block2Part);

    disconnectBlockVec.push_back(stk::tools::impl::get_block_pair(block1Part, block2Part));
  }
  infile.close();

  double meshReadTime = stk::wall_time();

  std::cout << "Starting disconnect block sequence" << std::endl;
  stk::tools::disconnect_user_blocks_partial(bulk, disconnectBlockVec, get_debug_level());

  double disconnectTime = stk::wall_time();

  std::string filename = get_basename(exodusFileName);
  output_mesh(bulk, "disconnected_" + filename);
  double meshWriteTime = stk::wall_time();

  if (bulk.parallel_rank() == 0) {
    std::cout << "Mesh read time = " << (meshReadTime - startTime) << " s" << std::endl;
    std::cout << "Mesh write time = " << (meshWriteTime - disconnectTime) << " s" << std::endl;
  }
}

typedef TestDisconnectFullAlgorithm TestDisconnectUserBlocks;

TEST_F(TestDisconnectUserBlocks, disconnect_user_blocks_8block_8hex_opposite_corners)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());

    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_7","block_8",0), BlockConnection("block_5","block_7",0)};
    test_user_block_disconnect(get_bulk(), blockPairsToDisconnectPairs, 13u);
  }
}

TEST_F(TestDisconnectUserBlocks, disconnect_named_user_blocks_8block_8hex_opposite_corners)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_8block_8hex_cube(get_bulk());

    BlockConnectionVector blockPairsToDisconnectPairs{BlockConnection("block_1","block_2",0), BlockConnection("block_2","block_4",0), BlockConnection("block_2","block_6",0),
                                                      BlockConnection("block_3","block_7",0), BlockConnection("block_7","block_8",0), BlockConnection("block_5","block_7",0)};
    test_named_user_block_disconnect(get_bulk(), blockPairsToDisconnectPairs, 13u);
  }
}

typedef TestDisconnectBlocks2D TestBlockPairCreation;

TEST_F(TestBlockPairCreation, test_block_pair_creation)
{
  setup_mesh_3block_4quad(get_bulk(), 1, 1);
  stk::mesh::PartVector parts;
  stk::tools::BlockPairVector partPairs;
  parts.push_back(get_meta().get_part("block_1"));
  parts.push_back(get_meta().get_part("block_2"));
  parts.push_back(get_meta().get_part("block_3"));
  stk::tools::impl::fill_ordered_block_pairs(parts, partPairs);

  EXPECT_EQ(3u, partPairs.size());
  EXPECT_EQ("block_1", partPairs[0].first->name());
  EXPECT_EQ("block_2", partPairs[0].second->name());
  EXPECT_EQ("block_1", partPairs[1].first->name());
  EXPECT_EQ("block_3", partPairs[1].second->name());
  EXPECT_EQ("block_2", partPairs[2].first->name());
  EXPECT_EQ("block_3", partPairs[2].second->name());
}
