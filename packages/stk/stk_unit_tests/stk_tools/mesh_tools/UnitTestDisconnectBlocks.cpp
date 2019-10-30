#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocks.hpp>
#include <stk_tools/mesh_tools/DisconnectBlocksImpl.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/environment/WallTime.hpp>
#include "stk_unit_test_utils/getOption.h"
#include "stk_util/parallel/ParallelReduce.hpp"
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

  bulk.mesh_meta_data().set_surface_to_block_mapping(&surfacePart, stk::mesh::ConstPartVector{&blockPart});
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

//  for (stk::mesh::Bucket * bucket : bulk.get_buckets(stk::topology::NODE_RANK, selector)) {
//    for (stk::mesh::Entity node : *bucket) {
//      std::cerr << "[p" << bulk.parallel_rank() << "] Common node: " << bulk.entity_key(node) << std::endl;
//    }
//  }

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

void output_mesh(const std::string & fileName, stk::mesh::BulkData & bulk)
{
  std::string writeOutput = stk::unit_test_util::get_option("--output", "off");
  if (writeOutput == "on") {
    stk::io::write_mesh(fileName, bulk);
  }
}

void output_mesh(stk::mesh::BulkData & bulk)
{
  const std::string fileName = std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()) + ".g";
  output_mesh(fileName, bulk);
}

int get_debug_level()
{
  int level = stk::unit_test_util::get_command_line_option("--debug", 0);
  return std::max(level, 0);
}

}//namespace



class TestDisconnectBlocks2D : public stk::unit_test_util::MeshFixture
{
protected:
  TestDisconnectBlocks2D() : stk::unit_test_util::MeshFixture(2)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
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

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), {&block1}));

    return {&block1};
  }

  stk::mesh::PartVector setup_mesh_2block_1quad()
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

    std::string meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1";
    std::vector<double> coordinates = { 0,0, 1,0, 0,1, 1,1 };
    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2}));
    EXPECT_EQ(4u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2};
  }

  stk::mesh::PartVector setup_mesh_2block_2quad()
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

    std::string meshDesc;
    std::vector<double> coordinates;
    if (get_bulk().parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2";
      coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
    }
    else {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_2";
      if (get_bulk().parallel_rank() == 0) {
        coordinates = { 0,0, 1,0, 0,1, 1,1 };
      }
      else {
        coordinates = { 1,0, 2,0, 1,1, 2,1 };
      }
    }
    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(2u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2}));
    EXPECT_EQ(6u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2};
  }

  stk::mesh::PartVector setup_mesh_2block_2quad_reversed()
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

    std::string meshDesc;
    std::vector<double> coordinates;
    if (get_bulk().parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1";
      coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
    }
    else {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "1,2,QUAD_4_2D,2,3,6,5,block_1";
      if (get_bulk().parallel_rank() == 0) {
        coordinates = { 0,0, 1,0, 0,1, 1,1 };
      }
      else if (get_bulk().parallel_rank() == 1) {
        coordinates = { 1,0, 2,0, 1,1, 2,1 };
      }
    }
    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(2u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2}));
    EXPECT_EQ(6u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2};
  }

  stk::mesh::PartVector setup_mesh_2block_4quad_corner(int decompPattern = 0)
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

    std::string meshDesc;
    std::vector<double> coordinates;
    if (get_bulk().parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_1";
      coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
    }
    else {
      if (decompPattern == 1) {
        // p0 for block_1 and p1 for block_2
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                   "0,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,0, 2,0, 1,1, 2,1 };
        }
      }
      else if (decompPattern == 2) {
        // p0 for bottom half and p1 for top half
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 3) {
        // p0 for non-face-adjacent block_1 element, p1 for block_2 and face-adjacent block_1 elements
        meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 4) {
        // p0 diagonal, p1 off-diagonal (checkerboard)
        meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else {
        std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
        exit(1);
      }
    }

    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(3u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2}));
    EXPECT_EQ(9u, get_num_total_nodes(get_bulk()));

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

    EXPECT_EQ(3u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2}));
    EXPECT_EQ(9u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2};
  }

  stk::mesh::PartVector setup_mesh_2block_4quad_swappedCorner(int decompPattern = 0)
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

    std::string meshDesc;
    std::vector<double> coordinates;
    if (get_bulk().parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
      coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
    }
    else if (get_bulk().parallel_size() == 2) {
      if (decompPattern == 1) {
        // p0 for block_1 and p1 for block_2
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "0,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,0, 2,0, 1,1, 2,1 };
        }
      }
      else if (decompPattern == 2) {
        // p0 for bottom half and p1 for top half
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 3) {
        // p0 for non-face-adjacent block_2 element, p1 for block_1 and face-adjacent block_2 elements
        meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 4) {
        // p0 diagonal, p1 off-diagonal (checkerboard)
        meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 1,2, 2,2 };
        }
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
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 2) {
        meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "2,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 1,0, 2,0, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 3) {
        meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_2\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "2,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,0, 1,0, 0,1, 1,1 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 4) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_2\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "2,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "0,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,0, 2,0, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 0,1, 1,1, 0,2, 1,2 };
        }
      }
    }

    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(3u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2}));
    EXPECT_EQ(9u, get_num_total_nodes(get_bulk()));

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

  void setup_mesh_3block_4quad_base(stk::mesh::PartVector & blocks, unsigned decompPattern)
  {
    std::string meshDesc;
    std::vector<double> coordinates;
    if (get_bulk().parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_3";
      coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
    }
    else if (get_bulk().parallel_size() == 2) {
      if (decompPattern == 1) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_3";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,0, 2,0, 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 2) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_3";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 3) {
        meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                   "0,4,QUAD_4_2D,5,6,9,8,block_3";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,0, 1,0, 0,1, 1,1 };
        }
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
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,1, 2,1, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 1,0, 2,0, 1,1, 2,1 };
        }
      }
      else if (decompPattern == 2) {
        meshDesc = "2,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                   "0,4,QUAD_4_2D,5,6,9,8,block_3";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 1,0, 2,0, 1,1, 2,1, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 0,0, 1,0, 0,1, 1,1 };
        }
      }
      else if (decompPattern == 3) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_1\n"
                   "2,4,QUAD_4_2D,5,6,9,8,block_3";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else {
        std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
        exit(1);
      }
    }

    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(4u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(9u, get_num_total_nodes(get_bulk()));
  }

  stk::mesh::PartVector setup_mesh_3block_4quad(unsigned blockOrder, unsigned decompPattern)
  {
    stk::mesh::PartVector blocks;
    if (blockOrder == 1) {
      blocks = create_3_blocks_order1();
    } else if (blockOrder == 2) {
      blocks = create_3_blocks_order2();
    } else if (blockOrder == 3) {
      blocks = create_3_blocks_order3();
    } else if (blockOrder == 4) {
      blocks = create_3_blocks_order4();
    } else if (blockOrder == 5) {
      blocks = create_3_blocks_order5();
    } else if (blockOrder == 6) {
      blocks = create_3_blocks_order6();
    } else {
      std::cerr << "ERROR: Unexpected part ordinal ordering!!!" << std::endl;
      exit(1);
    }

    setup_mesh_3block_4quad_base(blocks, decompPattern);
    return blocks;
  }

  void test_mesh_3block_4quad(unsigned blockOrder, unsigned decompPattern) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4quad(blockOrder, decompPattern);

    output_mesh("disconnect_3block_4quad_blockOrder" + std::to_string(blockOrder) + "_decomp" + std::to_string(decompPattern) + "_init.g", get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_3block_4quad_blockOrder" + std::to_string(blockOrder) + "_decomp" + std::to_string(decompPattern) + ".g", get_bulk());
  }

  stk::mesh::PartVector setup_mesh_3block_4quad_keepLowerRight(unsigned decompPattern)
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);
    stk::mesh::Part & block3 = create_part(stk::topology::QUAD_4_2D, "block_3", 3);

    std::string meshDesc;
    std::vector<double> coordinates;
    if (get_bulk().parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_2";
      coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
    }
    else if (get_bulk().parallel_size() == 2) {
      if (decompPattern == 1) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,0, 2,0, 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 2) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 3) {
        meshDesc = "1,1,QUAD_4_2D,1,2,5,4,block_3\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "0,3,QUAD_4_2D,4,5,8,7,block_3\n"
                   "0,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,0, 1,0, 0,1, 1,1 };
        }
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
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,1, 2,1, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 1,0, 2,0, 1,1, 2,1 };
        }
      }
      else if (decompPattern == 2) {
        meshDesc = "2,1,QUAD_4_2D,1,2,5,4,block_3\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
                   "0,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 1,0, 2,0, 1,1, 2,1, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 0,0, 1,0, 0,1, 1,1 };
        }
      }
      else if (decompPattern == 3) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_3\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_3\n"
                   "2,4,QUAD_4_2D,5,6,9,8,block_2";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 0,2, 1,2 };
        }
        else if (get_bulk().parallel_rank() == 2) {
          coordinates = { 1,1, 2,1, 1,2, 2,2 };
        }
      }
      else {
        std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
        exit(1);
      }
    }
    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(4u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2, &block3}));
    EXPECT_EQ(9u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2, &block3};
  }

  stk::mesh::PartVector setup_mesh_2block_4quad_checkerboard(unsigned decompPattern)
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

    std::string meshDesc;
    std::vector<double> coordinates;
    if (get_bulk().parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_1";
      coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
    }
    else {
      if (decompPattern == 1) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 2) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "0,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2 };
        }
      }
      else if (decompPattern == 3) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_2\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
      }
      else {
        std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
        exit(1);
      }
    }

    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(5u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2}));
    EXPECT_EQ(9u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2};
  }

  stk::mesh::PartVector setup_mesh_3block_4quad_checkerboard(unsigned decompPattern)
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);
    stk::mesh::Part & block3 = create_part(stk::topology::QUAD_4_2D, "block_3", 3);

    std::string meshDesc;
    std::vector<double> coordinates;
    if (get_bulk().parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                 "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
                 "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
                 "0,4,QUAD_4_2D,5,6,9,8,block_1";
      coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
    }
    else {
      if (decompPattern == 1) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
      }
      else if (decompPattern == 2) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "0,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 1,2, 2,2 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2 };
        }
      }
      else if (decompPattern == 3) {
        meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                   "1,2,QUAD_4_2D,2,3,6,5,block_3\n"
                   "1,3,QUAD_4_2D,4,5,8,7,block_2\n"
                   "1,4,QUAD_4_2D,5,6,9,8,block_1";
        if (get_bulk().parallel_rank() == 0) {
          coordinates = { 0,0, 1,0, 0,1, 1,1 };
        }
        else if (get_bulk().parallel_rank() == 1) {
          coordinates = { 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
        }
      }
      else {
        std::cerr << "ERROR: Unexpected decomposition pattern (" << decompPattern << ")!" << std::endl;
        exit(1);
      }
    }

    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(5u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2, &block3}));
    EXPECT_EQ(9u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2, &block3};
  }

  stk::mesh::PartVector setup_mesh_2block_2quad_diagonal()
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);

    std::string meshDesc;
    std::vector<double> coordinates;
    if (get_bulk().parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                 "0,2,QUAD_4_2D,4,5,7,6,block_2";
      coordinates = { 0,0, 1,0, 0,1, 1,1, 2,1, 1,2, 2,2 };
    }
    else {
      meshDesc = "0,1,QUAD_4_2D,1,2,4,3,block_1\n"
                 "1,2,QUAD_4_2D,4,5,7,6,block_2";
      if (get_bulk().parallel_rank() == 0) {
        coordinates = { 0,0, 1,0, 0,1, 1,1 };
      }
      else if (get_bulk().parallel_rank() == 1) {
        coordinates = { 1,1, 2,1, 1,2, 2,2 };
      }
    }

    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(1u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2}));
    EXPECT_EQ(7u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2};
  }

  stk::mesh::PartVector setup_mesh_3block_4quad_bowtie()
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);
    stk::mesh::Part & block3 = create_part(stk::topology::QUAD_4_2D, "block_3", 3);

    std::string meshDesc = "0,1,QUAD_4_2D,1,2, 6, 5,block_1\n"
                           "0,2,QUAD_4_2D,3,4, 7, 6,block_2\n"
                           "0,3,QUAD_4_2D,8,6,11,10,block_3\n"
                           "0,4,QUAD_4_2D,6,9,13,12,block_1";
    std::vector<double> coordinates = { 0,0, 0.9,0, 1.1,0, 2,0, 0,0.9, 1,1, 2,0.9, 0,1.1, 2,1.1, 0,2, 0.9,2, 1.1,2, 2,2 };
    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ( 1u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2, &block3}));
    EXPECT_EQ(13u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2, &block3};
  }

  void fill_mesh_description_4block_4quad_np1(stk::mesh::BulkData& bulk,
                                              std::string& meshDesc, std::vector<double>& coordinates)
  {
    ThrowRequire(bulk.parallel_size() == 1);
    meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_4\n"
               "0,2,QUAD_4_2D,2,3,6,5,block_3\n"
               "0,3,QUAD_4_2D,4,5,8,7,block_2\n"
               "0,4,QUAD_4_2D,5,6,9,8,block_1";
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

    if (get_bulk().parallel_rank() == 0) {
      coordinates = { 0,0, 1,0, 2,0, 0,1, 1,1, 2,1 };
    }
    else if (get_bulk().parallel_rank() == 1) {
      coordinates = { 0,1, 1,1, 2,1, 0,2, 1,2, 2,2 };
    }
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

    if (get_bulk().parallel_rank() == 0) {
      coordinates = { 0,0, 1,0, 0,1, 1,1 };
    }
    else if (get_bulk().parallel_rank() == 1) {
      coordinates = { 1,0, 2,0, 1,0, 1,1, 2,1, 0,2, 1,2 };
    }
    else if (get_bulk().parallel_rank() == 2) {
      coordinates = { 1,1, 2,1, 1,2, 2,2 };
    }
  }

  stk::mesh::PartVector setup_mesh_4block_4quad(unsigned blockOrder)
  {
    stk::mesh::Part & block1 = create_part(stk::topology::QUAD_4_2D, "block_1", 1);
    stk::mesh::Part & block2 = create_part(stk::topology::QUAD_4_2D, "block_2", 2);
    stk::mesh::Part & block3 = create_part(stk::topology::QUAD_4_2D, "block_3", 3);
    stk::mesh::Part & block4 = create_part(stk::topology::QUAD_4_2D, "block_4", 4);

    std::string meshDesc;
    std::vector<double> coordinates;
    switch(get_bulk().parallel_size()) {
      case 1:
        fill_mesh_description_4block_4quad_np1(get_bulk(), meshDesc, coordinates);
        break;
      case 2:
        fill_mesh_description_4block_4quad_np2(get_bulk(), blockOrder, meshDesc, coordinates);
        break;
      case 3:
        fill_mesh_description_4block_4quad_np3(get_bulk(), blockOrder, meshDesc, coordinates);
        break;
      default:
        ThrowRequireMsg(false, "Unexpected proc count for this test\n");
        break;
    }

    stk::unit_test_util::fill_mesh_using_text_mesh_with_coordinates(meshDesc, coordinates, get_bulk());

    EXPECT_EQ(5u, get_num_intersecting_nodes(get_bulk(), {&block1, &block2, &block3, &block4}));
    EXPECT_EQ(9u, get_num_total_nodes(get_bulk()));

    return {&block1, &block2, &block3, &block4};
  }

};

TEST_F(TestDisconnectBlocks2D, disconnect_1block_1quad)
{
  stk::mesh::PartVector blocks = setup_mesh_1block_1quad();

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(4u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_1block_1quad.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_1quad)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_1quad();

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(4u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_1quad.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad();

  output_mesh("disconnect_2block_2quad_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_2quad.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad_updateGraph)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad();

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
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_reversed();

  output_mesh("disconnect_2block_2quad_reversed_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_2quad_reversed.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(1);

  output_mesh("disconnect_2block_4quad_corner_decomp1_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_corner_decomp1.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(2);

  output_mesh("disconnect_2block_4quad_corner_decomp2_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_corner_decomp2.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(3);

  output_mesh("disconnect_2block_4quad_corner_decomp3_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_corner_decomp3.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_corner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(4);

  output_mesh("disconnect_2block_4quad_corner_decomp4_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_corner_decomp4.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(1);

  output_mesh("disconnect_2block_4quad_swappedCorner_decomp1_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_swappedCorner_decomp1.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(2);

  output_mesh("disconnect_2block_4quad_swappedCorner_decomp2_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_swappedCorner_decomp2.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(3);

  output_mesh("disconnect_2block_4quad_swappedCorner_decomp3_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_swappedCorner_decomp3.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_swappedCorner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(4);

  output_mesh("disconnect_2block_4quad_swappedCorner_decomp4_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_swappedCorner_decomp4.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder1_decomp1)
{
  test_mesh_3block_4quad(1, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder1_decomp2)
{
  test_mesh_3block_4quad(1, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder1_decomp3)
{
  test_mesh_3block_4quad(1, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder2_decomp1)
{
  test_mesh_3block_4quad(2, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder2_decomp2)
{
  test_mesh_3block_4quad(2, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder2_decomp3)
{
  test_mesh_3block_4quad(2, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder3_decomp1)
{
  test_mesh_3block_4quad(3, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder3_decomp2)
{
  test_mesh_3block_4quad(3, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder3_decomp3)
{
  test_mesh_3block_4quad(3, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder4_decomp1)
{
  test_mesh_3block_4quad(4, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder4_decomp2)
{
  test_mesh_3block_4quad(4, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder4_decomp3)
{
  test_mesh_3block_4quad(4, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder5_decomp1)
{
  test_mesh_3block_4quad(5, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder5_decomp2)
{
  test_mesh_3block_4quad(5, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder5_decomp3)
{
  test_mesh_3block_4quad(5, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder6_decomp1)
{
  test_mesh_3block_4quad(6, 1);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder6_decomp2)
{
  test_mesh_3block_4quad(6, 2);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_blockOrder6_decomp3)
{
  test_mesh_3block_4quad(6, 3);
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(1);

  output_mesh("disconnect_3block_4quad_keepLowerRight_decomp1_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_3block_4quad_keepLowerRight_decomp1.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(2);

  output_mesh("disconnect_3block_4quad_keepLowerRight_decomp2_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_3block_4quad_keepLowerRight_decomp2.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_keepLowerRight_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(3);

  output_mesh("disconnect_3block_4quad_keepLowerRight_decomp3_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_3block_4quad_keepLowerRight_decomp3.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(1);

  output_mesh("disconnect_2block_4quad_checkerboard_decomp1_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_checkerboard_decomp1.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(2);

  output_mesh("disconnect_2block_4quad_checkerboard_decomp2_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_checkerboard_decomp2.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(3);

  output_mesh("disconnect_2block_4quad_checkerboard_decomp3_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(14u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_4quad_checkerboard_decomp3.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(1);

  output_mesh("disconnect_3block_4quad_checkerboard_decomp1_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_3block_4quad_checkerboard_decomp1.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(2);

  output_mesh("disconnect_3block_4quad_checkerboard_decomp2_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_3block_4quad_checkerboard_decomp2.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(3);

  output_mesh("disconnect_3block_4quad_checkerboard_decomp3_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_3block_4quad_checkerboard_decomp3.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_2block_2quad_diagonal)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_diagonal();

  output_mesh("disconnect_2block_2quad_diagonal_init.g", get_bulk());

  stk::tools::disconnect_all_blocks(get_bulk());

  EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
  EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
  EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

  output_mesh("disconnect_2block_2quad_diagonal.g", get_bulk());
}

TEST_F(TestDisconnectBlocks2D, disconnect_3block_4quad_bowtie)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4quad_bowtie();

    output_mesh("disconnect_3block_4quad_bowtie_init.g", get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ( 0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(15u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_3block_4quad_bowtie.g", get_bulk());
  }
}


class TestDisconnectBlocks : public stk::unit_test_util::MeshFixture
{
protected:
  TestDisconnectBlocks() : stk::unit_test_util::MeshFixture(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  }

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
    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x1x1", get_bulk());

    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), {block1}));
    EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1};
  }

  stk::mesh::PartVector setup_mesh_2block_1hex()
  {
    stk::mesh::Part * block2 = &create_part(stk::topology::HEX_8, "block_2", 2);
    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x1x1", get_bulk());
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), {block1, block2}));
    EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1, block2};
  }

  stk::mesh::PartVector setup_mesh_2block_2hex()
  {
    stk::mesh::Part * block2 = &create_part(stk::topology::HEX_8, "block_2", 2);
    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x1x2", get_bulk());
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

    EXPECT_EQ( 4u, get_num_intersecting_nodes(get_bulk(), {block1, block2}));
    EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1, block2};
  }

  stk::mesh::PartVector setup_mesh_2block_2hex_withSides()
  {
    stk::mesh::Part * block2 = &create_part(stk::topology::HEX_8, "block_2", 2);
    create_part(stk::topology::QUAD_4, "surface_1", 1);

    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x1x2", get_bulk());
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");

    create_sideset(get_bulk(), "surface_1", "block_1");
    create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");

    EXPECT_EQ( 4u, get_num_intersecting_nodes(get_bulk(), {block1, block2}));
    EXPECT_EQ( 1u, get_num_common_sides(get_bulk(), {block1, block2}));
    EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
    EXPECT_EQ( 1u, get_num_total_sides(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1, block2};
  }

  stk::mesh::PartVector setup_mesh_3block_4hex()
  {
    stk::mesh::Part * block2 = &create_part(stk::topology::HEX_8, "block_2", 2);
    stk::mesh::Part * block3 = &create_part(stk::topology::HEX_8, "block_3", 3);

    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x2x2", get_bulk());
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{4}, "block_1", "block_3");

    EXPECT_EQ( 8u, get_num_intersecting_nodes(get_bulk(), {block1, block2, block3}));
    EXPECT_EQ(18u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1, block2, block3};
  }

  stk::mesh::PartVector setup_mesh_3block_4hex_withSides()
  {
    stk::mesh::Part * block2 = &create_part(stk::topology::HEX_8, "block_2", 2);
    stk::mesh::Part * block3 = &create_part(stk::topology::HEX_8, "block_3", 3);
    create_part(stk::topology::QUAD_4, "surface_1", 1);
    create_part(stk::topology::QUAD_4, "surface_2", 2);

    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x2x2", get_bulk());
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{2}, "block_1", "block_2");
    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{4}, "block_1", "block_3");

    create_sideset(get_bulk(), "surface_1", "block_1");
    create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");
    create_sides_between_blocks(get_bulk(), "block_1", "block_3", "surface_1");
    create_sideset(get_bulk(), "surface_2", "block_2");
    create_sides_between_blocks(get_bulk(), "block_2", "block_3", "surface_2");

    EXPECT_EQ( 8u, get_num_intersecting_nodes(get_bulk(), {block1, block2, block3}));
    EXPECT_EQ( 3u, get_num_common_sides(get_bulk(), {block1, block2, block3}));
    EXPECT_EQ(18u, get_num_total_nodes(get_bulk()));
    EXPECT_EQ( 3u, get_num_total_sides(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1, block2, block3};
  }

  stk::mesh::PartVector setup_mesh_2block_2cubeOfTet()
  {
    stk::mesh::Part * block2 = &create_part(stk::topology::TET_4, "block_2", 2);

    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x1x2|tets", get_bulk());
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12},
                                   "block_1", "block_2");

    EXPECT_EQ( 4u, get_num_intersecting_nodes(get_bulk(), {block1, block2}));
    EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1, block2};
  }

  stk::mesh::PartVector setup_mesh_2block_2cubeOfTet_withSides()
  {
    stk::mesh::Part * block2 = &create_part(stk::topology::TET_4, "block_2", 2);
    create_part(stk::topology::TRI_3, "surface_1", 1);

    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x1x2|tets", get_bulk());
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12},
                                   "block_1", "block_2");

    create_sideset(get_bulk(), "surface_1", "block_1");
    create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");

    EXPECT_EQ( 4u, get_num_intersecting_nodes(get_bulk(), {block1, block2}));
    EXPECT_EQ( 2u, get_num_common_sides(get_bulk(), {block1, block2}));
    EXPECT_EQ(12u, get_num_total_nodes(get_bulk()));
    EXPECT_EQ( 2u, get_num_total_sides(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1, block2};
  }

  stk::mesh::PartVector setup_mesh_3block_4cubeOfTet()
  {
    stk::mesh::Part * block2 = &create_part(stk::topology::TET_4, "block_2", 2);
    stk::mesh::Part * block3 = &create_part(stk::topology::TET_4, "block_3", 3);

    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x2x2|tets", get_bulk());
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12}, "block_1", "block_2");
    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{19, 20, 21, 22, 23, 24}, "block_1", "block_3");

    EXPECT_EQ( 8u, get_num_intersecting_nodes(get_bulk(), {block1, block2, block3}));
    EXPECT_EQ(18u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1, block2, block3};
  }

  stk::mesh::PartVector setup_mesh_3block_4cubeOfTet_withSides()
  {
    stk::mesh::Part * block2 = &create_part(stk::topology::TET_4, "block_2", 2);
    stk::mesh::Part * block3 = &create_part(stk::topology::TET_4, "block_3", 3);
    create_part(stk::topology::TRI_3, "surface_1", 1);
    create_part(stk::topology::TRI_3, "surface_2", 2);

    //allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::fill_mesh("generated:1x2x2|tets", get_bulk());
    stk::mesh::Part * block1 = get_meta().get_part("block_1");

    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{7, 8, 9, 10, 11, 12}, "block_1", "block_2");
    move_elems_from_block_to_block(get_bulk(), std::vector<stk::mesh::EntityId>{19, 20, 21, 22, 23, 24}, "block_1", "block_3");

    create_sideset(get_bulk(), "surface_1", "block_1");
    create_sides_between_blocks(get_bulk(), "block_1", "block_2", "surface_1");
    create_sides_between_blocks(get_bulk(), "block_1", "block_3", "surface_1");
    create_sideset(get_bulk(), "surface_2", "block_2");
    create_sides_between_blocks(get_bulk(), "block_2", "block_3", "surface_2");

    EXPECT_EQ( 8u, get_num_intersecting_nodes(get_bulk(), {block1, block2, block3}));
    EXPECT_EQ( 6u, get_num_common_sides(get_bulk(), {block1, block2, block3}));
    EXPECT_EQ(18u, get_num_total_nodes(get_bulk()));
    EXPECT_EQ( 6u, get_num_total_sides(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    return {block1, block2, block3};
  }
};

TEST_F(TestDisconnectBlocks, disconnect_1block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_1block_1hex();

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_1block_1hex.g", get_bulk());
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_2block_1hex();

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u, get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(8u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_2block_1hex.g", get_bulk());
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex();

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_2block_2hex.g", get_bulk());
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_2block_2hex_withSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex_withSides();

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_2block_2hex_withSides.g", get_bulk());
  }
}

TEST_F(TestDisconnectBlocks, disconnect_2block_2cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet();

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_2block_2cubeOfTet.g", get_bulk());
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_2block_2cubeOfTet_withSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet_withSides();

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(16u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_2block_2cubeOfTet_withSides.g", get_bulk());
  }
}

TEST_F(TestDisconnectBlocks, disconnect_3block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex();

    output_mesh("disconnect_3block_4hex_init.g", get_bulk());

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_3block_4hex.g", get_bulk());
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_3block_4hex_withSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex_withSides();

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_3block_4hex_withSides.g", get_bulk());
  }
}

TEST_F(TestDisconnectBlocks, disconnect_3block_4cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet();

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_3block_4cubeOfTet.g", get_bulk());
  }
}

TEST_F(TestDisconnectBlocks, DISABLED_disconnect_3block_4cubeOfTet_withSides)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet_withSides();

    stk::tools::disconnect_all_blocks(get_bulk());

    EXPECT_EQ(0u,  get_num_intersecting_nodes(get_bulk(), blocks));
    EXPECT_EQ(28u, get_num_total_nodes(get_bulk()));
    EXPECT_FALSE(check_orphaned_nodes(get_bulk()));

    output_mesh("disconnect_3block_4cubeOfTet_withSides.g", get_bulk());
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

  output_mesh("disconnected_" + exodusFileName, bulk);
  double meshWriteTime = stk::wall_time();

  if (bulk.parallel_rank() == 0) {
    std::cout << " Mesh read time = " << (meshReadTime - startTime) << " s" << std::endl;
    std::cout << "Disconnect time = " << (disconnectTime - meshReadTime) << " s" << std::endl;
    std::cout << "Mesh write time = " << (meshWriteTime - disconnectTime) << " s" << std::endl;
  }
}

typedef TestDisconnectBlocks2D TestReconnectBlocks2D;
typedef TestDisconnectBlocks   TestReconnectBlocks;

void insert_block_pair(stk::mesh::Part* block1, stk::mesh::Part* block2,
                       std::vector<stk::tools::impl::BlockPairType>& blockPairs)
{
  ThrowRequire(nullptr != block1 && nullptr != block2);
  ThrowRequire(block1 != block2);

  if(block2->mesh_meta_data_ordinal() > block1->mesh_meta_data_ordinal()) {
    blockPairs.push_back(std::make_pair(block1, block2));
  } else {
    blockPairs.push_back(std::make_pair(block2, block1));
  }
}

void test_one_block_reconnect(stk::mesh::BulkData& bulk)
{
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  ThrowRequire(allBlocksInMesh.size() == 1u);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;
  info.debugLevel = get_debug_level();

  std::vector<stk::tools::impl::BlockPairType> blockPairs;

  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::disconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::reconnect_block_pairs(bulk, blockPairs, info);
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

  std::vector<stk::tools::impl::BlockPairType> blockPairs;

  insert_block_pair(blocks[0], blocks[1], blockPairs);

  EXPECT_EQ(expectedCommonNodeCount, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::disconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::reconnect_block_pairs(bulk, blockPairs, info);
  EXPECT_EQ(expectedCommonNodeCount, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  output_mesh(bulk);
}

TEST_F(TestReconnectBlocks2D, reconnect_1block_1quad)
{
  stk::mesh::PartVector blocks = setup_mesh_1block_1quad();
  test_one_block_reconnect(get_bulk());
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_2quad)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad();
  test_two_block_reconnect(get_bulk(), blocks, 2u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_2quad_diagonal)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_diagonal();
  test_two_block_reconnect(get_bulk(), blocks, 1u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_2quad_reversed)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_2quad_reversed();
  test_two_block_reconnect(get_bulk(), blocks, 2u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(1);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(2);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(3);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_corner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_corner(4);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(1);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(2);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(3);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_swappedCorner_decomp4)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_swappedCorner(4);
  test_two_block_reconnect(get_bulk(), blocks, 3u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(1);
  test_two_block_reconnect(get_bulk(), blocks, 5u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(2);
  test_two_block_reconnect(get_bulk(), blocks, 5u);
}

TEST_F(TestReconnectBlocks2D, reconnect_2block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_2block_4quad_checkerboard(3);
  test_two_block_reconnect(get_bulk(), blocks, 5u);
}

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

void test_reconnect_block_pairs(stk::mesh::BulkData& bulk, stk::mesh::PartVector& blocks,
                                const std::vector<BlockConnection>& reconnectPairs,
                                unsigned expectedGlobalInitialCommonNodes,
                                unsigned expectedGlobalConnectCommonNodes)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::PartVector allBlocksInMesh;
  stk::tools::impl::get_all_blocks_in_mesh(bulk, allBlocksInMesh);
  stk::tools::impl::LinkInfo info;
  info.preserveOrphans = true;
  info.debugLevel = get_debug_level();

  std::vector<stk::tools::impl::BlockPairType> disconnectBlockPairs;
  for(unsigned i = 0; i < allBlocksInMesh.size() - 1; i++) {
    for(unsigned j = i+1; j < allBlocksInMesh.size(); j++) {
      insert_block_pair(allBlocksInMesh[i], allBlocksInMesh[j], disconnectBlockPairs);
    }
  }

  std::vector<stk::tools::impl::BlockPairType> reconnectBlockPairs;
  for(const BlockConnection& connectPair : reconnectPairs) {
    ThrowRequire(connectPair.block1 != connectPair.block2);

    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    ThrowRequire(block1 != nullptr && block2 != nullptr);

    insert_block_pair(block1, block2, reconnectBlockPairs);
  }

  EXPECT_EQ(expectedGlobalInitialCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::disconnect_block_pairs(bulk, disconnectBlockPairs, info);
  output_mesh("test.g", bulk);
  EXPECT_EQ(0u, get_num_intersecting_nodes(bulk, allBlocksInMesh));
  stk::tools::reconnect_block_pairs(bulk, reconnectBlockPairs, info);
  EXPECT_EQ(expectedGlobalConnectCommonNodes, get_num_intersecting_nodes(bulk, allBlocksInMesh));

  stk::mesh::PartVector connectedBlockPair(2);
  for(const BlockConnection& connectPair : reconnectPairs) {
    stk::mesh::Part* block1 = meta.get_part(connectPair.block1);
    stk::mesh::Part* block2 = meta.get_part(connectPair.block2);

    connectedBlockPair[0] = block1;
    connectedBlockPair[1] = block2;

    EXPECT_EQ(connectPair.numExpectedIntersectingNodes, get_num_intersecting_nodes(bulk, connectedBlockPair));
  }
  output_mesh(bulk);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(1, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp1_v2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(1, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(1, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder1_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(1, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(2, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2), BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(2, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2), BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp2_v2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(2, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp2_v3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(2, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 2);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder2_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(2, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",2), BlockConnection("block_2","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder3_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(3, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_3","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder3_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(3, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_3","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder3_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(3, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_1",2), BlockConnection("block_3","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder4_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(4, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_2",2), BlockConnection("block_3","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder4_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(4, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_2",2), BlockConnection("block_3","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder4_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(4, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_3","block_2",2), BlockConnection("block_3","block_1",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder5_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(5, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder5_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(5, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder5_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(5, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",2), BlockConnection("block_2","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder6_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(6, 1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_3",2), BlockConnection("block_1","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder6_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(6, 2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_3",2), BlockConnection("block_1","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_blockOrder6_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad(6, 3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_3",2), BlockConnection("block_1","block_2",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_keepLowerRight_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_keepLowerRight_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_keepLowerRight_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_keepLowerRight(3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",2), BlockConnection("block_1","block_3",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 4, 3);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_checkerboard_decomp1)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",3), BlockConnection("block_1","block_3",3)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 5);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_checkerboard_decomp2)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(2);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",3), BlockConnection("block_1","block_3",3)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 5);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_checkerboard_decomp3)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_checkerboard(3);
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",3), BlockConnection("block_1","block_3",3)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 5);
}

TEST_F(TestReconnectBlocks2D, reconnect_3block_4quad_bowtie)
{
  stk::mesh::PartVector blocks = setup_mesh_3block_4quad_bowtie();
  BlockConnectionVector reconnectPairs{BlockConnection("block_1","block_2",1), BlockConnection("block_1","block_3",1)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 1, 1);
}

TEST_F(TestReconnectBlocks2D, reconnect_4block_4quad_blockorder1)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 3) { return; }

  stk::mesh::PartVector blocks = setup_mesh_4block_4quad(1);
  BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_4",2)};
  test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 5, 2);
}

TEST_F(TestReconnectBlocks, reconnect_1block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_1block_1hex();
    test_one_block_reconnect(get_bulk());
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_1hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 1) {
    stk::mesh::PartVector blocks = setup_mesh_2block_1hex();
    test_two_block_reconnect(get_bulk(), blocks, 0u);
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_2hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2hex();
    test_two_block_reconnect(get_bulk(), blocks, 4u);
  }
}

TEST_F(TestReconnectBlocks, reconnect_2block_2cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_2block_2cubeOfTet();
    test_two_block_reconnect(get_bulk(), blocks, 4u);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4hex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex();
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",4), BlockConnection("block_2","block_1",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 6);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4hex_v2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex();
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4hex_v3)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4hex();
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_1",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 4);
  }
}

TEST_F(TestReconnectBlocks, reconnect_3block_4cubeOfTet)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) <= 2) {
    stk::mesh::PartVector blocks = setup_mesh_3block_4cubeOfTet();
    BlockConnectionVector reconnectPairs{BlockConnection("block_2","block_3",4), BlockConnection("block_2","block_1",4)};
    test_reconnect_block_pairs(get_bulk(), blocks, reconnectPairs, 8, 6);
  }
}
