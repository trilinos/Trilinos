// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <gtest/gtest.h>

#include <algorithm>
#include <stk_io/FillMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/CompositeRank.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_topology/topology.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>

#include "Assembly.hpp"

namespace stk
{
namespace io
{
namespace unit_test
{

class GetAssemblyEntities : public stk::io::unit_test::Assembly
{
 protected:
  struct BucketTestData {
    BucketTestData() : topology(stk::topology::INVALID_TOPOLOGY) {}

    BucketTestData(stk::topology topology_, const std::vector<size_t>& entitiesPerBucket_ = {})
        : topology(topology_), entitiesPerBucket(entitiesPerBucket_)
    {
    }

    stk::topology topology;
    std::vector<size_t> entitiesPerBucket;
  };

  struct EntityTestData {
    EntityTestData() : numUniqueLeaves(0) {}

    EntityTestData(size_t numUniqueLeaves_, const stk::mesh::EntityIdVector& ids_)
        : numUniqueLeaves(numUniqueLeaves_), ids(ids_)
    {
    }

    size_t numUniqueLeaves;
    stk::mesh::EntityIdVector ids;
  };

  void verify_num_elements(size_t goldCount)
  {
    std::vector<size_t> counts;
    stk::mesh::count_entities(get_meta().universal_part(), get_bulk(), counts);
    EXPECT_EQ(goldCount, counts[stk::topology::ELEM_RANK]);
  }

  void verify_single_element(
      stk::mesh::EntityId elemId, stk::topology topology, const stk::mesh::EntityIdVector& nodeIds)
  {
    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elemId);
    EXPECT_TRUE(get_bulk().is_valid(element));
    EXPECT_EQ(topology, get_bulk().bucket(element).topology());
    verify_nodes_on_element(element, nodeIds);
  }

  void verify_nodes_on_element(stk::mesh::Entity element, const stk::mesh::EntityIdVector& goldNodeIds)
  {
    stk::mesh::EntityVector nodes(get_bulk().begin_nodes(element), get_bulk().end_nodes(element));
    EXPECT_EQ(goldNodeIds, get_node_ids(nodes));
  }

  stk::mesh::EntityIdVector get_node_ids(const stk::mesh::EntityVector& nodes)
  {
    stk::mesh::EntityIdVector nodeIds;
    for (const stk::mesh::Entity& node : nodes) {
      nodeIds.emplace_back(get_bulk().identifier(node));
    }
    return nodeIds;
  }

  void test_get_buckets(stk::mesh::Part& assemblyPart, stk::mesh::EntityRank rank, const BucketTestData& gold)
  {
    stk::mesh::PartVector parts = stk::io::get_unique_leaf_parts(get_meta(), assemblyPart.name());
    EXPECT_EQ(1u, parts.size());

    stk::mesh::Selector assemblySelector(assemblyPart);
    stk::mesh::BucketVector buckets = get_bulk().get_buckets(rank, assemblySelector);
    EXPECT_EQ(gold.entitiesPerBucket.size(), buckets.size());

    for (size_t i = 0; i < buckets.size(); ++i) {
      auto bucket = buckets[i];
      EXPECT_EQ(gold.topology, bucket->topology());
      EXPECT_EQ(gold.entitiesPerBucket[i], bucket->size());
      EXPECT_TRUE(bucket->member_any(parts));
    }
  }

  stk::mesh::EntityVector get_sorted_entities_from_ids(stk::mesh::EntityRank rank, const stk::mesh::EntityIdVector& ids)
  {
    stk::mesh::EntityVector entities;
    entities.reserve(ids.size());

    for (auto id : ids) {
      stk::mesh::Entity entity = get_bulk().get_entity(rank, id);
      EXPECT_TRUE(get_bulk().is_valid(entity)) << "Invalid entity id: " << id << " of rank: " << rank;
      entities.push_back(entity);
    }
    std::sort(entities.begin(), entities.end(), stk::mesh::EntityLess(get_bulk()));

    return entities;
  }

  void compare_selected_entities(
      stk::mesh::Part& assemblyPart, stk::mesh::EntityRank rank, const stk::mesh::EntityIdVector& goldIds)
  {
    stk::mesh::EntityVector goldEntities = get_sorted_entities_from_ids(rank, goldIds);

    stk::mesh::Selector assemblySelector(assemblyPart);
    stk::mesh::EntityVector entities;
    stk::mesh::get_entities(get_bulk(), rank, assemblySelector, entities);
    EXPECT_EQ(goldEntities.size(), entities.size());

    for (auto entity : entities) {
      stk::mesh::EntityKey key = get_bulk().entity_key(entity);
      EXPECT_TRUE(
          std::binary_search(goldEntities.begin(), goldEntities.end(), entity, stk::mesh::EntityLess(get_bulk())))
          << "Could not find selected entity: " << key << " in gold list";
    }
  }

  void test_get_entities(stk::mesh::Part& assemblyPart, stk::mesh::EntityRank rank, const EntityTestData& gold)
  {
    stk::mesh::PartVector parts = stk::io::get_unique_leaf_parts(get_meta(), assemblyPart.name());
    EXPECT_EQ(gold.numUniqueLeaves, parts.size());
    compare_selected_entities(assemblyPart, rank, gold.ids);
  }

  stk::mesh::Part& setup_single_particle_mesh()
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

    const std::string assemblyName("simpleAssembly");
    stk::mesh::Part& assemblyPart = create_single_particle_assembly(assemblyName, 10);
    create_single_particle_mesh();
    return assemblyPart;
  }

  void create_single_particle_mesh()
  {
    std::string meshDesc = "textmesh:0,1,PARTICLE,1,block_1";
    stk::io::fill_mesh(meshDesc, get_bulk());
  }

  stk::mesh::Part& create_single_particle_assembly(const std::string& assemblyName, int id)
  {
    const std::vector<std::string> partNames{"block_1"};

    stk::mesh::Part& assemblyPart = create_assembly(assemblyName, id);
    stk::mesh::Part& block1Part = create_io_part(partNames[0], 1, stk::topology::PARTICLE);
    declare_subsets(assemblyPart, {&block1Part});

    return assemblyPart;
  }

  stk::mesh::Part& setup_single_hex_mesh(stk::mesh::EntityRank assemblyRank = stk::topology::INVALID_RANK)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

    const std::string assemblyName("simpleAssembly");
    stk::mesh::Part& assemblyPart = create_single_hex_assembly(assemblyName, 10, assemblyRank);
    create_single_hex_mesh();
    return assemblyPart;
  }

  void create_single_hex_mesh()
  {
    std::string meshDesc = "textmesh:0,1,HEX_8,1,2,3,4,5,6,7,8,block_1";
    stk::io::fill_mesh(meshDesc, get_bulk());
  }

  stk::mesh::Part& create_single_hex_assembly(const std::string& assemblyName, int id,
                                              stk::mesh::EntityRank assemblyRank = stk::topology::INVALID_RANK)
  {
    const std::vector<std::string> partNames{"block_1"};

    stk::mesh::Part& assemblyPart = create_assembly(assemblyName, id, assemblyRank);
    stk::mesh::Part& block1Part = create_io_part(partNames[0], 1);
    declare_subsets(assemblyPart, {&block1Part});

    return assemblyPart;
  }

  stk::mesh::Part& setup_two_hex_mesh_with_block_in_assembly(const std::string& assemblyBlock)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

    const std::string assemblyName("simpleAssembly");
    stk::mesh::Part& assemblyPart = create_two_hex_blocks_and_assembly_from_block(assemblyName, 10, assemblyBlock);
    create_two_hex_mesh();

    return assemblyPart;
  }

  void create_two_hex_mesh()
  {
    std::string meshDesc =
        "textmesh:"
        "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n"
        "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
    stk::io::fill_mesh(meshDesc, get_bulk());

    verify_num_elements(2);
    verify_single_element(1, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2, stk::topology::HEX_8, stk::mesh::EntityIdVector{5, 6, 7, 8, 9, 10, 11, 12});
  }

  stk::mesh::Part& create_two_hex_blocks_and_assembly_from_block(
      const std::string& assemblyName, int id, const std::string& assemblyBlock)
  {
    return create_two_blocks_and_assembly_from_block(
        assemblyName, id, assemblyBlock, {stk::topology::HEX_8, stk::topology::HEX_8});
  }

  stk::mesh::Part& setup_hex_pyramid_mesh_with_block_in_assembly(const std::string& assemblyBlock)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

    const std::string assemblyName("simpleAssembly");
    stk::mesh::Part& assemblyPart = create_hex_pyramid_blocks_and_assembly_from_block(assemblyName, 10, assemblyBlock);
    create_hex_pyramid_mesh();

    return assemblyPart;
  }

  void create_hex_pyramid_mesh()
  {
    std::string meshDesc =
        "textmesh:"
        "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
        "0,2,PYRAMID_5,5,6,7,8,9,block_2";
    stk::io::fill_mesh(meshDesc, get_bulk());

    verify_num_elements(2);
    verify_single_element(1, stk::topology::HEX_8, stk::mesh::EntityIdVector{1, 2, 3, 4, 5, 6, 7, 8});
    verify_single_element(2, stk::topology::PYRAMID_5, stk::mesh::EntityIdVector{5, 6, 7, 8, 9});
  }

  stk::mesh::Part& create_hex_pyramid_blocks_and_assembly_from_block(
      const std::string& assemblyName, int id, const std::string& assemblyBlock)
  {
    return create_two_blocks_and_assembly_from_block(
        assemblyName, id, assemblyBlock, {stk::topology::HEX_8, stk::topology::PYRAMID_5});
  }

  stk::mesh::Part& create_two_blocks_and_assembly_from_block(const std::string& assemblyName,
      int id,
      const std::string& assemblyBlock,
      const std::vector<stk::topology>& blockTopologies)
  {
    const std::vector<std::string> partNames{"block_1", "block_2"};
    assert((assemblyBlock == partNames[0]) || (assemblyBlock == partNames[1]));
    assert(blockTopologies.size() == 2);

    stk::mesh::Part& assemblyPart = create_assembly(assemblyName, id);
    stk::mesh::Part& block1Part = create_io_part(partNames[0], 1, blockTopologies[0]);
    stk::mesh::Part& block2Part = create_io_part(partNames[1], 2, blockTopologies[1]);

    if (assemblyBlock == partNames[0]) {
      declare_subsets(assemblyPart, {&block1Part});
    } else if (assemblyBlock == partNames[1]) {
      declare_subsets(assemblyPart, {&block2Part});
    }

    return assemblyPart;
  }
};

TEST_F(GetAssemblyEntities, noEntities)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  const std::string assemblyName("myAssembly");
  stk::mesh::Part& assemblyPart = create_assembly(assemblyName, 100);

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::ELEM_RANK, EntityTestData(0, {}));
  test_get_entities(assemblyPart, stk::topology::FACE_RANK, EntityTestData(0, {}));
  test_get_entities(assemblyPart, stk::topology::EDGE_RANK, EntityTestData(0, {}));
  test_get_entities(assemblyPart, stk::topology::NODE_RANK, EntityTestData(0, {}));
}

TEST_F(GetAssemblyEntities, singleParticle_getElementBucket)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  stk::mesh::Part& assemblyPart = setup_single_particle_mesh();

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::ELEM_RANK, BucketTestData(stk::topology::PARTICLE, {1}));
}

TEST_F(GetAssemblyEntities, singleParticle_getElement)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  stk::mesh::Part& assemblyPart = setup_single_particle_mesh();

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::ELEM_RANK, EntityTestData(1, {1}));
}

TEST_F(GetAssemblyEntities, singleParticle_getNodeBucket)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  stk::mesh::Part& assemblyPart = setup_single_particle_mesh();

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::NODE_RANK, BucketTestData(stk::topology::NODE, {1}));
}

TEST_F(GetAssemblyEntities, singleParticle_getNode)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  stk::mesh::Part& assemblyPart = setup_single_particle_mesh();

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::NODE_RANK, EntityTestData(1, {1}));
}

TEST_F(GetAssemblyEntities, singleHex_getElementBucket)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  stk::mesh::Part& assemblyPart = setup_single_hex_mesh();

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::ELEM_RANK, BucketTestData(stk::topology::HEX_8, {1}));
}

TEST_F(GetAssemblyEntities, singleHex_getElementBucket_rankedAssemblyPart)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  stk::mesh::Part& assemblyPart = setup_single_hex_mesh(stk::topology::ELEM_RANK);

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::ELEM_RANK, BucketTestData(stk::topology::HEX_8, {1}));
}

TEST_F(GetAssemblyEntities, singleHex_getElement)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  stk::mesh::Part& assemblyPart = setup_single_hex_mesh();

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::ELEM_RANK, EntityTestData(1, {1}));
}

TEST_F(GetAssemblyEntities, singleHex_getNodeBucket)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  stk::mesh::Part& assemblyPart = setup_single_hex_mesh();

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::NODE_RANK, BucketTestData(stk::topology::NODE, {8}));
}

TEST_F(GetAssemblyEntities, singleHex_getNodes)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  stk::mesh::Part& assemblyPart = setup_single_hex_mesh();

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::NODE_RANK, EntityTestData(1, {1, 2, 3, 4, 5, 6, 7, 8}));
}

TEST_F(GetAssemblyEntities, twoHexDifferentBlocks_getElementBucket)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  const std::string assemblyBlock("block_1");
  stk::mesh::Part& assemblyPart = setup_two_hex_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::ELEM_RANK, BucketTestData(stk::topology::HEX_8, {1}));
}

TEST_F(GetAssemblyEntities, twoHexDifferentBlocks_getElement)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  const std::string assemblyBlock("block_1");
  stk::mesh::Part& assemblyPart = setup_two_hex_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::ELEM_RANK, EntityTestData(1, {1}));
}

TEST_F(GetAssemblyEntities, twoHexDifferentBlocks_getNodeBuckets)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_1");
  stk::mesh::Part& assemblyPart = setup_two_hex_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::NODE_RANK, BucketTestData(stk::topology::NODE, {4, 4}));
}

TEST_F(GetAssemblyEntities, twoHexDifferentBlocks_getNodes)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_2");
  stk::mesh::Part& assemblyPart = setup_two_hex_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::NODE_RANK, EntityTestData(1, {5, 6, 7, 8, 9, 10, 11, 12}));
}

TEST_F(GetAssemblyEntities, hexPyramidMesh_getElementBucketOnHex)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_1");
  stk::mesh::Part& assemblyPart = setup_hex_pyramid_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::ELEM_RANK, BucketTestData(stk::topology::HEX_8, {1}));
}

TEST_F(GetAssemblyEntities, hexPyramidMesh_getElementBucketOnPyramid)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_2");
  stk::mesh::Part& assemblyPart = setup_hex_pyramid_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::ELEM_RANK, BucketTestData(stk::topology::PYRAMID_5, {1}));
}

TEST_F(GetAssemblyEntities, hexPyramidMesh_getElementOnHex)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_1");
  stk::mesh::Part& assemblyPart = setup_hex_pyramid_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::ELEM_RANK, EntityTestData(1, {1}));
}

TEST_F(GetAssemblyEntities, hexPyramidMesh_getElementOnPyramid)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_2");
  stk::mesh::Part& assemblyPart = setup_hex_pyramid_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::ELEM_RANK, EntityTestData(1, {2}));
}

TEST_F(GetAssemblyEntities, hexPyramidMesh_getNodeBucketsOnHex)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_1");
  stk::mesh::Part& assemblyPart = setup_hex_pyramid_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::NODE_RANK, BucketTestData(stk::topology::NODE, {4, 4}));
}

TEST_F(GetAssemblyEntities, hexPyramidMesh_getNodeBucketsOnPyramid)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_2");
  stk::mesh::Part& assemblyPart = setup_hex_pyramid_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_buckets(assemblyPart, stk::topology::NODE_RANK, BucketTestData(stk::topology::NODE, {1, 4}));
}

TEST_F(GetAssemblyEntities, hexPyramidMesh_getNodesOnHex)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_1");
  stk::mesh::Part& assemblyPart = setup_hex_pyramid_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::NODE_RANK, EntityTestData(1, {1, 2, 3, 4, 5, 6, 7, 8}));
}

TEST_F(GetAssemblyEntities, hexPyramidMesh_getNodesOnPyramid)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  const std::string assemblyBlock("block_2");
  stk::mesh::Part& assemblyPart = setup_hex_pyramid_mesh_with_block_in_assembly(assemblyBlock);

  test_assembly_part_attributes(assemblyPart);
  test_get_entities(assemblyPart, stk::topology::NODE_RANK, EntityTestData(1, {5, 6, 7, 8, 9}));
}

//  (ID) = nodeset assembly nodes
//  [ID] = sideset assembly faces
//  {ID} = element block assembly elements
//
//                      3--------------7------------(11)-------------15
//                     /|             /|             /|             /|
//                    / |            / |            / |            / |
//                   /  |           /  |           /  |           /  |
//                  4--------------8------------(12)-------------16  |
//                  |   |          |   |          |   |          |   |
//             [15] |   |  {1}     |   |   2      |   |   3      |   | [36]
//                  |   |          |   |          |   |          |   |
//                  |   2----------|---6----------|-(10)---------|---14
//     Y            |  /           |  /           |  /           |  /
//     |  X         | /            | /            | /            | /
//     | /          |/             |/             |/             |/
//     |/           1--------------5-------------(9)-------------13
//     ------Z

TEST_F(GetAssemblyEntities, heterogeneousAssemblies)
{
  if (stk::parallel_machine_size(get_comm()) != 1) {
    return;
  }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  const std::string assemblyName("simpleAssembly");

  stk::mesh::Part& block1Assembly = create_assembly("block1Assembly", 10);
  stk::mesh::Part& surface1Assembly = create_assembly("surface1Assembly", 11);
  stk::mesh::Part& surface2Assembly = create_assembly("surface2Assembly", 12);
  stk::mesh::Part& surface1And2Assembly = create_assembly("surface1And2Assembly", 13);
  stk::mesh::Part& block1AndSurface1And2Assembly = create_assembly("block1AndSurface1And2Assembly", 14);
  stk::mesh::Part& nodeset1Assembly = create_assembly("nodeset1Assembly", 15);
  stk::mesh::Part& block1AndSurface1And2AndNodeset1Assembly =
      create_assembly("block1AndSurface1And2AndNodeset1Assembly", 16);

  stk::mesh::Part& block1 = create_io_part("block_1", 1, stk::topology::HEX_8);
  stk::mesh::Part& block2 = create_io_part("block_2", 2, stk::topology::HEX_8);
  stk::mesh::Part& surface1 = create_io_part("surface_1", 1, stk::topology::QUAD_4);
  stk::mesh::Part& surface2 = create_io_part("surface_2", 2, stk::topology::QUAD_4);
  stk::mesh::Part& nodeset1 = create_io_part("nodelist_1", 1, stk::topology::NODE);

  stk::io::fill_mesh("generated:1x1x3|sideset:zZ", get_bulk());

  move_element(2, block1, block2);
  move_element(3, block1, block2);

  add_nodes({9, 10, 11, 12}, nodeset1);

  declare_subsets(block1Assembly, {&block1});
  declare_subsets(surface1Assembly, {&surface1});
  declare_subsets(surface2Assembly, {&surface2});
  declare_subsets(surface1And2Assembly, {&surface1Assembly, &surface2Assembly});
  declare_subsets(block1AndSurface1And2Assembly, {&block1Assembly, &surface1Assembly, &surface2Assembly});
  declare_subsets(nodeset1Assembly, {&nodeset1});
  declare_subsets(
      block1AndSurface1And2AndNodeset1Assembly, {&block1Assembly, &surface1Assembly, &surface2Assembly, &nodeset1});

  EXPECT_EQ(stk::topology::ELEM_RANK, stk::mesh::CompositeRank::get_rank(block1Assembly));
  EXPECT_EQ(stk::topology::FACE_RANK, stk::mesh::CompositeRank::get_rank(surface1Assembly));
  EXPECT_EQ(stk::topology::FACE_RANK, stk::mesh::CompositeRank::get_rank(surface2Assembly));
  EXPECT_EQ(stk::topology::FACE_RANK, stk::mesh::CompositeRank::get_rank(surface1And2Assembly));
  EXPECT_EQ(stk::topology::INVALID_RANK, stk::mesh::CompositeRank::get_rank(block1AndSurface1And2Assembly));
  EXPECT_EQ(stk::topology::NODE_RANK, stk::mesh::CompositeRank::get_rank(nodeset1Assembly));
  EXPECT_EQ(stk::topology::INVALID_RANK, stk::mesh::CompositeRank::get_rank(block1AndSurface1And2AndNodeset1Assembly));

  test_assembly_part_attributes(block1Assembly);
  test_assembly_part_attributes(surface1Assembly);
  test_assembly_part_attributes(surface2Assembly);
  test_assembly_part_attributes(surface1And2Assembly);
  test_assembly_part_attributes(block1AndSurface1And2Assembly);
  test_assembly_part_attributes(nodeset1Assembly);
  test_assembly_part_attributes(block1AndSurface1And2AndNodeset1Assembly);

  compare_selected_entities(block1Assembly, stk::topology::ELEM_RANK, {1});
  compare_selected_entities(block1Assembly, stk::topology::FACE_RANK, {15});
  compare_selected_entities(block1Assembly, stk::topology::EDGE_RANK, {});
  compare_selected_entities(block1Assembly, stk::topology::NODE_RANK, {1, 2, 3, 4, 5, 6, 7, 8});

  compare_selected_entities(surface1Assembly, stk::topology::FACE_RANK, {15});
  compare_selected_entities(surface1Assembly, stk::topology::EDGE_RANK, {});
  compare_selected_entities(surface1Assembly, stk::topology::NODE_RANK, {1, 2, 3, 4});

  compare_selected_entities(surface2Assembly, stk::topology::FACE_RANK, {36});
  compare_selected_entities(surface2Assembly, stk::topology::EDGE_RANK, {});
  compare_selected_entities(surface2Assembly, stk::topology::NODE_RANK, {13, 14, 15, 16});

  compare_selected_entities(surface1And2Assembly, stk::topology::FACE_RANK, {15, 36});
  compare_selected_entities(surface1And2Assembly, stk::topology::EDGE_RANK, {});
  compare_selected_entities(surface1And2Assembly, stk::topology::NODE_RANK, {1, 2, 3, 4, 13, 14, 15, 16});

  compare_selected_entities(block1AndSurface1And2Assembly, stk::topology::ELEM_RANK, {1});
  compare_selected_entities(block1AndSurface1And2Assembly, stk::topology::FACE_RANK, {15, 36});
  compare_selected_entities(block1AndSurface1And2Assembly, stk::topology::EDGE_RANK, {});
  compare_selected_entities(
      block1AndSurface1And2Assembly, stk::topology::NODE_RANK, {1, 2, 3, 4, 5, 6, 7, 8, 13, 14, 15, 16});

  compare_selected_entities(nodeset1Assembly, stk::topology::NODE_RANK, {9, 10, 11, 12});

  compare_selected_entities(block1AndSurface1And2AndNodeset1Assembly, stk::topology::ELEM_RANK, {1});
  compare_selected_entities(block1AndSurface1And2AndNodeset1Assembly, stk::topology::FACE_RANK, {15, 36});
  compare_selected_entities(block1AndSurface1And2AndNodeset1Assembly, stk::topology::EDGE_RANK, {});
  compare_selected_entities(block1AndSurface1And2AndNodeset1Assembly, stk::topology::NODE_RANK,
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16});
}

}  // namespace unit_test
}  // namespace io
}  // namespace stk
