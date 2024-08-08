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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_THROW, etc
#include <stdexcept>                    // for runtime_error
#include <stk_unit_test_utils/ioUtils.hpp>  // for fill_mesh_using_stk_io
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <string>                       // for string
#include <cstdlib>                   // for unsetenv, setenv

#include "mpi.h"                        // for MPI_COMM_WORLD, MPI_Comm, etc
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_mesh/baseImpl/Partition.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc

namespace
{

std::shared_ptr<stk::mesh::BulkData> build_mesh(unsigned spatialDim,
                                                stk::ParallelMachine comm,
                                                stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA)
{
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(spatialDim);
  builder.set_aura_option(auraOption);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
  return bulk;
}

TEST(UnitTestChangeParts, test_throw_on_internal_part_change)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulkData = *bulkPtr;

  std::string generatedMeshSpec = "generated:1x1x4";
  stk::io::fill_mesh(generatedMeshSpec, bulkData);

  stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, 1);

  stk::mesh::PartVector addParts;
  stk::mesh::PartVector removeParts;

  addParts.push_back(&metaData.locally_owned_part());
  EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);

  addParts.clear();
  addParts.push_back(&metaData.globally_shared_part());
  EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);

  addParts.clear();
  removeParts.push_back(&metaData.locally_owned_part());
  EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);

  removeParts.clear();
  removeParts.push_back(&metaData.globally_shared_part());
  EXPECT_THROW(bulkData.change_entity_parts(node, addParts, removeParts), std::runtime_error);
}

TEST(UnitTestChangeParts, test_batch_part_change)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 1) {
    return;
  }

  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulkData = *bulkPtr;

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::mesh::Part& part = metaData.declare_part_with_topology("new_part", stk::topology::NODE);
  stk::unit_test_util::setup_text_mesh(bulkData, meshDesc);

  stk::mesh::Entity elem1 = bulkData.get_entity(stk::topology::ELEM_RANK, 1u);
  EXPECT_TRUE(bulkData.is_valid(elem1));

  stk::mesh::EntityVector nodes(bulkData.begin_nodes(elem1), bulkData.begin_nodes(elem1)+bulkData.num_nodes(elem1));
  EXPECT_EQ(8u, nodes.size());

  for(stk::mesh::Entity node : nodes) {
    EXPECT_FALSE(bulkData.bucket(node).member(part));
  }

  stk::mesh::PartVector add_parts(1, &part);
  bulkData.batch_change_entity_parts(nodes, add_parts, {});

  for(stk::mesh::Entity node : nodes) {
    EXPECT_TRUE(bulkData.bucket(node).member(part));
  }

  stk::mesh::PartVector remove_parts(1, &part);
  bulkData.batch_change_entity_parts(nodes, {}, remove_parts);

  for(stk::mesh::Entity node : nodes) {
    EXPECT_FALSE(bulkData.bucket(node).member(part));
  }
}

TEST(UnitTestChangeParts, test_superset_and_subset_part_change)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size( pm );

  if (p_size != 1) {
    return;
  }

  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, pm, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& metaData = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulkData = *bulkPtr;

  stk::mesh::Part& supersetPart = metaData.declare_part_with_topology("parent", stk::topology::NODE);
  stk::mesh::Part& subsetPart1  = metaData.declare_part_with_topology("child 1", stk::topology::NODE);
  stk::mesh::Part& subsetPart2  = metaData.declare_part_with_topology("child 2", stk::topology::NODE);

  metaData.declare_part_subset(supersetPart, subsetPart1);
  metaData.declare_part_subset(supersetPart, subsetPart2);

  stk::mesh::FieldBase &field = metaData.declare_field<double>(stk::topology::NODE_RANK, "sam");
  stk::mesh::put_field_on_mesh(field, supersetPart, nullptr);

  std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
  stk::unit_test_util::setup_text_mesh(bulkData, meshDesc);

  stk::mesh::Entity node1 = bulkData.get_entity(stk::topology::NODE_RANK, 1u);
  EXPECT_TRUE(bulkData.is_valid(node1));

  stk::mesh::Entity node2 = bulkData.get_entity(stk::topology::NODE_RANK, 2u);
  EXPECT_TRUE(bulkData.is_valid(node2));

  EXPECT_FALSE(bulkData.bucket(node1).member(supersetPart));
  EXPECT_FALSE(bulkData.bucket(node1).member(subsetPart1));
  EXPECT_FALSE(bulkData.bucket(node1).member(subsetPart2));

  EXPECT_FALSE(bulkData.bucket(node2).member(supersetPart));
  EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart1));
  EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart2));

  bulkData.modification_begin();
  stk::mesh::PartVector add_parts(1, &subsetPart1);
  bulkData.change_entity_parts(node1, add_parts, {});
  bulkData.modification_end();

  EXPECT_TRUE(bulkData.bucket(node1).member(supersetPart));
  EXPECT_TRUE(bulkData.bucket(node1).member(subsetPart1));
  EXPECT_FALSE(bulkData.bucket(node1).member(subsetPart2));

  bulkData.modification_begin();
  add_parts[0] = &supersetPart;
  bulkData.change_entity_parts(node2, add_parts, {});
  bulkData.modification_end();

  EXPECT_TRUE(bulkData.bucket(node2).member(supersetPart));
  EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart1));
  EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart2));

  stk::mesh::PartVector partVector{&subsetPart1, &subsetPart2};
  stk::mesh::Selector selector = supersetPart & !stk::mesh::selectUnion(partVector);
  stk::mesh::EntityVector nodes;
  stk::mesh::get_selected_entities(selector, bulkData.buckets(stk::topology::NODE_RANK), nodes);

  EXPECT_EQ(1u, nodes.size());
  EXPECT_EQ(2u, bulkData.identifier(nodes[0]));

  bulkData.modification_begin();
  add_parts[0] = &subsetPart2;
  bulkData.change_entity_parts(node2, add_parts, {});
  bulkData.modification_end();

  EXPECT_TRUE(bulkData.bucket(node2).member(supersetPart));
  EXPECT_FALSE(bulkData.bucket(node2).member(subsetPart1));
  EXPECT_TRUE(bulkData.bucket(node2).member(subsetPart2));

  double* node1Data = (double*) stk::mesh::field_data(field, node1);
  double* node2Data = (double*) stk::mesh::field_data(field, node1);

  EXPECT_TRUE(node1Data != nullptr);
  EXPECT_TRUE(node2Data != nullptr);


  for(unsigned i=3; i<8u; ++i) {
    stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, i);
    EXPECT_TRUE(bulkData.is_valid(node));
    double* nodeData = (double*) stk::mesh::field_data(field, node);
    EXPECT_TRUE(nodeData == nullptr);
  }
}

TEST(ChangeElemParts, DontMarkSharedNodesModifiedUnnecessarily)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if(numProcs > 2) { GTEST_SKIP(); }

  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::Part& myPart = meta.declare_part("myPart");
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh("generated:1x1x2", bulk);

  bulk.modification_begin();
  stk::mesh::Part& block1 = *meta.get_part("block_1");
  stk::mesh::PartVector parts = {&block1, &myPart};
  stk::mesh::EntityId elemId = bulk.parallel_rank() + 1;
  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
  bulk.change_entity_parts(elem, parts);
  bulk.modification_end();

  //elem should be "modified" because it was not already in myPart
  EXPECT_EQ(stk::mesh::Modified, bulk.state(elem));
  //node 5 should not be "modified" because it was already in block_1 and
  //shouldn't be induced into myPart
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_FALSE(bulk.bucket(node5).member(myPart));
  EXPECT_EQ(stk::mesh::Unchanged, bulk.state(node5));
}

TEST(ChangeElemParts, addThenRemoveElemPart_checkSharedNode)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if(numProcs > 2) { GTEST_SKIP(); }

  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::Part& myPart = meta.declare_part_with_topology("myPart", stk::topology::HEX_8);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh("generated:1x1x2", bulk);

  stk::mesh::Part& block1 = *meta.get_part("block_1");
  stk::mesh::PartVector parts = {&block1, &myPart}, empty;
  stk::mesh::EntityId elemId = bulk.parallel_rank() + 1;
  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, elemId);

  bulk.modification_begin();
  bulk.change_entity_parts(elem, parts, empty);
  bulk.modification_end();

  //elem should be "modified" because it was not already in myPart
  EXPECT_EQ(stk::mesh::Modified, bulk.state(elem));
  EXPECT_TRUE(bulk.bucket(elem).member(myPart));
  //node 5 should be "modified" because it was also added to myPart
  stk::mesh::Entity node5 = bulk.get_entity(stk::topology::NODE_RANK, 5);
  EXPECT_TRUE(bulk.bucket(node5).member(myPart));
  EXPECT_EQ(stk::mesh::Modified, bulk.state(node5));

  bulk.modification_begin();
  bulk.change_entity_parts(elem, empty, parts);
  bulk.modification_end();

  //elem should be "modified" because it is no longer in myPart
  EXPECT_EQ(stk::mesh::Modified, bulk.state(elem));
  EXPECT_FALSE(bulk.bucket(elem).member(myPart));
  //node 5 should be "modified" because it is no longer in myPart
  EXPECT_FALSE(bulk.bucket(node5).member(myPart));
  EXPECT_EQ(stk::mesh::Modified, bulk.state(node5));
}

class TestChangePartsWithSelector : public stk::unit_test_util::MeshFixture
{
public:
  TestChangePartsWithSelector()
    : stk::unit_test_util::MeshFixture(3)
  {
    setenv("STK_MESH_RUN_CONSISTENCY_CHECK", "ON", 1);
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void setup_n_unranked_parts_n_elements(unsigned numBlockParts, unsigned numElem, stk::mesh::PartVector& parts)
  {
    std::string meshSpec = "generated:1x1x" + std::to_string(numElem);

    for(unsigned i = 1; i < numBlockParts; i++) {
      std::string partName = "block_" + std::to_string(i+1);
      parts.push_back(&get_meta().declare_part(partName));
    }
    stk::io::fill_mesh(meshSpec, get_bulk());
    parts.insert(parts.begin(), get_meta().get_part("block_1"));
  }

  void setup_n_ranked_parts_n_elements(unsigned numBlockParts, unsigned numElem, stk::mesh::PartVector& parts, std::string opt = "")
  {
    std::string meshSpec = "generated:1x1x" + std::to_string(numElem) + opt;

    for(unsigned i = 1; i < numBlockParts; i++) {
      std::string partName = "block_" + std::to_string(i+1);
      parts.push_back(&get_meta().declare_part(partName, stk::topology::ELEM_RANK));
    }
    stk::io::fill_mesh(meshSpec, get_bulk());
    parts.insert(parts.begin(), get_meta().get_part("block_1"));
  }

  void setup_n_ranked_parts_n_elements_with_faces(unsigned numBlockParts, unsigned numElem, stk::mesh::PartVector& parts)
  {
    setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts, "|sideset:XxYyZz");
  }

  void insert_into_add_part(stk::mesh::Part* part)
  {
    addParts.push_back(part);
  }

  void insert_into_remove_part(stk::mesh::Part* part)
  {
    removeParts.push_back(part);
  }

  void clear_mod_part_vectors()
  {
    addParts.clear();
    removeParts.clear();
  }

  void change_entity_parts(stk::mesh::Entity entity)
  {
    get_bulk().modification_begin();
    EXPECT_NO_THROW(get_bulk().change_entity_parts(entity, addParts, removeParts));
    get_bulk().modification_end();

    clear_mod_part_vectors();
  }

  void batch_change_entity_parts(stk::mesh::Entity entity)
  {
    get_bulk().modification_begin();
    EXPECT_NO_THROW(get_bulk().change_entity_parts(entity, addParts, removeParts));
    get_bulk().modification_end();

    clear_mod_part_vectors();
  }

  void batch_change_entity_parts_with_selector(stk::mesh::Selector selector, stk::mesh::EntityRank rank)
  {
    EXPECT_NO_THROW(get_bulk().batch_change_entity_parts(selector, rank, addParts, removeParts));
    clear_mod_part_vectors();
  }

  void change_entity_parts_with_selector(stk::mesh::Selector selector, stk::mesh::EntityRank rank)
  {
    EXPECT_NO_THROW(get_bulk().modification_begin());
    EXPECT_NO_THROW(get_bulk().change_entity_parts(selector, rank, addParts, removeParts));
    EXPECT_NO_THROW(get_bulk().modification_end());
    clear_mod_part_vectors();
  }

  void test_entity_counts_in_parts(stk::mesh::Part* part, stk::mesh::EntityRank rank, unsigned expectedCountInPart)
  {
    unsigned entityCount = stk::mesh::count_entities(get_bulk(), rank, stk::mesh::Selector(*part));
    EXPECT_EQ(expectedCountInPart, entityCount) << " unexpected entity count in " << part->name() << " " << rank << std::endl;
  }

  void test_bucket_count(stk::mesh::EntityRank, unsigned expectedBucketCount)
  {
    unsigned bucketCount = get_bulk().get_buckets(stk::topology::ELEM_RANK, get_meta().universal_part()).size();
    EXPECT_EQ(expectedBucketCount, bucketCount);
  }

  void test_partitions_equality(stk::mesh::EntityRank rank, stk::mesh::Bucket const* bucket1, stk::mesh::Bucket const* bucket2)
  {
    stk::mesh::impl::Partition* partition1 = bucket1->getPartition();
    stk::mesh::impl::Partition* partition2 = bucket2->getPartition();

    EXPECT_EQ(partition1->size(), partition2->size());
    EXPECT_EQ(partition1->num_buckets(), partition2->num_buckets());

    const unsigned* key1 = partition1->key();
    const unsigned* key2 = partition2->key();

    EXPECT_EQ(key1[0], key2[0]);

    if(key1[0] == key2[0]) {
      for(unsigned i = 0; i < key1[0]; i++) {
        EXPECT_EQ(key1[i], key2[i]);
      }
    }

    for(auto it1 = partition1->begin(), it2 = partition2->begin();
        it1 != partition1->end(); ++it1, ++it2) {

      EXPECT_EQ((*it1)->bucket_id(), (*it2)->bucket_id());
    }
  }

  stk::mesh::PartVector addParts, removeParts;
};

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_remove_initial_part)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 1;
  unsigned numElem = 1;
  stk::mesh::PartVector parts;
  setup_n_unranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];

  insert_into_remove_part(block1Part);
  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_to_new_unranked_part)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 1;
  stk::mesh::PartVector parts;
  setup_n_unranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);
  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 0);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_to_new_unranked_part_in_two_procs)
{
  if(get_bulk().parallel_size() != 2) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_unranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);
  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 0);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_to_new_unranked_part_in_two_procs_non_batch)
{
  if(get_bulk().parallel_size() != 2) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_unranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);
  change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 0);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_to_new_ranked_part_in_two_procs)
{
  if(get_bulk().parallel_size() != 2) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);
  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 8);
}

TEST_F(TestChangePartsWithSelector, two_elements_3_blocks_move_each_to_different_block_in_serial)
{
  if(get_bulk().parallel_size() != 1) { GTEST_SKIP(); }

  unsigned numBlockParts = 3;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];
  stk::mesh::Part* block3Part = parts[2];

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);
  stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  batch_change_entity_parts(elem1);

  insert_into_add_part(block3Part);
  insert_into_remove_part(block1Part);
  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 8);
  test_entity_counts_in_parts(block3Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block3Part, stk::topology::NODE_RANK, 8);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_to_new_ranked_part_in_two_procs_non_batch)
{
  if(get_bulk().parallel_size() != 2) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);
  change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 8);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_add_new_unranked_part)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 1;
  stk::mesh::PartVector parts;
  setup_n_unranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);
  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 0);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_remove_unranked_part)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 1;
  stk::mesh::PartVector parts;
  setup_n_unranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  {
    insert_into_add_part(block2Part);

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elements);
    change_entity_parts(elements[0]);

    test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

    test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 0);
  }

  insert_into_remove_part(block2Part);
  batch_change_entity_parts_with_selector(stk::mesh::selectUnion(stk::mesh::PartVector{block1Part, block2Part}), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 0);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_add_new_ranked_part)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 1;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);

  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 8);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_add_new_ranked_part_non_batch)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 1;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);

  change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 8);
}

TEST_F(TestChangePartsWithSelector, element_in_with_faces_ranked_part_add_new_ranked_part_remove_ranked_part)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements_with_faces(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);

  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::FACE_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 2);
  test_entity_counts_in_parts(block2Part, stk::topology::FACE_RANK, 10);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 12);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_remove_ranked_part)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 1;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  {
    insert_into_add_part(block2Part);

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elements);
    change_entity_parts(elements[0]);

    test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

    test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 8);
  }

  insert_into_remove_part(block1Part);
  batch_change_entity_parts_with_selector(stk::mesh::selectUnion(stk::mesh::PartVector{block1Part, block2Part}), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 8);
}

TEST_F(TestChangePartsWithSelector, elements_in_different_ranked_parts_but_shared_nodes)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  {
    stk::mesh::Entity elemToMove = get_bulk().get_entity(stk::topology::ELEM_RANK, 2u);
    EXPECT_TRUE(get_bulk().is_valid(elemToMove));
    insert_into_add_part(block2Part);
    insert_into_remove_part(block1Part);
    change_entity_parts(elemToMove);

    test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

    test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 8);

    const stk::mesh::BucketVector& nodeBucketsUnion = get_bulk().get_buckets(stk::topology::NODE_RANK, stk::mesh::selectUnion(parts));
    EXPECT_EQ(3u, nodeBucketsUnion.size());

    const stk::mesh::BucketVector& nodeBucketsIntersection = get_bulk().get_buckets(stk::topology::NODE_RANK, stk::mesh::selectIntersection(parts));
    EXPECT_EQ(1u, nodeBucketsIntersection.size());
  }

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);
  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 2);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 12);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_add_new_to_ranked_part_no_conflicting_partitions)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  {
    insert_into_add_part(block2Part);
    insert_into_remove_part(block1Part);

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elements);
    change_entity_parts(elements[1]);

    test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

    test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 8);

    test_bucket_count(stk::topology::ELEM_RANK, 2u);
  }

  insert_into_add_part(block2Part);
  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_bucket_count(stk::topology::ELEM_RANK, 2u);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 2);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 12);
}

TEST_F(TestChangePartsWithSelector, element_in_ranked_part_move_to_unranked_part_conflicting_partitions)
{
  if(get_bulk().parallel_size() > 1) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_unranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  {
    insert_into_add_part(block2Part);
    insert_into_remove_part(block1Part);

    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elements);
    change_entity_parts(elements[1]);

    test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 8);

    test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 1);
    test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 0);

    test_bucket_count(stk::topology::ELEM_RANK, 2u);
  }

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);
  batch_change_entity_parts_with_selector(stk::mesh::Selector(*block1Part), stk::topology::ELEM_RANK);

  test_entity_counts_in_parts(block1Part, stk::topology::ELEM_RANK, 0);
  test_entity_counts_in_parts(block1Part, stk::topology::NODE_RANK, 0);

  test_entity_counts_in_parts(block2Part, stk::topology::ELEM_RANK, 2);
  test_entity_counts_in_parts(block2Part, stk::topology::NODE_RANK, 0);

  test_bucket_count(stk::topology::ELEM_RANK, 2u);

  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);
  test_partitions_equality(stk::topology::ELEM_RANK, buckets[0], buckets[1]);
}

TEST_F(TestChangePartsWithSelector, two_procs_with_different_selectors)
{
  if(get_bulk().parallel_size() != 2) { return; }

  unsigned numBlockParts = 2;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];

  insert_into_add_part(block2Part);
  insert_into_remove_part(block1Part);

  stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;

  if(get_bulk().parallel_rank() == 0) {
    EXPECT_ANY_THROW(get_bulk().batch_change_entity_parts(stk::mesh::Selector(*block1Part), rank, addParts, removeParts));
  } else {
    EXPECT_ANY_THROW(get_bulk().batch_change_entity_parts(stk::mesh::Selector(*block2Part), rank, addParts, removeParts));
  }
}

TEST_F(TestChangePartsWithSelector, two_procs_with_different_parts)
{
  if(get_bulk().parallel_size() != 2) { return; }

  unsigned numBlockParts = 3;
  unsigned numElem = 2;
  stk::mesh::PartVector parts;
  setup_n_ranked_parts_n_elements(numBlockParts, numElem, parts);

  stk::mesh::Part* block1Part = parts[0];
  stk::mesh::Part* block2Part = parts[1];
  stk::mesh::Part* block3Part = parts[2];

  insert_into_remove_part(block1Part);
  if(get_bulk().parallel_rank() == 0) {
    insert_into_add_part(block2Part);
  } else {
    insert_into_add_part(block3Part);
  }

  stk::mesh::EntityRank rank = stk::topology::ELEM_RANK;
  EXPECT_ANY_THROW(get_bulk().batch_change_entity_parts(stk::mesh::Selector(*block1Part), rank, addParts, removeParts));
}

}
