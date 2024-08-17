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

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <string.h>                     // for memcpy, memmove
#include <algorithm>                    // for binary_search, sort
#include <iostream>                     // for basic_ostream::operator<<, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/util/AdjustForAlignment.hpp>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/FieldBase.hpp>  // for FieldMetaDataVector, etc
#include <stk_mesh/base/FieldDataManager.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, entity_rank_names, etc
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/environment/CPUTime.hpp>  // for cpu_time
#include <stk_util/util/string_utils.hpp>
#include <vector>                       // for vector, vector<>::iterator
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityId, etc

namespace stk { namespace mesh { class Part; } }

namespace
{

constexpr const char * fieldPrefix {"testField"};

std::string field_name(int fieldNumber) {
  return std::string(fieldPrefix) + std::to_string(fieldNumber);
}

bool is_test_field(const stk::mesh::FieldBase & field) {
  return stk::string_starts_with(field.name(), fieldPrefix);
}

stk::mesh::FieldVector extract_test_fields(const stk::mesh::FieldVector & allFields) {
  stk::mesh::FieldVector fields;
  std::copy_if(allFields.begin(), allFields.end(), std::back_inserter(fields),
               [&](stk::mesh::FieldBase * f) { return is_test_field(*f); });
  return fields;
}

void createPart(stk::mesh::MetaData& meshMetaData)
{
  meshMetaData.declare_part("part1");
}

template <typename T>
void initializeTestField(stk::mesh::MetaData& meshMetaData)
{
  stk::mesh::Field<T> &field1 = meshMetaData.declare_field<T>(stk::topology::NODE_RANK, field_name(1));
  T initial_value1 = 13;
  stk::mesh::put_field_on_entire_mesh_with_initial_value(field1, &initial_value1);
  stk::mesh::Field<T> &field2 = meshMetaData.declare_field<T>(stk::topology::NODE_RANK, field_name(2));
  T initial_value2 = 4;
  stk::mesh::put_field_on_entire_mesh_with_initial_value(field2, &initial_value2);
  stk::mesh::Field<T> &field3 = meshMetaData.declare_field<T>(stk::topology::NODE_RANK, field_name(3));
  T initial_value3[3] = {1, 2, 3};
  stk::mesh::put_field_on_mesh(field3, meshMetaData.universal_part(), 3, initial_value3);
  meshMetaData.commit();
}

template <typename T>
void testAllocateFieldData(stk::mesh::BulkData& bulkData, stk::mesh::FieldDataManager * fieldDataManager,
                           const size_t extraCapacityInBytes, const size_t numNodes)
{
  const stk::mesh::MetaData& meshMetaData = bulkData.mesh_meta_data();

  bulkData.deactivate_field_updating();
  bulkData.modification_begin();
  std::vector<stk::mesh::EntityId> node_ids(numNodes);
  for (size_t i = 0; i < node_ids.size(); i++) {
    node_ids[i] = stk::mesh::EntityId(i + 1);
    bulkData.declare_node(node_ids[i]);
  }
  bulkData.modification_end();
  bulkData.allocate_field_data();

  const stk::mesh::FieldVector &fields = meshMetaData.get_fields();
  for (size_t field_index = 0; field_index < fields.size(); field_index++) {
    const stk::mesh::FieldBase &field = *fields[field_index];
    if (is_test_field(field)) {
      const T *initial_value = reinterpret_cast<const T*>(field.get_initial_value());

      size_t totalBytesAllocatedForField = 0;
      for (const stk::mesh::Bucket * bucket : bulkData.buckets(stk::topology::NODE_RANK)) {
        size_t bytesPerEntity = field.get_meta_data_for_field()[bucket->bucket_id()].m_bytesPerEntity;
        size_t numEntitiesAllocated = (dynamic_cast<stk::mesh::DefaultFieldDataManager*>(fieldDataManager)) ? bucket->capacity()
                                                                                                            : bucket->size();
        totalBytesAllocatedForField += stk::adjust_up_to_alignment_boundary(numEntitiesAllocated*bytesPerEntity,
                                                                            fieldDataManager->get_alignment_bytes());
        size_t bucketCapacityOrSize = bucket->size();
        T* field_data_ptr = reinterpret_cast<T *>(stk::mesh::field_data(field, *bucket));
        for (size_t i=0; i<bucketCapacityOrSize; i++) {
          unsigned field_max_size = field.max_size();
          for (unsigned k=0; k<field_max_size; k++) {
            EXPECT_EQ(initial_value[k], field_data_ptr[field_max_size*i+k]);
          }
        }
      }
      EXPECT_EQ(totalBytesAllocatedForField + extraCapacityInBytes,
                fieldDataManager->get_num_bytes_allocated_on_field(field.mesh_meta_data_ordinal()));
    }
  }
}

template <typename T>
void testReorderBucketFieldData(stk::mesh::BulkData& bulkData, stk::mesh::FieldDataManager * fieldDataManager,
                                stk::mesh::EntityRank rank, const stk::mesh::FieldVector& fields,
                                const std::vector<unsigned> &reorderedBucketIds)
{
  const stk::mesh::BucketVector& buckets = bulkData.buckets(rank);
  for (size_t b=0; b<buckets.size(); ++b) {
    const stk::mesh::Bucket& bucket = *buckets[b];
    for (size_t i=0; i<fields.size(); ++i) {
      T value = (b+1)*1000 + (i+1)*100;

      T* data = reinterpret_cast<T*>(fields[i]->get_meta_data_for_field()[b].m_data);
      unsigned field_size = fields[i]->max_size();
      for (size_t offset_into_bucket=0; offset_into_bucket<bucket.size(); ++offset_into_bucket) {
        for (unsigned j=0; j<field_size; ++j) {
          data[offset_into_bucket*field_size + j] = value;
        }
      }
    }
  }

  fieldDataManager->reorder_bucket_field_data(rank, fields, reorderedBucketIds);

  for (size_t bucketIndex=0; bucketIndex<buckets.size(); ++bucketIndex) {
    size_t oldBucketIndex = reorderedBucketIds[bucketIndex];
    size_t oldBucketSize = buckets[oldBucketIndex]->size();

    for (size_t i=0; i<fields.size(); ++i) {
      T expected_value = (oldBucketIndex+1)*1000 + (i+1)*100;

      T* dataReorderedBucket = reinterpret_cast<T*>(fields[i]->get_meta_data_for_field()[bucketIndex].m_data);
      unsigned field_size = fields[i]->max_size();
      for (size_t offset_into_bucket=0; offset_into_bucket<oldBucketSize; ++offset_into_bucket) {
        for (unsigned j=0; j<field_size; ++j) {
          EXPECT_EQ(expected_value, dataReorderedBucket[offset_into_bucket*field_size + j]);
        }
      }
    }
  }
}

void testTwoEntitiesTwoBuckets(stk::mesh::BulkData &bulkData, stk::mesh::FieldDataManager * fieldDataManager)
{
  bulkData.deactivate_field_updating();

  // =======
  bulkData.modification_begin();

  const stk::mesh::MetaData& meshMetaData = bulkData.mesh_meta_data();
  stk::mesh::Part& part1 = *meshMetaData.get_part("part1");

  stk::mesh::EntityId nodeId1 = 1;
  stk::mesh::EntityId nodeId2 = 2;

  stk::mesh::Entity node1 = bulkData.declare_node(nodeId1, stk::mesh::ConstPartVector{&meshMetaData.universal_part()});
  stk::mesh::Entity node2 = bulkData.declare_node(nodeId2, stk::mesh::ConstPartVector{&part1});

  bulkData.modification_end();

  //=======

  bulkData.allocate_field_data();

  size_t expectedNumBuckets = 2;
  EXPECT_EQ(expectedNumBuckets, bulkData.buckets(stk::topology::NODE_RANK).size());

  //=======

  bulkData.modification_begin();

  stk::mesh::PartVector addToParts;//empty
  stk::mesh::PartVector rmFromParts(1, &part1);

  bulkData.change_entity_parts(node2, addToParts, rmFromParts); // node2 from part1 to universal_part

  bulkData.modification_end();

  //=======

  expectedNumBuckets = 1;
  EXPECT_EQ(expectedNumBuckets, bulkData.buckets(stk::topology::NODE_RANK).size());

  //=======

  bulkData.modification_begin();

  bulkData.destroy_entity(node1);
  bulkData.destroy_entity(node2);

  bulkData.modification_end();

  //=======

  expectedNumBuckets = 0;
  EXPECT_EQ(expectedNumBuckets, bulkData.buckets(stk::topology::NODE_RANK).size());
}

std::shared_ptr<stk::mesh::BulkData> build_mesh(unsigned spatialDim,
                                                stk::ParallelMachine comm,
                                                std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager)
{
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(spatialDim);
  builder.set_field_data_manager(std::move(fieldDataManager));
  return builder.create();
}

using TestTypes = ::testing::Types<char, unsigned char, signed char, short, unsigned short, int, unsigned int,
long, unsigned long, long long, unsigned long long, float, double, long double,
std::complex<float>, std::complex<double>>;

template <typename T>
class TestDefaultFieldDataManager : public testing::Test {};
TYPED_TEST_SUITE(TestDefaultFieldDataManager, TestTypes,);

TYPED_TEST(TestDefaultFieldDataManager, AllocateFieldData)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  const size_t numRanks = stk::mesh::entity_rank_names().size();
  auto fieldDataManager = std::make_unique<stk::mesh::DefaultFieldDataManager>(numRanks);
  auto * localFieldDataManager = fieldDataManager.get();
  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD, std::move(fieldDataManager));
  stk::mesh::MetaData& meshMetaData = bulkDataPtr->mesh_meta_data();
  initializeTestField<FieldDataType>(meshMetaData);

  size_t numNodes = 20;
  size_t extraCapacity = 0;
  testAllocateFieldData<FieldDataType>(*bulkDataPtr, localFieldDataManager, extraCapacity, numNodes);
}

TYPED_TEST(TestDefaultFieldDataManager, AllocateFieldDataTwoBuckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  const size_t numRanks = stk::mesh::entity_rank_names().size();
  auto fieldDataManager = std::make_unique<stk::mesh::DefaultFieldDataManager>(numRanks);
  auto * localFieldDataManager = fieldDataManager.get();
  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD, std::move(fieldDataManager));
  stk::mesh::MetaData& meshMetaData = bulkDataPtr->mesh_meta_data();
  initializeTestField<FieldDataType>(meshMetaData);

  const size_t numNodes = 700;
  const size_t extraCapacity = 0;
  testAllocateFieldData<FieldDataType>(*bulkDataPtr, localFieldDataManager, extraCapacity, numNodes);
}

TYPED_TEST(TestDefaultFieldDataManager, TwoEntitiesTwoBuckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  const size_t numRanks = stk::mesh::entity_rank_names().size();
  auto fieldDataManager = std::make_unique<stk::mesh::DefaultFieldDataManager>(numRanks);
  auto * localFieldDataManager = fieldDataManager.get();
  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD, std::move(fieldDataManager));
  stk::mesh::MetaData& meshMetaData = bulkDataPtr->mesh_meta_data();
  createPart(meshMetaData);
  initializeTestField<FieldDataType>(meshMetaData);

  testTwoEntitiesTwoBuckets(*bulkDataPtr, localFieldDataManager);
}

template <typename T>
class TestContiguousFieldDataManager : public testing::Test {};
TYPED_TEST_SUITE(TestContiguousFieldDataManager, TestTypes,);

TYPED_TEST(TestContiguousFieldDataManager, AllocateFieldData)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  auto * localFieldDataManager = fieldDataManager.get();
  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD, std::move(fieldDataManager));
  stk::mesh::MetaData& meshMetaData = bulkDataPtr->mesh_meta_data();
  initializeTestField<FieldDataType>(meshMetaData);
  size_t numNodes = 20;
  const size_t extraCapacity = localFieldDataManager->get_extra_capacity();

  testAllocateFieldData<FieldDataType>(*bulkDataPtr, localFieldDataManager, extraCapacity, numNodes);
}

TYPED_TEST(TestContiguousFieldDataManager, AllocateFieldDataAndReorderBuckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  auto * localFieldDataManager = fieldDataManager.get();
  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD, std::move(fieldDataManager));
  stk::mesh::MetaData& meshMetaData = bulkDataPtr->mesh_meta_data();
  initializeTestField<FieldDataType>(meshMetaData);
  size_t numNodes = 10000;
  const size_t extraCapacity = localFieldDataManager->get_extra_capacity();
  testAllocateFieldData<FieldDataType>(*bulkDataPtr, localFieldDataManager, extraCapacity, numNodes);

  const int num_buckets = static_cast<int>(numNodes/stk::mesh::get_default_maximum_bucket_capacity() + 1);
  std::vector<unsigned> reorderedBucketIds(num_buckets,0);
  for (size_t i=0; i<reorderedBucketIds.size(); i++) {
    reorderedBucketIds[i] = reorderedBucketIds.size()-i-1;
  }
  testReorderBucketFieldData<FieldDataType>(*bulkDataPtr, localFieldDataManager, stk::topology::NODE_RANK,
                                            meshMetaData.get_fields(), reorderedBucketIds);
}

TYPED_TEST(TestContiguousFieldDataManager, TwoEntitiesTwoBuckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  auto * localFieldDataManager = fieldDataManager.get();
  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD, std::move(fieldDataManager));
  stk::mesh::MetaData& meshMetaData = bulkDataPtr->mesh_meta_data();
  createPart(meshMetaData);
  initializeTestField<FieldDataType>(meshMetaData);

  testTwoEntitiesTwoBuckets(*bulkDataPtr, localFieldDataManager);
}

template <typename T>
void initialize2Parts2Fields(stk::mesh::MetaData &meshMetaData)
{
  stk::mesh::Part &part1 = meshMetaData.declare_part("part1");
  stk::mesh::Part &part2 = meshMetaData.declare_part("part2");

  stk::mesh::Field<T> &field1 = meshMetaData.declare_field<T>(stk::topology::NODE_RANK, field_name(1));
  T initial_value1 = 13;
  stk::mesh::put_field_on_mesh( field1, part1, &initial_value1);

  stk::mesh::Field<T> &field2 = meshMetaData.declare_field<T>(stk::topology::NODE_RANK, field_name(2));
  T initial_value2 = 4;
  stk::mesh::put_field_on_mesh( field2, part2, &initial_value2);

  meshMetaData.commit();
}

template <typename T>
void testPartToNodeMapping(const stk::mesh::BulkData &bulkData, const stk::mesh::Part& part,
                           const std::vector<stk::mesh::EntityId> &entityIds, const stk::mesh::Field<T> &goodField,
                           const stk::mesh::Field<T> &nullField)
{
  const stk::mesh::BucketVector &buckets = bulkData.get_buckets(stk::topology::NODE_RANK, part);
  for (size_t i=0; i<buckets.size(); i++) {
    stk::mesh::Bucket &bucket = *buckets[i];
    size_t expectedNumNodes = entityIds.size();
    ASSERT_EQ(expectedNumNodes, bucket.size());
    T *fieldData = stk::mesh::field_data(goodField, bucket);
    T* nullPtr = nullptr;
    ASSERT_NE(nullPtr, fieldData);
    T *nullFieldData = stk::mesh::field_data(nullField, bucket);
    EXPECT_EQ(nullPtr, nullFieldData);
    for (size_t j=0; j<bucket.size(); j++) {
      stk::mesh::Entity node = bucket[j];
      fieldData[j] = bulkData.identifier(node);
      EXPECT_EQ(entityIds[j], bulkData.identifier(node));
    }
  }
}

TYPED_TEST(TestContiguousFieldDataManager, nodalFieldNotOnAllNodeBuckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD, std::move(fieldDataManager));
  stk::mesh::BulkData& bulkData = *bulkDataPtr;
  stk::mesh::MetaData& meshMetaData = bulkData.mesh_meta_data();
  initialize2Parts2Fields<FieldDataType>(meshMetaData);

  bulkData.deactivate_field_updating();

  // =======
  bulkData.modification_begin();

  std::vector<stk::mesh::EntityId> part1Nodes(4,0);

  part1Nodes[0] = 1;
  part1Nodes[1] = 2;
  part1Nodes[2] = 3;
  part1Nodes[3] = 4;

  std::vector<stk::mesh::EntityId> part2Nodes(4,0);

  part2Nodes[0] = 5;
  part2Nodes[1] = 6;
  part2Nodes[2] = 7;
  part2Nodes[3] = 8;

  stk::mesh::Part *part1_ptr = meshMetaData.get_part("part1");
  stk::mesh::Part *part2_ptr = meshMetaData.get_part("part2");

  EXPECT_TRUE(part1_ptr != nullptr);
  EXPECT_TRUE(part2_ptr != nullptr);

  stk::mesh::Part &part1 = *part1_ptr;
  stk::mesh::Part &part2 = *part2_ptr;

  bulkData.declare_node(part1Nodes[0], stk::mesh::ConstPartVector{&part1});
  bulkData.declare_node(part1Nodes[1], stk::mesh::ConstPartVector{&part1});
  bulkData.declare_node(part1Nodes[2], stk::mesh::ConstPartVector{&part1});
  bulkData.declare_node(part1Nodes[3], stk::mesh::ConstPartVector{&part1});

  bulkData.declare_node(part2Nodes[0], stk::mesh::ConstPartVector{&part2});
  bulkData.declare_node(part2Nodes[1], stk::mesh::ConstPartVector{&part2});
  bulkData.declare_node(part2Nodes[2], stk::mesh::ConstPartVector{&part2});
  bulkData.declare_node(part2Nodes[3], stk::mesh::ConstPartVector{&part2});

  bulkData.modification_end();

  bulkData.allocate_field_data();

  stk::mesh::Field<FieldDataType> &field1 = *meshMetaData.get_field<FieldDataType>(stk::topology::NODE_RANK, field_name(1));
  stk::mesh::Field<FieldDataType> &field2 = *meshMetaData.get_field<FieldDataType>(stk::topology::NODE_RANK, field_name(2));
  testPartToNodeMapping(bulkData, part1, part1Nodes, field1, field2);
  testPartToNodeMapping(bulkData, part2, part2Nodes, field2, field1);
}

TYPED_TEST(TestContiguousFieldDataManager, allocate_bucket_field_data)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  auto * localFieldDataManager = fieldDataManager.get();
  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD, std::move(fieldDataManager));
  stk::mesh::BulkData& bulkData = *bulkDataPtr;
  stk::mesh::MetaData& meshMetaData = bulkData.mesh_meta_data();
  initialize2Parts2Fields<FieldDataType>(meshMetaData);

  const stk::mesh::FieldVector &allFields = meshMetaData.get_fields();
  const stk::mesh::FieldVector fields = extract_test_fields(allFields);

  const stk::mesh::FieldBase &fieldOnPart1 = *fields[0];
  const stk::mesh::FieldBase &fieldOnPart2 = *fields[1];
  const stk::mesh::FieldMetaDataVector &part1FieldMetaDataVector = fieldOnPart1.get_meta_data_for_field();
  const stk::mesh::FieldMetaDataVector &part2FieldMetaDataVector = fieldOnPart2.get_meta_data_for_field();
  const size_t expectedNumInitialBucketsInField = 0;
  ASSERT_EQ(expectedNumInitialBucketsInField, part1FieldMetaDataVector.size());


  stk::mesh::PartVector parts(1, meshMetaData.get_part("part1"));
  size_t unusedSize = 0;
  size_t unusedCapacity = 0;
  localFieldDataManager->allocate_bucket_field_data(stk::topology::NODE_RANK, allFields, parts, unusedSize, unusedCapacity);

  size_t expectedNumBucketsInField = 1;
  ASSERT_EQ(expectedNumBucketsInField, part1FieldMetaDataVector.size());
  ASSERT_EQ(expectedNumBucketsInField, part2FieldMetaDataVector.size());
  const unsigned char *nullPointer = nullptr;
  EXPECT_EQ(nullPointer, part1FieldMetaDataVector[0].m_data);
  EXPECT_EQ(nullPointer, part2FieldMetaDataVector[0].m_data);
  int expectedNumBytesPerEntity = sizeof(FieldDataType);
  EXPECT_EQ(expectedNumBytesPerEntity, part1FieldMetaDataVector[0].m_bytesPerEntity);
  expectedNumBytesPerEntity = 0;
  EXPECT_EQ(expectedNumBytesPerEntity, part2FieldMetaDataVector[0].m_bytesPerEntity);


  localFieldDataManager->allocate_bucket_field_data(stk::topology::NODE_RANK, allFields, parts, unusedSize, unusedCapacity);

  expectedNumBucketsInField = 2;
  ASSERT_EQ(expectedNumBucketsInField, part1FieldMetaDataVector.size());
  ASSERT_EQ(expectedNumBucketsInField, part2FieldMetaDataVector.size());
  for (size_t i=0; i<expectedNumBucketsInField; i++) {
    EXPECT_EQ(nullPointer, part1FieldMetaDataVector[i].m_data);
    EXPECT_EQ(nullPointer, part2FieldMetaDataVector[i].m_data);
  }
  expectedNumBytesPerEntity = sizeof(FieldDataType);
  EXPECT_EQ(expectedNumBytesPerEntity, part1FieldMetaDataVector[0].m_bytesPerEntity);
  EXPECT_EQ(expectedNumBytesPerEntity, part1FieldMetaDataVector[1].m_bytesPerEntity);
  expectedNumBytesPerEntity = 0;
  EXPECT_EQ(expectedNumBytesPerEntity, part2FieldMetaDataVector[0].m_bytesPerEntity);
  EXPECT_EQ(expectedNumBytesPerEntity, part2FieldMetaDataVector[1].m_bytesPerEntity);



  parts.push_back(meshMetaData.get_part("part2"));
  localFieldDataManager->allocate_bucket_field_data(stk::topology::NODE_RANK, allFields, parts, unusedSize, unusedCapacity);

  expectedNumBucketsInField = 3;
  ASSERT_EQ(expectedNumBucketsInField, part1FieldMetaDataVector.size());
  ASSERT_EQ(expectedNumBucketsInField, part2FieldMetaDataVector.size());

  for (size_t i=0; i<expectedNumBucketsInField; i++) {
    EXPECT_EQ(nullPointer, part1FieldMetaDataVector[i].m_data);
    EXPECT_EQ(nullPointer, part2FieldMetaDataVector[i].m_data);
  }

  expectedNumBytesPerEntity = sizeof(FieldDataType);
  EXPECT_EQ(expectedNumBytesPerEntity, part1FieldMetaDataVector[0].m_bytesPerEntity);
  EXPECT_EQ(expectedNumBytesPerEntity, part1FieldMetaDataVector[1].m_bytesPerEntity);
  EXPECT_EQ(expectedNumBytesPerEntity, part1FieldMetaDataVector[2].m_bytesPerEntity);

  expectedNumBytesPerEntity = 0;
  EXPECT_EQ(expectedNumBytesPerEntity, part2FieldMetaDataVector[0].m_bytesPerEntity);
  EXPECT_EQ(expectedNumBytesPerEntity, part2FieldMetaDataVector[1].m_bytesPerEntity);
  expectedNumBytesPerEntity = sizeof(FieldDataType);
  EXPECT_EQ(expectedNumBytesPerEntity, part2FieldMetaDataVector[2].m_bytesPerEntity);
}


size_t allocateAndTestNodeBucketFieldData(const std::vector<stk::mesh::PartVector> &partsTable,
                                          const std::vector<std::vector<int> > &bytesPerEntityForField,
                                          stk::mesh::FieldDataManager & fieldDataManager,
                                          const stk::mesh::FieldVector &allFields,
                                          const stk::mesh::FieldVector &fields)
{
  const stk::mesh::FieldBase &fieldOnPart1 = *fields[0];
  const stk::mesh::FieldBase &fieldOnPart2 = *fields[1];
  const stk::mesh::FieldMetaDataVector &part1FieldMetaDataVector = fieldOnPart1.get_meta_data_for_field();
  const stk::mesh::FieldMetaDataVector &part2FieldMetaDataVector = fieldOnPart2.get_meta_data_for_field();

  const size_t bucketSize = 123;
  const size_t bucketCapacity = 123;
  for (size_t i=0; i<partsTable.size(); i++) {
    fieldDataManager.allocate_bucket_field_data(stk::topology::NODE_RANK, allFields, partsTable[i],
                                                bucketSize, bucketCapacity);

    size_t expectedNumBucketsInField = i+1;
    EXPECT_EQ(expectedNumBucketsInField, part1FieldMetaDataVector.size());
    EXPECT_EQ(expectedNumBucketsInField, part2FieldMetaDataVector.size());
    for (size_t j=0; j<expectedNumBucketsInField; j++) {
      int expectedNumBytesPerEntity = bytesPerEntityForField[0][j];
      EXPECT_EQ(expectedNumBytesPerEntity, part1FieldMetaDataVector[j].m_bytesPerEntity);
      expectedNumBytesPerEntity = bytesPerEntityForField[1][j];
      EXPECT_EQ(expectedNumBytesPerEntity, part2FieldMetaDataVector[j].m_bytesPerEntity);
    }
  }
  return bucketCapacity;
}


void deallocateNodeBucketFieldData(const stk::mesh::FieldVector & fields,
                                   stk::mesh::FieldDataManager & fieldDataManager,
                                   size_t numberOfBuckets,
                                   size_t bucketCapacity)
{
  for (size_t bucket_id=0 ; bucket_id<numberOfBuckets ; ++bucket_id) {
    fieldDataManager.deallocate_bucket_field_data(stk::topology::NODE_RANK,bucket_id,bucketCapacity,fields);
  }
}

template <typename T>
void allocate_bucket_field_data_tableBased(stk::mesh::FieldDataManager & fieldDataManager)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const size_t spatialDim = 3;
  stk::mesh::MetaData meshMetaData(spatialDim, stk::mesh::entity_rank_names());
  initialize2Parts2Fields<T>(meshMetaData);


  const stk::mesh::FieldVector &allFields = meshMetaData.get_fields();
  const stk::mesh::FieldVector fields = extract_test_fields(allFields);

  stk::mesh::Part *part1 = meshMetaData.get_part("part1");
  stk::mesh::Part *part2 = meshMetaData.get_part("part2");
  int numBucketsToCreate = 3;
  std::vector<stk::mesh::PartVector> partsTable(numBucketsToCreate);

  partsTable[0].push_back(part1);

  partsTable[1].push_back(part2);

  partsTable[2].push_back(part1);
  partsTable[2].push_back(part2);

  std::vector<std::vector<int> > bytesPerEntityForField(fields.size());
  for (size_t i=0; i<partsTable.size(); i++) {
    int bytesPerEntityForField1 = 0;
    int bytesPerEntityForField2 = 0;
    for (size_t j=0; j<partsTable[i].size(); j++) {
      if (partsTable[i][j] == part1) {
        bytesPerEntityForField1 = sizeof(T);
      }
      if (partsTable[i][j] == part2) {
        bytesPerEntityForField2 = sizeof(T);
      }
    }
    bytesPerEntityForField[0].push_back(bytesPerEntityForField1);
    bytesPerEntityForField[1].push_back(bytesPerEntityForField2);
  }

  size_t bucketCapacity = allocateAndTestNodeBucketFieldData(partsTable, bytesPerEntityForField, fieldDataManager,
                                                             allFields, fields);
  deallocateNodeBucketFieldData(fields,fieldDataManager,partsTable.size(),bucketCapacity);
}

TYPED_TEST(TestContiguousFieldDataManager, allocate_bucket_field_data_tableBased)
{
  using FieldDataType = TypeParam;
  stk::mesh::ContiguousFieldDataManager fieldDataManager;
  allocate_bucket_field_data_tableBased<FieldDataType>(fieldDataManager);
}

TYPED_TEST(TestDefaultFieldDataManager, allocate_bucket_field_data_tableBased)
{
  using FieldDataType = TypeParam;
  const int numRanks = 5;
  stk::mesh::DefaultFieldDataManager fieldDataManager(numRanks);
  allocate_bucket_field_data_tableBased<FieldDataType>(fieldDataManager);
}

template <typename T>
void testAddingSingleEntity(stk::mesh::MetaData &meshMetaData,
                            stk::mesh::ContiguousFieldDataManager & fieldDataManager)
{
  initialize2Parts2Fields<T>(meshMetaData);

  const stk::mesh::FieldVector &allFields = meshMetaData.get_fields();
  const stk::mesh::FieldVector fields = extract_test_fields(allFields);

  const unsigned field1Ordinal = fields[0]->mesh_meta_data_ordinal();
  const unsigned field2Ordinal = fields[1]->mesh_meta_data_ordinal();
  stk::mesh::PartVector parts(1, meshMetaData.get_part("part1"));
  size_t unusedSize = 0;
  size_t unusedCapacity = 0;
  fieldDataManager.allocate_bucket_field_data(stk::topology::NODE_RANK, allFields, parts, unusedSize, unusedCapacity);

  size_t expectedNumBuckets = 1;
  const stk::mesh::FieldMetaDataVector &part1FieldMetaDataVector = fields[field1Ordinal]->get_meta_data_for_field();
  const stk::mesh::FieldMetaDataVector &part2FieldMetaDataVector = fields[field2Ordinal]->get_meta_data_for_field();
  ASSERT_EQ(expectedNumBuckets, part1FieldMetaDataVector.size());
  ASSERT_EQ(expectedNumBuckets, part2FieldMetaDataVector.size());

  int destinationBucketId = 0;
  int destinationBucketOffset = 0;
  fieldDataManager.add_field_data_for_entity(allFields, stk::topology::NODE_RANK, destinationBucketId, destinationBucketOffset);

  const std::vector<size_t> &numBytesAllocated = fieldDataManager.get_num_bytes_allocated_per_field_array();
  ASSERT_EQ(fields.size(), numBytesAllocated.size());
  const unsigned firstBucketIndex = 0;
  const unsigned bytesPerEntity = part1FieldMetaDataVector[firstBucketIndex].m_bytesPerEntity;
  const unsigned alignment = fieldDataManager.get_alignment_bytes();

  std::vector<size_t> expectedNumBytesAllocated {stk::adjust_up_to_alignment_boundary(bytesPerEntity, alignment) +
        fieldDataManager.get_extra_capacity(), 0};
  EXPECT_EQ(expectedNumBytesAllocated[field1Ordinal], numBytesAllocated[field1Ordinal]);
  EXPECT_EQ(expectedNumBytesAllocated[field2Ordinal], numBytesAllocated[field2Ordinal]);

  const std::vector<unsigned char*> &fieldRawData = fieldDataManager.get_field_raw_data();
  ASSERT_EQ(fields.size(), fieldRawData.size());
  EXPECT_EQ(part1FieldMetaDataVector[firstBucketIndex].m_data, fieldRawData[field1Ordinal]);
  unsigned char *nullPointer = nullptr;
  EXPECT_EQ(nullPointer, part2FieldMetaDataVector[firstBucketIndex].m_data);
  EXPECT_EQ(nullPointer, fieldRawData[field2Ordinal]);
}

TYPED_TEST(TestContiguousFieldDataManager, add_field_data_for_entity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  stk::mesh::ContiguousFieldDataManager fieldDataManager;
  const size_t spatialDim = 3;
  stk::mesh::MetaData meshMetaData(spatialDim, stk::mesh::entity_rank_names());

  testAddingSingleEntity<FieldDataType>(meshMetaData, fieldDataManager);

  const stk::mesh::FieldVector &allFields = meshMetaData.get_fields();
  const stk::mesh::FieldVector fields = extract_test_fields(allFields);

  const unsigned field1Ordinal = fields[0]->mesh_meta_data_ordinal();
  const unsigned field2Ordinal = fields[1]->mesh_meta_data_ordinal();

  int destinationBucketId = 0;
  int destinationBucketOffset = 0;

  fieldDataManager.remove_field_data_for_entity(stk::topology::NODE_RANK, destinationBucketId, destinationBucketOffset,
                                                allFields);

  size_t expectedNumBuckets = 1;
  const stk::mesh::FieldMetaDataVector &part1FieldMetaDataVector = fields[field1Ordinal]->get_meta_data_for_field();
  const stk::mesh::FieldMetaDataVector &part2FieldMetaDataVector = fields[field2Ordinal]->get_meta_data_for_field();
  ASSERT_EQ(expectedNumBuckets, part1FieldMetaDataVector.size());

  const std::vector<size_t> &numBytesUsed = fieldDataManager.get_num_bytes_used_per_field_array();
  ASSERT_EQ(fields.size(), numBytesUsed.size());
  std::vector<size_t> expectedNumBytesUsed {0, 0};
  EXPECT_EQ(expectedNumBytesUsed[field1Ordinal], numBytesUsed[field1Ordinal]);
  EXPECT_EQ(expectedNumBytesUsed[field1Ordinal], numBytesUsed[field2Ordinal]);

  const std::vector<size_t> &numBytesAllocated = fieldDataManager.get_num_bytes_allocated_per_field_array();
  ASSERT_EQ(fields.size(), numBytesAllocated.size());
  const unsigned alignment = fieldDataManager.get_alignment_bytes();
  std::vector<size_t> expectedNumBytesAllocated {stk::adjust_up_to_alignment_boundary(sizeof(FieldDataType), alignment) +
        fieldDataManager.get_extra_capacity(), 0};
  EXPECT_EQ(expectedNumBytesAllocated[field1Ordinal], numBytesAllocated[field1Ordinal]);
  EXPECT_EQ(expectedNumBytesAllocated[field2Ordinal], numBytesAllocated[field2Ordinal]);

  const std::vector<unsigned char*> &fieldRawData = fieldDataManager.get_field_raw_data();
  ASSERT_EQ(fields.size(), fieldRawData.size());
  const unsigned firstBucketIndex = 0;
  EXPECT_EQ(part1FieldMetaDataVector[firstBucketIndex].m_data, fieldRawData[field1Ordinal]);
  unsigned char *nullPointer = nullptr;
  EXPECT_EQ(nullPointer, fieldRawData[field2Ordinal]);

  size_t unusedCapacity = 0;
  fieldDataManager.deallocate_bucket_field_data(stk::topology::NODE_RANK, destinationBucketId, unusedCapacity, fields);

  EXPECT_EQ(nullPointer, part1FieldMetaDataVector[firstBucketIndex].m_data);
  EXPECT_EQ(nullPointer, part2FieldMetaDataVector[firstBucketIndex].m_data);
  EXPECT_EQ(0, part1FieldMetaDataVector[firstBucketIndex].m_bytesPerEntity);
  EXPECT_EQ(0, part2FieldMetaDataVector[firstBucketIndex].m_bytesPerEntity);

  std::vector<unsigned> reorderBucketIds;
  fieldDataManager.reorder_bucket_field_data(stk::topology::NODE_RANK, fields, reorderBucketIds);

  expectedNumBytesAllocated = {0 + fieldDataManager.get_extra_capacity(), 0 + fieldDataManager.get_extra_capacity()};
  EXPECT_EQ(expectedNumBytesAllocated[field1Ordinal], numBytesAllocated[field1Ordinal]);
  EXPECT_EQ(expectedNumBytesAllocated[field2Ordinal], numBytesAllocated[field2Ordinal]);

  expectedNumBuckets = 0;
  ASSERT_EQ(expectedNumBuckets, part1FieldMetaDataVector.size());
  ASSERT_EQ(expectedNumBuckets, part2FieldMetaDataVector.size());
  ASSERT_EQ(fields.size(), fieldRawData.size());
  EXPECT_NE(nullPointer, fieldRawData[field1Ordinal]);
  EXPECT_NE(nullPointer, fieldRawData[field2Ordinal]);
}

TYPED_TEST(TestContiguousFieldDataManager, deallocate_nonempty_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  stk::mesh::ContiguousFieldDataManager fieldDataManager;
  const size_t spatialDim = 3;
  stk::mesh::MetaData meshMetaData(spatialDim, stk::mesh::entity_rank_names());

  testAddingSingleEntity<FieldDataType>(meshMetaData, fieldDataManager);

  const stk::mesh::FieldVector &allFields = meshMetaData.get_fields();
  const stk::mesh::FieldVector fields = extract_test_fields(allFields);

  const unsigned field1Ordinal = fields[0]->mesh_meta_data_ordinal();
  const unsigned field2Ordinal = fields[1]->mesh_meta_data_ordinal();

  int destinationBucketId = 0;

  size_t unusedCapacity = 0;
  fieldDataManager.deallocate_bucket_field_data(stk::topology::NODE_RANK, destinationBucketId, unusedCapacity, fields);

  const unsigned firstBucketIndex = 0;
  const stk::mesh::FieldMetaDataVector &part1FieldMetaDataVector = fields[field1Ordinal]->get_meta_data_for_field();
  const stk::mesh::FieldMetaDataVector &part2FieldMetaDataVector = fields[field2Ordinal]->get_meta_data_for_field();
  unsigned char *nullPointer = nullptr;
  EXPECT_EQ(nullPointer, part1FieldMetaDataVector[firstBucketIndex].m_data);
  EXPECT_EQ(nullPointer, part2FieldMetaDataVector[firstBucketIndex].m_data);
  EXPECT_EQ(0, part1FieldMetaDataVector[firstBucketIndex].m_bytesPerEntity);
  EXPECT_EQ(0, part2FieldMetaDataVector[firstBucketIndex].m_bytesPerEntity);

  std::vector<unsigned> reorderBucketIds;
  fieldDataManager.reorder_bucket_field_data(stk::topology::NODE_RANK, fields, reorderBucketIds);

  size_t expectedNumBuckets = 0;
  ASSERT_EQ(expectedNumBuckets, part1FieldMetaDataVector.size());
  ASSERT_EQ(expectedNumBuckets, part2FieldMetaDataVector.size());
  const std::vector<unsigned char*> &fieldRawData = fieldDataManager.get_field_raw_data();
  ASSERT_EQ(fields.size(), fieldRawData.size());
  EXPECT_NE(nullPointer, fieldRawData[field1Ordinal]);
  EXPECT_NE(nullPointer, fieldRawData[field2Ordinal]);
}

template <typename T>
void testForCheating(T *data, size_t sizeOfData, std::vector<T>& valuesErased)
{
  for (size_t i=0; i<valuesErased.size(); i++) {
    ASSERT_FALSE(std::binary_search(data, data+sizeOfData, valuesErased[i]));
  }
}

template <typename T>
void setUpData(std::vector<T>&field, std::vector<size_t>& itemsToErase, std::vector<T> &valuesErased)
{
  size_t numItems=100000;
  field.resize(numItems,0);
  for (size_t i=0; i<field.size(); i++) {
    field[i] = i;
  }

  size_t numItemsToErase=10000;
  itemsToErase.resize(numItemsToErase);
  valuesErased.resize(numItemsToErase,0.0);
  for (size_t i=0; i<itemsToErase.size(); i++) {
    itemsToErase[i] = 10*(numItemsToErase-i-1);
    valuesErased[i] = field[itemsToErase[i]];
  }
}

TEST(ContiguousFieldDataManagerTest, algorithmExploration_eraseOneEntry)
{
  std::vector<size_t> itemsToErase;
  std::vector<double> field;
  std::vector<double> valuesErased;
  setUpData(field, itemsToErase, valuesErased);

  size_t numItems = field.size();
  size_t numItemsToErase = itemsToErase.size();

  std::vector<double>::iterator fbegin = field.begin();
  // Erase
  double startTime = stk::cpu_time();
  for (size_t i=0; i<numItemsToErase; i++) {
    field.erase(fbegin+itemsToErase[i]);
  }
  double totalTime = stk::cpu_time() - startTime;
  EXPECT_EQ(field.size(), numItems-numItemsToErase);
  testForCheating(field.data(), numItems - numItemsToErase, valuesErased);
  std::cerr << "Time = " << totalTime << " s" << std::endl;
}

TEST(ContiguousFieldDataManagerTest, algorithmExploration_fasterEraseOneEntry)
{
  std::vector<size_t> itemsToErase;
  std::vector<double> field;
  std::vector<double> valuesErased;
  setUpData(field, itemsToErase, valuesErased);

  size_t numItems = field.size();
  size_t numItemsToErase = itemsToErase.size();

  double *field_array = field.data();
  size_t field_array_length = field.size();

  // Erase
  double startTime = stk::cpu_time();
  for (size_t i=0; i<numItemsToErase; i++) {
    double* destination = field_array+itemsToErase[i];
    std::memmove(destination, destination+1, (field_array_length-itemsToErase[i]-1)*sizeof(double));
    --field_array_length;
  }

  double totalTime = stk::cpu_time() - startTime;
  EXPECT_EQ(field_array_length, numItems-numItemsToErase);
  testForCheating(field.data(), numItems - numItemsToErase, valuesErased);
  std::cerr << "Time = " << totalTime << " s" << std::endl;
}

TEST(ContiguousFieldDataManagerTest, algorithmExploration_slowestEraseOneEntry)
{
  std::vector<size_t> itemsToErase;
  std::vector<double> field;
  std::vector<double> valuesErased;
  setUpData(field, itemsToErase, valuesErased);

  size_t numItems = field.size();
  size_t numItemsToErase = itemsToErase.size();

  double *field_array = field.data();
  size_t field_array_length = field.size();

  // Erase
  double startTime = stk::cpu_time();
  for (size_t i=0; i<numItemsToErase; i++) {
    double* destination = field_array+itemsToErase[i];
    double* source = destination+1;
    int numItemsToSlide = field_array_length-itemsToErase[i]-1;
    for (int item=0; item<numItemsToSlide; ++item) {
      *destination++ = *source++;
    }
    --field_array_length;
  }

  double totalTime = stk::cpu_time() - startTime;
  EXPECT_EQ(field_array_length, numItems-numItemsToErase);
  testForCheating(field.data(), numItems - numItemsToErase, valuesErased);
  std::cerr << "Time = " << totalTime << " s" << std::endl;
}

TEST(ContiguousFieldDataManagerTest, algorithmExploration_memcpyishEraseOneEntry)
{
  std::vector<size_t> itemsToErase;
  std::vector<double> field;
  std::vector<double> valuesErased;
  setUpData(field, itemsToErase, valuesErased);

  size_t numItems = field.size();
  size_t numItemsToErase = itemsToErase.size();

  std::vector<double> scratchField(numItems,0);
  double *field_array = field.data();
  size_t field_array_length = field.size();
  double *scratchData = scratchField.data();

  // Erase
  double startTime = stk::cpu_time();
  for (size_t i=0; i<numItemsToErase; i++) {
    size_t numBytesToCopy = (field_array_length-itemsToErase[i]-1)*sizeof(double);
    double* destination = field_array+itemsToErase[i];
    std::memcpy(scratchData, destination+1, numBytesToCopy);
    std::memcpy(destination, scratchData, numBytesToCopy);
    --field_array_length;
  }

  double totalTime = stk::cpu_time() - startTime;
  EXPECT_EQ(field_array_length, numItems-numItemsToErase);
  testForCheating(field.data(), numItems - numItemsToErase, valuesErased);
  std::cerr << "Time = " << totalTime << " s" << std::endl;
}

TEST(ContiguousFieldDataManagerTest, algorithmExploration_batchDeletion)
{
  std::vector<size_t> itemsToErase;
  std::vector<double> field;
  std::vector<double> valuesErased;
  setUpData(field, itemsToErase, valuesErased);

  size_t numItems = field.size();
  size_t numItemsToErase = itemsToErase.size();
  std::vector<double> scratchField(field.size(),0);

  double startTime = stk::cpu_time();
  std::sort(itemsToErase.begin(), itemsToErase.end());

  std::vector<double *> startingPtr(numItemsToErase+1,nullptr);
  std::vector<int> itemsThisChunk(numItemsToErase+1, 0);
  std::vector<int> distanceToSlideLeft(numItemsToErase+1,0);

  startingPtr[0] = field.data();
  itemsThisChunk[0] = itemsToErase[0];
  distanceToSlideLeft[0] = 0;

  for (size_t i=1; i<itemsToErase.size(); i++) {
    startingPtr[i] = &field[itemsToErase[i-1]+1];
    itemsThisChunk[i] = itemsToErase[i]-itemsToErase[i-1]-1;
    distanceToSlideLeft[i] = i;
  }

  int lastIndex = itemsToErase.size();
  startingPtr[lastIndex] = &field[itemsToErase[lastIndex-1]] + 1;
  itemsThisChunk[lastIndex] = field.size()-itemsToErase[lastIndex-1]-1;
  distanceToSlideLeft[lastIndex] = lastIndex;

  size_t field_array_length = field.size();

  for (size_t i=1; i<numItemsToErase+1; i++) {
    double* destination = startingPtr[i]-distanceToSlideLeft[i];
    std::memmove(destination, startingPtr[i], (itemsThisChunk[i])*sizeof(double));
    field_array_length--;
  }

  double totalTime = stk::cpu_time() - startTime;
  EXPECT_EQ(field_array_length, numItems-numItemsToErase);
  testForCheating(field.data(), numItems - numItemsToErase, valuesErased);
  std::cerr << "Time = " << totalTime << " s" << std::endl;
}

}


