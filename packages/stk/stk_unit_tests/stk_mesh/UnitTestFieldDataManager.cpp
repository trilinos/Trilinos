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
void testAllocateFieldData(stk::mesh::BulkData& bulkData, const size_t extraCapacityInBytes, const size_t numNodes)
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

  const stk::mesh::FieldDataManager& fieldDataManager = bulkData.get_field_data_manager();

  const stk::mesh::FieldVector &fields = meshMetaData.get_fields();
  for (size_t field_index = 0; field_index < fields.size(); field_index++) {
    const stk::mesh::FieldBase &field = *fields[field_index];
    if (is_test_field(field)) {
      const T *initial_value = reinterpret_cast<const T*>(field.get_initial_value());

      size_t totalBytesAllocatedForField = 0;
      for (const stk::mesh::Bucket * bucket : bulkData.buckets(stk::topology::NODE_RANK)) {
        size_t bytesPerEntity = field.get_meta_data_for_field()[bucket->bucket_id()].m_bytesPerEntity;
        size_t numEntitiesAllocated = bucket->capacity();
        totalBytesAllocatedForField += stk::adjust_up_to_alignment_boundary(numEntitiesAllocated*bytesPerEntity,
                                                                            fieldDataManager.get_alignment_bytes());
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
                fieldDataManager.get_num_bytes_allocated_on_field(field.mesh_meta_data_ordinal()));
    }
  }
}

template <typename T>
void testReorderBucketFieldData(stk::mesh::BulkData& bulkData, stk::mesh::EntityRank rank,
                                const stk::mesh::FieldVector& fields, const std::vector<unsigned> &reorderedBucketIds)
{
  const stk::mesh::BucketVector& buckets = bulkData.buckets(rank);
  for (size_t b=0; b<buckets.size(); ++b) {
    const stk::mesh::Bucket& bucket = *buckets[b];
    for (size_t i=0; i<fields.size(); ++i) {
      if (fields[i]->data_traits().size_of == sizeof(T)) {
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
  }

  stk::mesh::FieldDataManager& fieldDataManager = bulkData.get_field_data_manager();
  fieldDataManager.reorder_bucket_field_data(rank, fields, reorderedBucketIds);

  for (size_t bucketIndex=0; bucketIndex<buckets.size(); ++bucketIndex) {
    size_t oldBucketIndex = reorderedBucketIds[bucketIndex];
    size_t oldBucketSize = buckets[oldBucketIndex]->size();

    for (size_t i=0; i<fields.size(); ++i) {
      if (fields[i]->data_traits().size_of == sizeof(T)) {
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
}

void testTwoEntitiesTwoBuckets(stk::mesh::BulkData &bulkData)
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

std::shared_ptr<stk::mesh::BulkData> build_mesh(unsigned spatialDim, stk::ParallelMachine comm)
{
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(spatialDim);
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

  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::MetaData& meshMetaData = bulkDataPtr->mesh_meta_data();
  initializeTestField<FieldDataType>(meshMetaData);

  size_t numNodes = 20;
  size_t extraCapacity = 0;
  testAllocateFieldData<FieldDataType>(*bulkDataPtr, extraCapacity, numNodes);
}

TYPED_TEST(TestDefaultFieldDataManager, AllocateFieldDataTwoBuckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::MetaData& meshMetaData = bulkDataPtr->mesh_meta_data();
  initializeTestField<FieldDataType>(meshMetaData);

  const size_t numNodes = 700;
  const size_t extraCapacity = 0;
  testAllocateFieldData<FieldDataType>(*bulkDataPtr, extraCapacity, numNodes);
}

TYPED_TEST(TestDefaultFieldDataManager, TwoEntitiesTwoBuckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  using FieldDataType = TypeParam;
  const size_t spatialDim = 3;

  std::shared_ptr<stk::mesh::BulkData> bulkDataPtr = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::MetaData& meshMetaData = bulkDataPtr->mesh_meta_data();
  createPart(meshMetaData);
  initializeTestField<FieldDataType>(meshMetaData);

  testTwoEntitiesTwoBuckets(*bulkDataPtr);
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

size_t allocateAndTestNodeBucketFieldData(const std::vector<stk::mesh::PartVector> &partsTable,
                                          const std::vector<std::vector<int> > &bytesPerEntityForField,
                                          const stk::mesh::FieldVector &allFields,
                                          const stk::mesh::FieldVector &fields)
{
  const stk::mesh::FieldBase &fieldOnPart1 = *fields[0];
  const stk::mesh::FieldBase &fieldOnPart2 = *fields[1];
  const stk::mesh::FieldMetaDataArrayType &part1FieldMetaDataVector = fieldOnPart1.get_meta_data_for_field();
  const stk::mesh::FieldMetaDataArrayType &part2FieldMetaDataVector = fieldOnPart2.get_meta_data_for_field();

  stk::mesh::FieldDataManager& fieldDataManager = fieldOnPart1.get_mesh().get_field_data_manager();

  const size_t bucketSize = 123;
  const size_t bucketCapacity = 123;
  for (size_t i=0; i<partsTable.size(); i++) {
    fieldDataManager.allocate_bucket_field_data(stk::topology::NODE_RANK, allFields, partsTable[i],
                                                bucketSize, bucketCapacity);

    size_t expectedNumBucketsInField = i+1;
    EXPECT_EQ(expectedNumBucketsInField, part1FieldMetaDataVector.extent(0));
    EXPECT_EQ(expectedNumBucketsInField, part2FieldMetaDataVector.extent(0));
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
                                   size_t numberOfBuckets,
                                   size_t bucketCapacity)
{
  stk::mesh::FieldDataManager& fieldDataManager = fields.front()->get_mesh().get_field_data_manager();

  for (size_t bucket_id=0 ; bucket_id<numberOfBuckets ; ++bucket_id) {
    fieldDataManager.deallocate_bucket_field_data(stk::topology::NODE_RANK, bucket_id, bucketCapacity, fields);
  }
}

template <typename T>
void allocate_bucket_field_data_tableBased()
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  const size_t spatialDim = 3;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  auto bulkData = builder.set_spatial_dimension(spatialDim).create();
  stk::mesh::MetaData& meshMetaData = bulkData->mesh_meta_data();

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

  size_t bucketCapacity = allocateAndTestNodeBucketFieldData(partsTable, bytesPerEntityForField, allFields, fields);
  deallocateNodeBucketFieldData(fields, partsTable.size(), bucketCapacity);
}

TYPED_TEST(TestDefaultFieldDataManager, allocate_bucket_field_data_tableBased)
{
  using FieldDataType = TypeParam;
  allocate_bucket_field_data_tableBased<FieldDataType>();
}

}


