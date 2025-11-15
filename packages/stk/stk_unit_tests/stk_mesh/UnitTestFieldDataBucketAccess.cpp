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
#include "FieldDataAccessFixture.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/Selector.hpp"   // for operator<<, Selector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/ConstFieldData.hpp>
#include <stk_mesh/base/FieldDataBase.hpp>
#include <stk_mesh/base/FieldIndexTypes.hpp>

namespace {

class FieldDataBucketAccess : public FieldDataAccessFixture {};

//==============================================================================
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_scalar(const stk::mesh::BulkData& bulk,
                      const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Bucket indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        bucketValues(entity) = ++value*10;
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        EXPECT_EQ(constBucketValues(entity), ++value*10);
      }
    }
  }

  // BucketId indexing (with Entity indexing from BucketValues)
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : bucketValues.entities()) {
        bucketValues(entity) = ++value*20;
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : constBucketValues.entities()) {
        EXPECT_EQ(constBucketValues(entity), ++value*20);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_scalar(get_bulk(),
                   field.data<stk::mesh::ReadWrite>(),
                   field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite>(),
                   fieldBase.data<int>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_field_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_scalar(get_bulk(),
                   field.data<stk::mesh::ReadWrite>(),
                   field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_fieldBase_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_field_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_scalar(get_bulk(),
                   field.data<stk::mesh::ReadWrite>(),
                   field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_fieldBase_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_fieldBase_layoutLeft_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_fieldBase_layoutAuto_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_fieldBase_layoutRight_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_fieldBase_layoutAuto_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_scalar_pointer(const stk::mesh::BulkData& bulk,
                              const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        bucketValues(entity) = ++value*10;
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        EXPECT_EQ(constBucketPtr[entity*entityStride], ++value*10);
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      int* bucketPtr = bucketValues.pointer();
      const int entityStride = bucketValues.entity_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        bucketPtr[entity*entityStride] = ++value*20;
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        EXPECT_EQ(constBucketValues(entity), ++value*20);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_scalar_pointer(get_bulk(),
                           field.data<stk::mesh::ReadWrite>(),
                           field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_scalar_pointer(get_bulk(),
                           field.data<stk::mesh::ReadWrite>(),
                           field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_scalar_pointer_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_scalar_pointer(get_bulk(),
                           field.data<stk::mesh::ReadWrite>(),
                           field.data());
}


//==============================================================================
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_multi_component(const stk::mesh::BulkData& bulk,
                               const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Bucket indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ComponentIdx component : bucketValues.components()) {
          bucketValues(entity, component) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
          EXPECT_EQ(constBucketValues(entity, component), ++value*10);
        }
      }
    }
  }

  // BucketId indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ComponentIdx component : bucketValues.components()) {
          bucketValues(entity, component) = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
          EXPECT_EQ(constBucketValues(entity, component), ++value*20);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_component(get_bulk(),
                            field.data<stk::mesh::ReadWrite>(),
                            field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite>(),
                            fieldBase.data<int>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_field_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_component(get_bulk(),
                            field.data<stk::mesh::ReadWrite>(),
                            field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_fieldBase_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_field_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_component(get_bulk(),
                            field.data<stk::mesh::ReadWrite>(),
                            field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_fieldBase_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_fieldBase_layoutLeft_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_fieldBase_layoutAuto_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_fieldBase_layoutRight_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_fieldBase_layoutAuto_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_multi_component_pointer(const stk::mesh::BulkData& bulk,
                                       const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ComponentIdx component : bucketValues.components()) {
          bucketValues(entity, component) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      const int componentStride = constBucketValues.component_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        for (int component = 0; component < constBucketValues.num_components(); ++component) {
          EXPECT_EQ(constBucketPtr[entity*entityStride + component*componentStride], ++value*10);
        }
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      int* bucketPtr = bucketValues.pointer();
      const int entityStride = bucketValues.entity_stride();
      const int componentStride = bucketValues.component_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        for (int component = 0; component < bucketValues.num_components(); ++component) {
          bucketPtr[entity*entityStride + component*componentStride] = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
          EXPECT_EQ(constBucketValues(entity, component), ++value*20);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_component_pointer(get_bulk(),
                                    field.data<stk::mesh::ReadWrite>(),
                                    field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_component_pointer(get_bulk(),
                                    field.data<stk::mesh::ReadWrite>(),
                                    field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponent_pointer_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_component_pointer(get_bulk(),
                                    field.data<stk::mesh::ReadWrite>(),
                                    field.data());
}


//==============================================================================
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_multi_copy(const stk::mesh::BulkData& bulk,
                          const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Bucket indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
          bucketValues(entity, copy) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
          EXPECT_EQ(constBucketValues(entity, copy), ++value*10);
        }
      }
    }
  }

  // BucketId indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
          bucketValues(entity, copy) = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
          EXPECT_EQ(constBucketValues(entity, copy), ++value*20);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_copy(get_bulk(),
                       field.data<stk::mesh::ReadWrite>(),
                       field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite>(),
                       fieldBase.data<int>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_field_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_copy(get_bulk(),
                       field.data<stk::mesh::ReadWrite>(),
                       field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_fieldBase_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_field_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_copy(get_bulk(),
                       field.data<stk::mesh::ReadWrite>(),
                       field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_fieldBase_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_fieldBase_layoutLeft_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_fieldBase_layoutAuto_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_fieldBase_layoutRight_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_fieldBase_layoutAuto_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}


//------------------------------------------------------------------------------
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_multi_copy_pointer(const stk::mesh::BulkData& bulk,
                                  const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
          bucketValues(entity, copy) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      const int copyStride = constBucketValues.copy_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        for (int copy = 0; copy < constBucketValues.num_copies(); ++copy) {
          EXPECT_EQ(constBucketPtr[entity*entityStride + copy*copyStride], ++value*10);
        }
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      int* bucketPtr = bucketValues.pointer();
      const int entityStride = bucketValues.entity_stride();
      const int copyStride = bucketValues.copy_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        for (int copy = 0; copy < bucketValues.num_copies(); ++copy) {
          bucketPtr[entity*entityStride + copy*copyStride] = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
          EXPECT_EQ(constBucketValues(entity, copy), ++value*20);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_copy_pointer(get_bulk(),
                               field.data<stk::mesh::ReadWrite>(),
                               field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_copy_pointer(get_bulk(),
                               field.data<stk::mesh::ReadWrite>(),
                               field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_pointer_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_copy_pointer(get_bulk(),
                               field.data<stk::mesh::ReadWrite>(),
                               field.data());
}


//==============================================================================
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_multi_copy_multi_component(const stk::mesh::BulkData& bulk,
                                          const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Bucket indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
          for (stk::mesh::ComponentIdx component : bucketValues.components()) {
            bucketValues(entity, copy, component) = ++value*10;
          }
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
          for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
            EXPECT_EQ(constBucketValues(entity, copy, component), ++value*10);
          }
        }
      }
    }
  }

  // BucketId indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
          for (stk::mesh::ComponentIdx component : bucketValues.components()) {
            bucketValues(entity, copy, component) = ++value*20;
          }
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
          for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
            EXPECT_EQ(constBucketValues(entity, copy, component), ++value*20);
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_copy_multi_component(
        get_bulk(),
        field.data<stk::mesh::ReadWrite>(),
        field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_host_multi_copy_multi_component(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite>(),
        fieldBase.data<int>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_field_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_copy_multi_component(
        get_bulk(),
        field.data<stk::mesh::ReadWrite>(),
        field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_fieldBase_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_copy_multi_component(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_field_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_copy_multi_component(
        get_bulk(),
        field.data<stk::mesh::ReadWrite>(),
        field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_fieldBase_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_copy_multi_component(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_fieldBase_layoutLeft_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_copy_multi_component(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_fieldBase_layoutAuto_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_copy_multi_component(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_fieldBase_layoutRight_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_copy_multi_component(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_fieldBase_layoutAuto_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_copy_multi_component(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  auto fieldData = field.data<stk::mesh::ReadWrite>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto bucketValues = fieldData.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity(0); entity < bucket->num_entities(); ++entity) {
      for (stk::mesh::CopyIdx copy(0); copy < bucketValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < bucketValues.num_components(); ++component) {
          bucketValues(entity, copy, component) = ++value;
        }
      }
    }
  }

  value = 0;
  auto constFieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto constBucketValues = constFieldData.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity(0); entity < bucket->num_entities(); ++entity) {
      for (stk::mesh::CopyIdx copy(0); copy < constBucketValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constBucketValues.num_components(); ++component) {
        EXPECT_EQ(constBucketValues(entity, copy, component), ++value);
        }
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto bucketValuesBase = fieldDataBase.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity(0); entity < bucket->num_entities(); ++entity) {
      for (stk::mesh::CopyIdx copy(0); copy < bucketValuesBase.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < bucketValuesBase.num_components(); ++component) {
        bucketValuesBase(entity, copy, component) = ++value*10;
        }
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto constBucketValuesBase = constFieldDataBase.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity(0); entity < bucket->num_entities(); ++entity) {
      for (stk::mesh::CopyIdx copy(0); copy < constBucketValuesBase.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constBucketValuesBase.num_components(); ++component) {
          EXPECT_EQ(constBucketValuesBase(entity, copy, component), ++value*10);
        }
      }
    }
  }
}


//------------------------------------------------------------------------------
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_multi_copy_multi_component_pointer(const stk::mesh::BulkData& bulk,
                                                  const FieldDataType& fieldData,
                                                  const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
          for (stk::mesh::ComponentIdx component : bucketValues.components()) {
            bucketValues(entity, copy, component) = ++value*10;
          }
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      const int copyStride = constBucketValues.copy_stride();
      const int componentStride = constBucketValues.component_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        for (int copy = 0; copy < constBucketValues.num_copies(); ++copy) {
          for (int component = 0; component < constBucketValues.num_components(); ++component) {
            EXPECT_EQ(constBucketPtr[entity*entityStride + copy*copyStride + component*componentStride], ++value*10);
          }
        }
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      int* bucketPtr = bucketValues.pointer();
      const int entityStride = bucketValues.entity_stride();
      const int copyStride = bucketValues.copy_stride();
      const int componentStride = bucketValues.component_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        for (int copy = 0; copy < bucketValues.num_copies(); ++copy) {
          for (int component = 0; component < bucketValues.num_components(); ++component) {
            bucketPtr[entity*entityStride + copy*copyStride + component*componentStride] = ++value*20;
          }
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
          for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
            EXPECT_EQ(constBucketValues(entity, copy, component), ++value*20);
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_copy_multi_component_pointer(get_bulk(),
                                               field.data<stk::mesh::ReadWrite>(),
                                               field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_copy_multi_component_pointer(get_bulk(),
                                               field.data<stk::mesh::ReadWrite>(),
                                               field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopy_multiComponent_pointer_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_copy_multi_component_pointer(get_bulk(),
                                               field.data<stk::mesh::ReadWrite>(),
                                               field.data());
}


//==============================================================================
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_multi_scalar(const stk::mesh::BulkData& bulk,
                            const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Bucket indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
          bucketValues(entity, scalar) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : constBucketValues.scalars()) {
          EXPECT_EQ(constBucketValues(entity, scalar), ++value*10);
        }
      }
    }
  }

  // BucketId indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
          bucketValues(entity, scalar) = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(bucket->bucket_id());
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : constBucketValues.scalars()) {
          EXPECT_EQ(constBucketValues(entity, scalar), ++value*20);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_scalar(
        get_bulk(),
        field.data<stk::mesh::ReadWrite>(),
        field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_host_multi_scalar(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite>(),
        fieldBase.data<int>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_field_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_scalar(
        get_bulk(),
        field.data<stk::mesh::ReadWrite>(),
        field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_fieldBase_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_scalar(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_field_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_scalar(
        get_bulk(),
        field.data<stk::mesh::ReadWrite>(),
        field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_fieldBase_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_scalar(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_fieldBase_layoutLeft_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_scalar(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_fieldBase_layoutAuto_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_scalar(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_fieldBase_layoutRight_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_scalar(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_fieldBase_layoutAuto_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_scalar(
        get_bulk(),
        fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
        fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  auto fieldData = field.data<stk::mesh::ReadWrite>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto bucketValues = fieldData.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity(0); entity < bucket->num_entities(); ++entity) {
      for (stk::mesh::ScalarIdx scalar(0); scalar < bucketValues.num_scalars(); ++scalar) {
        bucketValues(entity, scalar) = ++value;
      }
    }
  }

  value = 0;
  auto constFieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto constBucketValues = constFieldData.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity(0); entity < bucket->num_entities(); ++entity) {
      for (stk::mesh::ScalarIdx scalar(0); scalar < constBucketValues.num_scalars(); ++scalar) {
        EXPECT_EQ(constBucketValues(entity, scalar), ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto bucketValuesBase = fieldDataBase.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity(0); entity < bucket->num_entities(); ++entity) {
      for (stk::mesh::ScalarIdx scalar(0); scalar < bucketValuesBase.num_scalars(); ++scalar) {
        bucketValuesBase(entity, scalar) = ++value*10;
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int>();
  for (stk::mesh::Bucket* bucket : buckets) {
    auto constBucketValuesBase = constFieldDataBase.bucket_values(*bucket);
    for (stk::mesh::EntityIdx entity(0); entity < bucket->num_entities(); ++entity) {
      for (stk::mesh::ScalarIdx scalar(0); scalar < constBucketValuesBase.num_scalars(); ++scalar) {
        EXPECT_EQ(constBucketValuesBase(entity, scalar), ++value*10);
      }
    }
  }
}

//------------------------------------------------------------------------------
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_multi_scalar_pointer(const stk::mesh::BulkData& bulk,
                                    const FieldDataType& fieldData,
                                    const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
          bucketValues(entity, scalar) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      const int scalarStride = constBucketValues.scalar_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        for (int scalar = 0; scalar < constBucketValues.num_scalars(); ++scalar) {
          EXPECT_EQ(constBucketPtr[entity*entityStride + scalar*scalarStride], ++value*10);
        }
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      int* bucketPtr = bucketValues.pointer();
      const int entityStride = bucketValues.entity_stride();
      const int scalarStride = bucketValues.scalar_stride();
      for (int entity = 0; entity < bucket->num_entities(); ++entity) {
        for (int scalar = 0; scalar < bucketValues.num_scalars(); ++scalar) {
          bucketPtr[entity*entityStride + scalar*scalarStride] = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : constBucketValues.scalars()) {
          EXPECT_EQ(constBucketValues(entity, scalar), ++value*20);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_scalar_pointer(get_bulk(),
                                 field.data<stk::mesh::ReadWrite>(),
                                 field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_scalar_pointer(get_bulk(),
                                 field.data<stk::mesh::ReadWrite>(),
                                 field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalar_pointer_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_scalar_pointer(get_bulk(),
                                 field.data<stk::mesh::ReadWrite>(),
                                 field.data());
}


//==============================================================================
template <typename FieldType, typename FieldDataType, typename ConstFieldDataType>
void test_device_scalar(const stk::mesh::BulkData& bulk, FieldType& field,
                        const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      stk::mesh::EntityIdx numEntities = bucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          bucketValues(entity) = bucketId*10 + entity();
        }
      );
    }
  );

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          NGP_EXPECT_EQ(constBucketValues(entity), bucketId*10 + entity());
        }
      );
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_scalar_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_scalar(get_bulk(), field,
                     field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                     field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_scalar_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_device_scalar(get_bulk(), fieldBase,
                     fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                     fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
void test_device_scalar_pointer(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      stk::mesh::EntityIdx numEntities = bucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          bucketValues(entity) = bucketId*10 + entity();
        }
      );
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          NGP_EXPECT_EQ(constBucketPtr[entity*entityStride], bucketId*10 + entity());
        }
      );
    }
  );
}

NGP_TEST_F(FieldDataBucketAccess, device_scalar_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();
  test_device_scalar_pointer(get_bulk(), *m_field);
}


//==============================================================================
template <typename FieldType, typename FieldDataType, typename ConstFieldDataType>
void test_device_multi_component(const stk::mesh::BulkData& bulk, FieldType& field,
                                 const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      stk::mesh::EntityIdx numEntities = bucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ComponentIdx component : bucketValues.components()) {
            bucketValues(entity, component) = bucketId*100 + entity()*10 + component();
          }
        }
      );
    }
  );

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
            NGP_EXPECT_EQ(constBucketValues(entity, component), (bucketId*100 + entity()*10 + component()));
          }
        }
      );
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiComponent_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_component(get_bulk(), field,
                              field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                              field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiComponent_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_device_multi_component(get_bulk(), fieldBase,
                              fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                              fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiComponent_field_async)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_component(get_bulk(), field,
                              field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(stk::ngp::ExecSpace()),
                              field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(stk::ngp::ExecSpace()));
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiComponent_fieldBase_async)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_device_multi_component(get_bulk(), fieldBase,
                              fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(stk::ngp::ExecSpace()),
                              fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(stk::ngp::ExecSpace()));
}

//------------------------------------------------------------------------------
void test_device_multi_component_pointer(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      stk::mesh::EntityIdx numEntities = bucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ComponentIdx component : bucketValues.components()) {
            bucketValues(entity, component) = bucketId*100 + entity()*10 + component();
          }
        }
      );
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      const int componentStride = constBucketValues.component_stride();
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (int component = 0; component < constBucketValues.num_components(); ++component) {
            NGP_EXPECT_EQ(constBucketPtr[entity*entityStride + component*componentStride],
                          (bucketId*100 + entity()*10 + component));
          }
        }
      );
    }
  );
}

NGP_TEST_F(FieldDataBucketAccess, device_multiComponent_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();
  test_device_multi_component_pointer(get_bulk(), *m_field);
}


//==============================================================================
template <typename FieldType, typename FieldDataType, typename ConstFieldDataType>
void test_device_multi_copy(const stk::mesh::BulkData& bulk, FieldType& field,
                            const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      stk::mesh::EntityIdx numEntities = bucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
            bucketValues(entity, copy) = bucketId*100 + entity()*10 + copy();
          }
        }
      );
    }
  );

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
            NGP_EXPECT_EQ(constBucketValues(entity, copy), (bucketId*100 + entity()*10 + copy()));
          }
        }
      );
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiCopy_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_copy(get_bulk(), field,
                         field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                         field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiCopy_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_device_multi_copy(get_bulk(), fieldBase,
                         fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                         fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
void test_device_multi_copy_pointer(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  // Write and read values from Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      stk::mesh::EntityIdx numEntities = bucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
            bucketValues(entity, copy) = bucketId*100 + entity()*10 + copy();
          }
        }
      );
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      const int copyStride = constBucketValues.copy_stride();
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (int copy = 0; copy < constBucketValues.num_copies(); ++copy) {
            NGP_EXPECT_EQ(constBucketPtr[entity*entityStride + copy*copyStride], (bucketId*100 + entity()*10 + copy));
          }
        }
      );
    }
  );
}

NGP_TEST_F(FieldDataBucketAccess, device_multiCopy_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();
  test_device_multi_copy_pointer(get_bulk(), *m_field);
}


//==============================================================================
template <typename FieldType, typename FieldDataType, typename ConstFieldDataType>
void test_device_multi_copy_multi_component(const stk::mesh::BulkData& bulk, FieldType& field,
                                            const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = bucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
            for (stk::mesh::ComponentIdx component : bucketValues.components()) {
              bucketValues(entity, copy, component) = bucketId*1000 + entity()*100 + copy()*10 + component();
            }
          }
        }
      );
    }
  );

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
            for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
              NGP_EXPECT_EQ(constBucketValues(entity, copy, component),
                            (bucketId*1000 + entity()*100 + copy()*10 + component()));
            }
          }
        }
      );
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiCopy_multiComponent_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_copy_multi_component(get_bulk(), field,
                                         field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                                         field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiCopy_multiComponent_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_device_multi_copy_multi_component(get_bulk(), fieldBase,
                                         fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                                         fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
void test_device_multi_copy_multi_component_pointer(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = bucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
            for (stk::mesh::ComponentIdx component : bucketValues.components()) {
              bucketValues(entity, copy, component) = bucketId*1000 + entity()*100 + copy()*10 + component();
            }
          }
        }
      );
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      const int copyStride = constBucketValues.copy_stride();
      const int componentStride = constBucketValues.component_stride();
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (int copy = 0; copy < constBucketValues.num_copies(); ++copy) {
            for (int component = 0; component < constBucketValues.num_components(); ++component) {
              NGP_EXPECT_EQ(constBucketPtr[entity*entityStride + copy*copyStride + component*componentStride],
                            (bucketId*1000 + entity()*100 + copy*10 + component));
            }
          }
        }
      );
    }
  );
}

NGP_TEST_F(FieldDataBucketAccess, device_multiCopy_multiComponent_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_multi_copy_multi_component_pointer(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_multi_copy_multi_component_traditional_for_loop(stk::mesh::BulkData& bulk,
                                                                 stk::mesh::Field<int>& field)
{
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(field);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = bucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy(0); copy < bucketValues.num_copies(); ++copy) {
            for (stk::mesh::ComponentIdx component(0); component < bucketValues.num_components(); ++component) {
              bucketValues(entity, copy, component) = bucketId*1000 + entity()*100 + copy()*10 + component();
            }
          }
        }
      );
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy(0); copy < constBucketValues.num_copies(); ++copy) {
            for (stk::mesh::ComponentIdx component(0); component < constBucketValues.num_components(); ++component) {
              NGP_EXPECT_EQ(constBucketValues(entity, copy, component),
                            (bucketId*1000 + entity()*100 + copy()*10 + component()));
            }
          }
        }
      );
    }
  );

  // Write values to FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValuesBase = fieldDataBase.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = bucketValuesBase.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy(0); copy < bucketValuesBase.num_copies(); ++copy) {
            for (stk::mesh::ComponentIdx component(0); component < bucketValuesBase.num_components(); ++component) {
              bucketValuesBase(entity, copy, component) = bucketId*10000 + entity()*1000 + copy()*100 + component()*10;
            }
          }
        }
      );
    }
  );

  // Read const values from FieldBase
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValuesBase = constFieldDataBase.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = constBucketValuesBase.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::CopyIdx copy(0); copy < constBucketValuesBase.num_copies(); ++copy) {
            for (stk::mesh::ComponentIdx component(0); component < constBucketValuesBase.num_components(); ++component) {
              NGP_EXPECT_EQ(constBucketValuesBase(entity, copy, component),
                            (bucketId*10000 + entity()*1000 + copy()*100 + component()*10));
            }
          }
        }
      );
    }
  );
}

NGP_TEST_F(FieldDataBucketAccess, device_multiCopy_multiComponent_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_multi_copy_multi_component_traditional_for_loop(get_bulk(), *m_field);
}


//==============================================================================
template <typename FieldType>
void test_host_to_device_scalar(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        bucketValues(entity) = bucket->bucket_id()*10 + entity();
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto constBucketValues = constFieldData.bucket_values(bucketId);
        const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            NGP_EXPECT_EQ(constBucketValues(entity), (bucketId*10 + entity()));
          }
        );
      }
    );
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();
  test_host_to_device_scalar(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_scalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_left_field();
  test_host_to_device_scalar(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_scalar_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_right_field();
  test_host_to_device_scalar(get_bulk(), *m_rightField);
}

//==============================================================================
template <typename FieldType, typename FieldDataType, typename ConstFieldDataType>
void test_device_multi_scalar(const stk::mesh::BulkData& bulk, FieldType& field,
                              const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = bucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
            bucketValues(entity, scalar) = bucketId*1000 + entity()*100 + scalar()*10;
          }
        }
      );
    }
  );

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ScalarIdx scalar : constBucketValues.scalars()) {
            NGP_EXPECT_EQ(constBucketValues(entity, scalar), (bucketId*1000 + entity()*100 + scalar()*10));
          }
        }
      );
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiScalar_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_scalar(get_bulk(), field,
                           field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                           field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataBucketAccess, device_multiScalar_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_device_multi_scalar(get_bulk(), fieldBase,
                           fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                           fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
void test_device_multi_scalar_pointer(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = bucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
            bucketValues(entity, scalar) = bucketId*1000 + entity()*100 + scalar()*10;
          }
        }
      );
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const int* constBucketPtr = constBucketValues.pointer();
      const int entityStride = constBucketValues.entity_stride();
      const int scalarStride = constBucketValues.scalar_stride();
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (int scalar = 0; scalar < constBucketValues.num_scalars(); ++scalar) {
            NGP_EXPECT_EQ(constBucketPtr[entity*entityStride + scalar*scalarStride],
                          (bucketId*1000 + entity()*100 + scalar*10));
          }
        }
      );
    }
  );
}

NGP_TEST_F(FieldDataBucketAccess, device_multiScalar_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_multi_scalar_pointer(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_multi_scalar_traditional_for_loop(stk::mesh::BulkData& bulk,
                                                   stk::mesh::Field<int>& field)
{
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(field);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = fieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = bucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ScalarIdx scalar(0); scalar < bucketValues.num_scalars(); ++scalar) {
            bucketValues(entity, scalar) = bucketId*1000 + entity()*100 + scalar()*10;
          }
        }
      );
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValues = constFieldData.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ScalarIdx scalar(0); scalar < constBucketValues.num_scalars(); ++scalar) {
            NGP_EXPECT_EQ(constBucketValues(entity, scalar), (bucketId*1000 + entity()*100 + scalar()*10));
          }
        }
      );
    }
  );

  // Write values to FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValuesBase = fieldDataBase.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = bucketValuesBase.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ScalarIdx scalar(0); scalar < bucketValuesBase.num_scalars(); ++scalar) {
            bucketValuesBase(entity, scalar) = bucketId*10000 + entity()*1000 + scalar()*100;
          }
        }
      );
    }
  );

  // Read const values from FieldBase
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto constBucketValuesBase = constFieldDataBase.bucket_values(bucketId);
      const stk::mesh::EntityIdx numEntities = constBucketValuesBase.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
        [&](stk::mesh::EntityIdx entity) {
          for (stk::mesh::ScalarIdx scalar(0); scalar < constBucketValuesBase.num_scalars(); ++scalar) {
            NGP_EXPECT_EQ(constBucketValuesBase(entity, scalar), (bucketId*10000 + entity()*1000 + scalar()*100));
          }
        }
      );
    }
  );
}

NGP_TEST_F(FieldDataBucketAccess, device_multiScalar_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_multi_scalar_traditional_for_loop(get_bulk(), *m_field);
}


//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_multi_component(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ComponentIdx component : bucketValues.components()) {
          bucketValues(entity, component) = bucket->bucket_id()*100 + entity()*10 + component();
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto constBucketValues = constFieldData.bucket_values(bucketId);
        const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
              NGP_EXPECT_EQ(constBucketValues(entity, component), (bucketId*100 + entity()*10 + component()));
            }
          }
        );
      }
    );
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();
  test_host_to_device_multi_component(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_left_field();
  test_host_to_device_multi_component(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiComponent_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_right_field();
  test_host_to_device_multi_component(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_multi_copy(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
          bucketValues(entity, copy) = bucket->bucket_id()*100 + entity()*10 + copy();
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto constBucketValues = constFieldData.bucket_values(bucketId);
        const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
              NGP_EXPECT_EQ(constBucketValues(entity, copy), (bucketId*100 + entity()*10 + copy()));
            }
          }
        );
      }
    );
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiCopy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();
  test_host_to_device_multi_copy(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiCopy_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_left_field();
  test_host_to_device_multi_copy(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiCopy_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_right_field();
  test_host_to_device_multi_copy(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_multi_copy_multi_component(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
          for (stk::mesh::ComponentIdx component : bucketValues.components()) {
            bucketValues(entity, copy, component) = bucket->bucket_id()*1000 + entity()*100 + copy()*10 + component();
          }
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto constBucketValues = constFieldData.bucket_values(bucketId);
        const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
              for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
                NGP_EXPECT_EQ(constBucketValues(entity, copy, component),
                              (bucketId*1000 + entity()*100 + copy()*10 + component()));
              }
            }
          }
        );
      }
    );
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiCopy_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_host_to_device_multi_copy_multi_component(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiCopy_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_host_to_device_multi_copy_multi_component(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiCopy_multiComponent_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_right_field();
  test_host_to_device_multi_copy_multi_component(get_bulk(), *m_rightField);
}


//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_multi_scalar(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
          bucketValues(entity, scalar) = bucket->bucket_id()*1000 + entity()*100 + scalar()*10;
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto constBucketValues = constFieldData.bucket_values(bucketId);
        const stk::mesh::EntityIdx numEntities = constBucketValues.num_entities();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            for (stk::mesh::ScalarIdx scalar : constBucketValues.scalars()) {
              NGP_EXPECT_EQ(constBucketValues(entity, scalar), (bucketId*1000 + entity()*100 + scalar()*10));
            }
          }
        );
      }
    );
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiScalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_host_to_device_multi_scalar(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiScalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_host_to_device_multi_scalar(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedHostToDevice_multiScalar_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_right_field();
  test_host_to_device_multi_scalar(get_bulk(), *m_rightField);
}


//==============================================================================
template <typename FieldType>
void test_device_to_host_scalar(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto bucketValues = fieldData.bucket_values(bucketId);
        stk::mesh::EntityIdx numEntities = bucketValues.num_entities();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            bucketValues(entity) = bucketId*10 + entity();
          }
        );
      }
    );
  }

  {
    auto constFieldData = field.data();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        EXPECT_EQ(constBucketValues(entity), static_cast<int>(bucket->bucket_id()*10 + entity()));
      }
    }
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();
  test_device_to_host_scalar(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_scalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_left_field();
  test_device_to_host_scalar(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_scalar_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_right_field();
  test_device_to_host_scalar(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_multi_component(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto bucketValues = fieldData.bucket_values(bucketId);
        stk::mesh::EntityIdx numEntities = bucketValues.num_entities();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            for (stk::mesh::ComponentIdx component : bucketValues.components()) {
              bucketValues(entity, component) = bucketId*100 + entity()*10 + component();
            }
          }
        );
      }
    );
  }

  {
    auto constFieldData = field.data();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
          EXPECT_EQ(constBucketValues(entity, component),
                    static_cast<int>(bucket->bucket_id()*100 + entity()*10 + component()));
        }
      }
    }
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();
  test_device_to_host_multi_component(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_left_field();
  test_device_to_host_multi_component(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiComponent_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_right_field();
  test_device_to_host_multi_component(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_multi_copy(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto bucketValues = fieldData.bucket_values(bucketId);
        stk::mesh::EntityIdx numEntities = bucketValues.num_entities();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
              bucketValues(entity, copy) = bucketId*100 + entity()*10 + copy();
            }
          }
        );
      }
    );
  }

  {
    auto constFieldData = field.template data<>();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
          EXPECT_EQ(constBucketValues(entity, copy), static_cast<int>(bucket->bucket_id()*100 + entity()*10 + copy()));
        }
      }
    }
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiCopy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();
  test_device_to_host_multi_copy(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiCopy_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_left_field();
  test_device_to_host_multi_copy(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiCopy_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_right_field();
  test_device_to_host_multi_copy(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_multi_copy_multi_component(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto bucketValues = fieldData.bucket_values(bucketId);
        const stk::mesh::EntityIdx numEntities = bucketValues.num_entities();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
              for (stk::mesh::ComponentIdx component : bucketValues.components()) {
                bucketValues(entity, copy, component) = bucketId*1000 + entity()*100 + copy()*10 + component();
              }
            }
          }
        );
      }
    );
  }

  {
    auto constFieldData = field.data();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
          for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
            EXPECT_EQ(constBucketValues(entity, copy, component),
                      static_cast<int>(bucket->bucket_id()*1000 + entity()*100 + copy()*10 + component()));
          }
        }
      }
    }
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiCopy_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_to_host_multi_copy_multi_component(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiCopy_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_device_to_host_multi_copy_multi_component(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiCopy_multiComponent_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_right_field();
  test_device_to_host_multi_copy_multi_component(get_bulk(), *m_rightField);
}


//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_multi_scalar(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::NODE_RANK, field);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::NODE_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto bucketValues = fieldData.bucket_values(bucketId);
        const stk::mesh::EntityIdx numEntities = bucketValues.num_entities();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numEntities),
          [&](stk::mesh::EntityIdx entity) {
            for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
              bucketValues(entity, scalar) = bucketId*1000 + entity()*100 + scalar()*10;
            }
          }
        );
      }
    );
  }

  {
    auto constFieldData = field.data();
    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : constBucketValues.scalars()) {
            EXPECT_EQ(constBucketValues(entity, scalar),
                      static_cast<int>(bucket->bucket_id()*1000 + entity()*100 + scalar()*10));
        }
      }
    }
  }
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiScalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_to_host_multi_scalar(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiScalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_device_to_host_multi_scalar(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataBucketAccess, mixedDeviceToHost_multiScalar_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_right_field();
  test_device_to_host_multi_scalar(get_bulk(), *m_rightField);
}


//==============================================================================
TEST_F(FieldDataBucketAccess, host_isFieldDefined)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField");
  stk::mesh::Part& part1 = get_meta().declare_part_with_topology("part1", stk::topology::HEX_8);
  stk::mesh::Part& part2 = get_meta().declare_part_with_topology("part2", stk::topology::HEX_8);
  stk::mesh::put_field_on_mesh(*m_field, part1, 3, nullptr);

  stk::io::fill_mesh("generated:2x1x1", get_bulk());
  stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  get_bulk().modification_begin();
  get_bulk().change_entity_parts(elem1, stk::mesh::PartVector{&part1}, stk::mesh::PartVector{});
  get_bulk().change_entity_parts(elem2, stk::mesh::PartVector{&part2}, stk::mesh::PartVector{});
  get_bulk().modification_end();

  stk::mesh::Field<int>& field = *m_field;

  auto fieldData = field.data();

  const stk::mesh::BucketVector& buckets1 = get_bulk().get_buckets(stk::topology::ELEM_RANK, part1);
  const stk::mesh::BucketVector& buckets2 = get_bulk().get_buckets(stk::topology::ELEM_RANK, part2);
  auto bucket1Values = fieldData.bucket_values(*buckets1.front());
  auto bucket2Values = fieldData.bucket_values(*buckets2.front());

  EXPECT_EQ(bucket1Values.is_field_defined(), true);
  EXPECT_EQ(bucket2Values.is_field_defined(), false);
}

//------------------------------------------------------------------------------
void device_is_field_defined(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field,
                             stk::mesh::Part& part1, stk::mesh::Part& part2)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, part1,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto bucketValues = fieldData.bucket_values(entity.bucket_id);
      NGP_EXPECT_EQ(bucketValues.is_field_defined(), true);
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, part2,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto bucketValues = fieldData.bucket_values(entity.bucket_id);
      NGP_EXPECT_EQ(bucketValues.is_field_defined(), false);
    }
  );
}

TEST_F(FieldDataBucketAccess, device_isFieldDefined)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField");
  stk::mesh::Part& part1 = get_meta().declare_part_with_topology("part1", stk::topology::HEX_8);
  stk::mesh::Part& part2 = get_meta().declare_part_with_topology("part2", stk::topology::HEX_8);
  stk::mesh::put_field_on_mesh(*m_field, part1, 3, nullptr);

  stk::io::fill_mesh("generated:2x1x1", get_bulk());
  stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);

  get_bulk().modification_begin();
  get_bulk().change_entity_parts(elem1, stk::mesh::PartVector{&part1}, stk::mesh::PartVector{});
  get_bulk().change_entity_parts(elem2, stk::mesh::PartVector{&part2}, stk::mesh::PartVector{});
  get_bulk().modification_end();

  stk::mesh::Field<int>& field = *m_field;

  device_is_field_defined(get_bulk(), field, part1, part2);
}


// Note: We can only test host-side bounds-checking because the device-side checking
// issues a Kokkos::abort() which cannot be trapped in a test.

#ifdef STK_FIELD_BOUNDS_CHECK

//==============================================================================
TEST_F(FieldDataBucketAccess, host_consistencyCheck_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  stk::mesh::Field<int>& elemField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField1");
  stk::mesh::Field<int>& nodeField = get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodeField1");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(nodeField, get_meta().universal_part(), nullptr);
  create_single_element_mesh();
  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);

  const stk::mesh::Bucket& bucket_elem1 = get_bulk().bucket(elem1);
  const stk::mesh::Bucket& bucket_node1 = get_bulk().bucket(node1);

  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadWrite>().bucket_values(bucket_node1));                       // Wrong rank Bucket
  EXPECT_ANY_THROW(elemField.data().bucket_values(bucket_node1));  // Wrong rank Bucket

  // Acquire FieldData before opening modification cycle
  {
    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().modification_begin();
    get_bulk().declare_node(100);
    get_bulk().modification_end();

    EXPECT_NO_THROW(elemFieldData.bucket_values(bucket_elem1));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.bucket_values(bucket_elem1));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.bucket_values(bucket_node1));       // Stale FieldData from before mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.bucket_values(bucket_node1));  // Stale FieldData from before mesh mod
  }

  // Acquire FieldData during modification cycle
  {
    get_bulk().modification_begin();

    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().declare_node(101);

    EXPECT_NO_THROW(elemFieldData.bucket_values(bucket_elem1));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.bucket_values(bucket_elem1));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.bucket_values(bucket_node1));       // Stale FieldData during mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.bucket_values(bucket_node1));  // Stale FieldData during mesh mod

    get_bulk().modification_end();
  }

  auto secondMesh = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_initial_bucket_capacity(1)
                                                          .set_maximum_bucket_capacity(1).create();
  stk::io::fill_mesh("generated:2x1x1", *secondMesh);  // Create two-element mesh
  const stk::mesh::Entity elem2 = secondMesh->get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::Bucket& bucket_elem2 = secondMesh->bucket(elem2);
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadWrite>().bucket_values(bucket_elem2));                       // Bucket from different mesh
  EXPECT_ANY_THROW(elemField.data().bucket_values(bucket_elem2));  // Bucket from different mesh
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_consistencyCheck_bucketId)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  stk::mesh::Field<int>& elemField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField1");
  stk::mesh::Field<int>& nodeField = get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodeField1");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(nodeField, get_meta().universal_part(), nullptr);
  create_single_element_mesh();

  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadWrite>().bucket_values(1));                       // Bad Bucket Id
  EXPECT_ANY_THROW(elemField.data().bucket_values(1));  // Bad Bucket Id

  const int bucketId_elem1 = 0;
  const int bucketId_node1 = 0;

  // Acquire FieldData before opening modification cycle
  {
    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().modification_begin();
    get_bulk().declare_node(100);
    get_bulk().modification_end();

    EXPECT_NO_THROW(elemFieldData.bucket_values(bucketId_elem1));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.bucket_values(bucketId_elem1));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.bucket_values(bucketId_node1));       // Stale FieldData from before mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.bucket_values(bucketId_node1));  // Stale FieldData from before mesh mod
  }

  // Acquire FieldData during modification cycle
  {
    get_bulk().modification_begin();

    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().declare_node(101);

    EXPECT_NO_THROW(elemFieldData.bucket_values(bucketId_elem1));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.bucket_values(bucketId_elem1));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.bucket_values(bucketId_node1));       // Stale FieldData during mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.bucket_values(bucketId_node1));  // Stale FieldData during mesh mod

    get_bulk().modification_end();
  }
}

//==============================================================================
TEST_F(FieldDataBucketAccess, host_scalarField_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);
  const stk::mesh::Bucket& bucket = get_bulk().bucket(node1);  // Bucket holds both nodes

  auto bucketValues = field.data<stk::mesh::ReadWrite>().bucket_values(bucket);

  const stk::mesh::EntityIdx goodEntity = 1_entity;  // Only 2 entities
  const stk::mesh::EntityIdx badEntity  = 2_entity;

  const stk::mesh::CopyIdx goodCopy = 0_copy;  // Only 1 copy
  const stk::mesh::CopyIdx badCopy  = 1_copy;

  const stk::mesh::ComponentIdx goodComponent = 0_comp;  // Only 1 components
  const stk::mesh::ComponentIdx badComponent  = 1_comp;

  EXPECT_NO_THROW(bucketValues(goodEntity));                           // Normal scalar access
  EXPECT_NO_THROW(bucketValues(goodEntity, goodComponent));            // In-bounds component
  EXPECT_NO_THROW(bucketValues(goodEntity, goodComponent));            // In-bounds copy
  EXPECT_NO_THROW(bucketValues(goodEntity, goodCopy, goodComponent));  // In-bounds copy/component

  EXPECT_ANY_THROW(bucketValues(badEntity));                           // Out-of-bounds entity in scalar access

  EXPECT_ANY_THROW(bucketValues(badEntity, goodComponent));            // Out-of-bounds entity in component access
  EXPECT_ANY_THROW(bucketValues(goodEntity, badComponent));            // Out-of-bounds component

  EXPECT_ANY_THROW(bucketValues(badEntity, goodCopy));                 // Out-of-bounds entity in copy access
  EXPECT_ANY_THROW(bucketValues(goodEntity, badCopy));                 // Out-of-bounds copy

  EXPECT_ANY_THROW(bucketValues(badEntity, goodCopy, goodComponent));  // Out-of-bounds entity in copy/component access
  EXPECT_ANY_THROW(bucketValues(goodEntity, badCopy, goodComponent));  // Out-of-bounds copy/component
  EXPECT_ANY_THROW(bucketValues(goodEntity, goodCopy, badComponent));  // Out-of-bounds copy/component
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiComponentField_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), 3, nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);
  const stk::mesh::Bucket& bucket = get_bulk().bucket(node1);  // Bucket holds both nodes

  auto bucketValues = field.data<stk::mesh::ReadWrite>().bucket_values(bucket);

  const stk::mesh::EntityIdx goodEntity = 1_entity;  // Only 2 entities
  const stk::mesh::EntityIdx badEntity  = 2_entity;

  const stk::mesh::CopyIdx goodCopy = 0_copy;  // Only 1 copy
  const stk::mesh::CopyIdx badCopy  = 1_copy;

  const stk::mesh::ComponentIdx goodComponent = 2_comp;  // Only 3 components
  const stk::mesh::ComponentIdx badComponent  = 3_comp;

  EXPECT_NO_THROW(bucketValues(goodEntity, goodComponent));            // In-bounds component
  EXPECT_NO_THROW(bucketValues(goodEntity, goodCopy, goodComponent));  // In-bounds copy/component

  EXPECT_ANY_THROW(bucketValues(goodEntity));                          // Mistaken scalar access

  EXPECT_ANY_THROW(bucketValues(badEntity, goodComponent));            // Out-of-bounds entity in component access
  EXPECT_ANY_THROW(bucketValues(goodEntity, badComponent));            // Out-of-bounds component

  EXPECT_ANY_THROW(bucketValues(badEntity, goodCopy));                 // Out-of-bounds entity in copy access
  EXPECT_ANY_THROW(bucketValues(goodEntity, goodCopy));                // Mistaken copy access

  EXPECT_ANY_THROW(bucketValues(badEntity, goodCopy, goodComponent));  // Out-of-bounds entity in copy/component access
  EXPECT_ANY_THROW(bucketValues(goodEntity, badCopy, goodComponent));  // Out-of-bounds copy
  EXPECT_ANY_THROW(bucketValues(goodEntity, goodCopy, badComponent));  // Out-of-bounds component
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopyField_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), 1, 8, nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);
  const stk::mesh::Bucket& bucket = get_bulk().bucket(node1);  // Bucket holds both nodes

  auto bucketValues = field.data<stk::mesh::ReadWrite>().bucket_values(bucket);

  const stk::mesh::EntityIdx goodEntity = 1_entity;  // Only 2 entities
  const stk::mesh::EntityIdx badEntity  = 2_entity;

  const stk::mesh::CopyIdx goodCopy = 7_copy;  // Only 8 copy
  const stk::mesh::CopyIdx badCopy  = 8_copy;

  const stk::mesh::ComponentIdx goodComponent = 0_comp;  // Only 1 component
  const stk::mesh::ComponentIdx badComponent  = 1_comp;

  EXPECT_NO_THROW(bucketValues(goodEntity, goodCopy));                 // In-bounds copy
  EXPECT_NO_THROW(bucketValues(goodEntity, goodCopy, goodComponent));  // In-bounds copy/component

  EXPECT_ANY_THROW(bucketValues(goodEntity));                          // Mistaken scalar access

  EXPECT_ANY_THROW(bucketValues(badEntity, goodComponent));            // Out-of-bounds entity in component access
  EXPECT_ANY_THROW(bucketValues(goodEntity, goodComponent));           // Mistaken component access

  EXPECT_ANY_THROW(bucketValues(badEntity, goodCopy));                 // Out-of-bounds entity in copy access
  EXPECT_ANY_THROW(bucketValues(goodEntity, badCopy));                 // Out-of-bounds copy

  EXPECT_ANY_THROW(bucketValues(badEntity, goodCopy, goodComponent));  // Out-of-bounds entity in copy/component access
  EXPECT_ANY_THROW(bucketValues(goodEntity, badCopy, goodComponent));  // Out-of-bounds copy
  EXPECT_ANY_THROW(bucketValues(goodEntity, goodCopy, badComponent));  // Out-of-bounds component
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiCopyMultiComponentField_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), 3, 8, nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);
  const stk::mesh::Bucket& bucket = get_bulk().bucket(node1);  // Bucket holds both nodes

  auto bucketValues = field.data<stk::mesh::ReadWrite>().bucket_values(bucket);

  const stk::mesh::EntityIdx goodEntity = 1_entity;  // Only 2 entities
  const stk::mesh::EntityIdx badEntity  = 2_entity;

  const stk::mesh::CopyIdx goodCopy = 7_copy;  // Only 8 copies
  const stk::mesh::CopyIdx badCopy  = 8_copy;

  const stk::mesh::ComponentIdx goodComponent = 2_comp;  // Only 3 components
  const stk::mesh::ComponentIdx badComponent  = 3_comp;

  EXPECT_NO_THROW(bucketValues(goodEntity, goodCopy, goodComponent));  // In-bounds copy/component

  EXPECT_ANY_THROW(bucketValues(goodEntity));                          // Mistaken scalar access

  EXPECT_ANY_THROW(bucketValues(badEntity, goodComponent));            // Out-of-bounds entity in component access
  EXPECT_ANY_THROW(bucketValues(goodEntity, goodComponent));           // Mistaken component acces

  EXPECT_ANY_THROW(bucketValues(badEntity, goodCopy));                 // Out-of-bounds entity in copy access
  EXPECT_ANY_THROW(bucketValues(goodEntity, goodCopy));                // Mistaken copy access

  EXPECT_ANY_THROW(bucketValues(badEntity, goodCopy, goodComponent));  // Out-of-bounds entity in copy/component access
  EXPECT_ANY_THROW(bucketValues(goodEntity, badCopy, goodComponent));  // Out-of-bounds copy
  EXPECT_ANY_THROW(bucketValues(goodEntity, goodCopy, badComponent));  // Out-of-bounds component
}

//------------------------------------------------------------------------------
TEST_F(FieldDataBucketAccess, host_multiScalarField_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), 3, 8, nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);
  const stk::mesh::Bucket& bucket = get_bulk().bucket(node1);  // Bucket holds both nodes

  auto bucketValues = field.data<stk::mesh::ReadWrite>().bucket_values(bucket);

  const stk::mesh::EntityIdx goodEntity = 1_entity;  // Only 2 entities

  const stk::mesh::ScalarIdx goodScalar = 23_scalar;  // Only 24 scalars
  const stk::mesh::ScalarIdx badScalar  = 24_scalar;

  EXPECT_NO_THROW(bucketValues(goodEntity, goodScalar));  // In-bounds scalar
  EXPECT_ANY_THROW(bucketValues(goodEntity, badScalar));  // Out-of-bounds scalar
}

//==============================================================================

#endif  // STK_FIELD_BOUNDS_CHECK

} //namespace <anonymous>
