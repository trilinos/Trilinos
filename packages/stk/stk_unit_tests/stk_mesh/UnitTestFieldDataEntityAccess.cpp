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
#include <stk_io/FillMesh.hpp>

namespace {

class FieldDataEntityAccess : public FieldDataAccessFixture {};

//==============================================================================
TEST_F(FieldDataEntityAccess, inconsistentTemplateParameters)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::HostSpace, stk::mesh::Layout::Left>()));  // Correct
  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>()));  // Correct

  // We only do a static_cast() in release builds and don't do any template parameter checking, so this will
  // not detect any problems and may lead to memory corruption.
#ifdef STK_FIELD_BOUNDS_CHECK
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadOnly,  stk::ngp::HostSpace, stk::mesh::Layout::Left>()));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>()));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::HostSpace, stk::mesh::Layout::Right>()));  // Wrong layout
  EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>()));  // Wrong layout
#endif

  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::DeviceSpace, stk::mesh::Layout::Left>()));  // Correct
  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace, stk::mesh::Layout::Left>()));  // Correct

#ifdef STK_FIELD_BOUNDS_CHECK
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadOnly,  stk::ngp::DeviceSpace, stk::mesh::Layout::Left>()));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::DeviceSpace, stk::mesh::Layout::Left>()));  // Wrong datatype
  //EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::DeviceSpace, stk::mesh::Layout::Right>()));  // Trapped by static_assert()
  //EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace, stk::mesh::Layout::Right>()));  // Trapped by static_assert()
#endif
}

TEST_F(FieldDataEntityAccess, inconsistentTemplateParameters_async)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::HostSpace,
                                  stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Correct
  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace,
                                  stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Correct

  // We only do a static_cast() in release builds and don't do any template parameter checking, so this will
  // not detect any problems and may lead to memory corruption.
#ifdef STK_FIELD_BOUNDS_CHECK
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadOnly, stk::ngp::HostSpace,
                                   stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::HostSpace,
                                   stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::HostSpace,
                                   stk::mesh::Layout::Right>(stk::ngp::ExecSpace())));  // Wrong layout
  EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace,
                                   stk::mesh::Layout::Right>(stk::ngp::ExecSpace())));  // Wrong layout
#endif

  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::DeviceSpace,
                                  stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Correct
  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace,
                                  stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Correct

#ifdef STK_FIELD_BOUNDS_CHECK
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadOnly,  stk::ngp::DeviceSpace,
                                   stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::DeviceSpace,
                                   stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Wrong datatype
  // EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::DeviceSpace,
  //                                  stk::mesh::Layout::Right>(stk::ngp::ExecSpace())));  // Trapped by static_assert()
  // EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace,
  //                                  stk::mesh::Layout::Right>(stk::ngp::ExecSpace())));  // Trapped by static_assert()
#endif
}

//==============================================================================
template <typename FieldDataType, typename ConstFieldDataType>
void test_host_scalar(const stk::mesh::BulkData& bulk,
                      const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  // Entity indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        entityValues() = ++value*10;
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        EXPECT_EQ(constEntityValues(), ++value*10);
      }
    }
  }

  // MeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto entityValues = fieldData.entity_values(mi);
        entityValues() = ++value*20;
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto constEntityValues = constFieldData.entity_values(mi);
        EXPECT_EQ(constEntityValues(), ++value*20);
      }
    }
  }

  // FastMeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto entityValues = fieldData.entity_values(fmi);
        entityValues() = ++value*30;
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto constEntityValues = constFieldData.entity_values(fmi);
        EXPECT_EQ(constEntityValues(), ++value*30);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_scalar(get_bulk(),
                   field.data<stk::mesh::ReadWrite>(),
                   field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite>(),
                   fieldBase.data<int>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_field_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_scalar(get_bulk(),
                   field.data<stk::mesh::ReadWrite>(),
                   field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_fieldBase_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_field_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_scalar(get_bulk(),
                   field.data<stk::mesh::ReadWrite>(),
                   field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_fieldBase_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_fieldBase_layoutLeft_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_fieldBase_layoutAuto_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_fieldBase_layoutRight_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_scalar(get_bulk(),
                   fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                   fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_fieldBase_layoutAuto_layoutRight)
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
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        entityValues() = ++value*10;
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        const int* constEntityPtr = constEntityValues.pointer();
        EXPECT_EQ(*constEntityPtr, ++value*10);
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        int* entityPtr = entityValues.pointer();
        *entityPtr = ++value*20;
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        EXPECT_EQ(constEntityValues(), ++value*20);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_scalar_pointer(get_bulk(),
                           field.data<stk::mesh::ReadWrite>(),
                           field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_scalar_pointer(get_bulk(),
                           field.data<stk::mesh::ReadWrite>(),
                           field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_pointer_layoutRight)
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

  // Entity indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(component) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          EXPECT_EQ(constEntityValues(component), ++value*10);
        }
      }
    }
  }

  // MeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto entityValues = fieldData.entity_values(mi);
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(component) = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto constEntityValues = constFieldData.entity_values(mi);
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          EXPECT_EQ(constEntityValues(component), ++value*20);
        }
      }
    }
  }

  // FastMeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto entityValues = fieldData.entity_values(fmi);
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(component) = ++value*30;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto constEntityValues = constFieldData.entity_values(fmi);
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          EXPECT_EQ(constEntityValues(component), ++value*30);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_component(get_bulk(),
                            field.data<stk::mesh::ReadWrite>(),
                            field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite>(),
                            fieldBase.data<int>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_field_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_component(get_bulk(),
                            field.data<stk::mesh::ReadWrite>(),
                            field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_fieldBase_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_field_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_component(get_bulk(),
                            field.data<stk::mesh::ReadWrite>(),
                            field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_fieldBase_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_fieldBase_layoutLeft_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_fieldBase_layoutAuto_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_fieldBase_layoutRight_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_component(get_bulk(),
                            fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                            fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_fieldBase_layoutAuto_layoutRight)
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
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(component) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        const int* constEntityPtr = constEntityValues.pointer();
        const int componentStride = constEntityValues.component_stride();
        for (int component = 0; component < constEntityValues.num_components(); ++component) {
          EXPECT_EQ(constEntityPtr[component*componentStride], ++value*10);
        }
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        int* entityPtr = entityValues.pointer();
        const int componentStride = entityValues.component_stride();
        for (int component = 0; component < entityValues.num_components(); ++component) {
          entityPtr[component*componentStride] = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          EXPECT_EQ(constEntityValues(component), ++value*20);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_component_pointer(get_bulk(),
                                    field.data<stk::mesh::ReadWrite>(),
                                    field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_component_pointer(get_bulk(),
                                    field.data<stk::mesh::ReadWrite>(),
                                    field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_pointer_layoutRight)
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

  // Entity indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          entityValues(copy) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          EXPECT_EQ(constEntityValues(copy), ++value*10);
        }
      }
    }
  }

  // MeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto entityValues = fieldData.entity_values(mi);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          entityValues(copy) = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto constEntityValues = constFieldData.entity_values(mi);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          EXPECT_EQ(constEntityValues(copy), ++value*20);
        }
      }
    }
  }

  // FastMeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto entityValues = fieldData.entity_values(fmi);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          entityValues(copy) = ++value*30;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto constEntityValues = constFieldData.entity_values(fmi);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          EXPECT_EQ(constEntityValues(copy), ++value*30);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_copy(get_bulk(),
                       field.data<stk::mesh::ReadWrite>(),
                       field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite>(),
                       fieldBase.data<int>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_field_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_copy(get_bulk(),
                       field.data<stk::mesh::ReadWrite>(),
                       field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_fieldBase_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_field_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_right_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Right>& field = *m_rightField;

  test_host_multi_copy(get_bulk(),
                       field.data<stk::mesh::ReadWrite>(),
                       field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_fieldBase_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Right>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_fieldBase_layoutLeft_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Left>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_fieldBase_layoutAuto_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Auto>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Left>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_fieldBase_layoutRight_layoutAuto)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_right_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_rightField);

  test_host_multi_copy(get_bulk(),
                       fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostSpace, stk::mesh::Layout::Right>(),
                       fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostSpace, stk::mesh::Layout::Auto>());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_fieldBase_layoutAuto_layoutRight)
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
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          entityValues(copy) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        const int* constEntityPtr = constEntityValues.pointer();
        const int copyStride = constEntityValues.copy_stride();
        for (int copy = 0; copy < constEntityValues.num_copies(); ++copy) {
          EXPECT_EQ(constEntityPtr[copy*copyStride], ++value*10);
        }
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        int* entityPtr = entityValues.pointer();
        const int copyStride = entityValues.copy_stride();
        for (int copy = 0; copy < entityValues.num_copies(); ++copy) {
          entityPtr[copy*copyStride] = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          EXPECT_EQ(constEntityValues(copy), ++value*20);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_copy_pointer(get_bulk(),
                               field.data<stk::mesh::ReadWrite>(),
                               field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_copy_pointer(get_bulk(),
                               field.data<stk::mesh::ReadWrite>(),
                               field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_pointer_layoutRight)
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

  // Entity indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            entityValues(copy, component) = ++value*10;
          }
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
            EXPECT_EQ(constEntityValues(copy, component), ++value*10);
          }
        }
      }
    }
  }

  // MeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto entityValues = fieldData.entity_values(mi);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            entityValues(copy, component) = ++value*20;
          }
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto constEntityValues = constFieldData.entity_values(mi);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
            EXPECT_EQ(constEntityValues(copy, component), ++value*20);
          }
        }
      }
    }
  }

  // FastMeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto entityValues = fieldData.entity_values(fmi);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            entityValues(copy, component) = ++value*30;
          }
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto constEntityValues = constFieldData.entity_values(fmi);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
            EXPECT_EQ(constEntityValues(copy, component), ++value*30);
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_field)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_fieldBase)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_field_layoutLeft)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_fieldBase_layoutLeft)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_field_layoutRight)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_fieldBase_layoutRight)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_fieldBase_layoutLeft_layoutAuto)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_fieldBase_layoutAuto_layoutLeft)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_fieldBase_layoutRight_layoutAuto)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_fieldBase_layoutAuto_layoutRight)
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  stk::mesh::Field<int>& field = *m_field;
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  auto fieldData = field.data<stk::mesh::ReadWrite>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < entityValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < entityValues.num_components(); ++component) {
          entityValues(copy, component) = ++value;
        }
      }
    }
  }

  value = 0;
  auto constFieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < constEntityValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constEntityValues.num_components(); ++component) {
          EXPECT_EQ(constEntityValues(copy, component), ++value);
        }
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < entityValuesBase.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < entityValuesBase.num_components(); ++component) {
          entityValuesBase(copy, component) = ++value*10;
        }
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < constEntityValuesBase.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constEntityValuesBase.num_components(); ++component) {
          EXPECT_EQ(constEntityValuesBase(copy, component), ++value*10);
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
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            entityValues(copy, component) = ++value*10;
          }
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        const int* constEntityPtr = constEntityValues.pointer();
        const int copyStride = constEntityValues.copy_stride();
        const int componentStride = constEntityValues.component_stride();
        for (int copy = 0; copy < constEntityValues.num_copies(); ++copy) {
          for (int component = 0; component < constEntityValues.num_components(); ++component) {
            EXPECT_EQ(constEntityPtr[copy*copyStride + component*componentStride], ++value*10);
          }
        }
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        int* entityPtr = entityValues.pointer();
        const int copyStride = entityValues.copy_stride();
        const int componentStride = entityValues.component_stride();
        for (int copy = 0; copy < entityValues.num_copies(); ++copy) {
          for (int component = 0; component < entityValues.num_components(); ++component) {
            entityPtr[copy*copyStride + component*componentStride] = ++value*20;
          }
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
            EXPECT_EQ(constEntityValues(copy, component), ++value*20);
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_copy_multi_component_pointer(get_bulk(),
                                               field.data<stk::mesh::ReadWrite>(),
                                               field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_copy_multi_component_pointer(get_bulk(),
                                               field.data<stk::mesh::ReadWrite>(),
                                               field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_pointer_layoutRight)
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

  // Entity indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::ScalarIdx scalar : entityValues.scalars()) {
          entityValues(scalar) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::ScalarIdx scalar : constEntityValues.scalars()) {
          EXPECT_EQ(constEntityValues(scalar), ++value*10);
        }
      }
    }
  }

  // MeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto entityValues = fieldData.entity_values(mi);
        for (stk::mesh::ScalarIdx scalar : entityValues.scalars()) {
          entityValues(scalar) = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        stk::mesh::MeshIndex mi{bucket, bucketOrd};
        auto constEntityValues = constFieldData.entity_values(mi);
        for (stk::mesh::ScalarIdx scalar : constEntityValues.scalars()) {
          EXPECT_EQ(constEntityValues(scalar), ++value*20);
        }
      }
    }
  }

  // FastMeshIndex indexing
  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto entityValues = fieldData.entity_values(fmi);
        for (stk::mesh::ScalarIdx scalar : entityValues.scalars()) {
          entityValues(scalar) = ++value*30;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (unsigned bucketOrd = 0; bucketOrd < bucket->size();  ++bucketOrd) {
        const stk::mesh::FastMeshIndex fmi{bucket->bucket_id(), bucketOrd};
        auto constEntityValues = constFieldData.entity_values(fmi);
        for (stk::mesh::ScalarIdx scalar : constEntityValues.scalars()) {
          EXPECT_EQ(constEntityValues(scalar), ++value*30);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiScalar_field)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_fieldBase)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_field_layoutLeft)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_fieldBase_layoutLeft)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_field_layoutRight)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_fieldBase_layoutRight)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_fieldBase_layoutLeft_layoutAuto)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_fieldBase_layoutAuto_layoutLeft)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_fieldBase_layoutRight_layoutAuto)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_fieldBase_layoutAuto_layoutRight)
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
TEST_F(FieldDataEntityAccess, host_multiScalar_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  stk::mesh::Field<int>& field = *m_field;
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  auto fieldData = field.data<stk::mesh::ReadWrite>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar(0); scalar < entityValues.num_scalars(); ++scalar) {
        entityValues(scalar) = ++value;
      }
    }
  }

  value = 0;
  auto constFieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar(0); scalar < constEntityValues.num_scalars(); ++scalar) {
        EXPECT_EQ(constEntityValues(scalar), ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar(0); scalar < entityValuesBase.num_scalars(); ++scalar) {
        entityValuesBase(scalar) = ++value*10;
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar(0); scalar < constEntityValuesBase.num_scalars(); ++scalar) {
        EXPECT_EQ(constEntityValuesBase(scalar), ++value*10);
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
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::ScalarIdx scalar : entityValues.scalars()) {
          entityValues(scalar) = ++value*10;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        const int* constEntityPtr = constEntityValues.pointer();
        const int scalarStride = constEntityValues.scalar_stride();
        for (int scalar = 0; scalar < constEntityValues.num_scalars(); ++scalar) {
          EXPECT_EQ(constEntityPtr[scalar*scalarStride], ++value*10);
        }
      }
    }
  }

  {
    int value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto entityValues = fieldData.entity_values(entity);
        int* entityPtr = entityValues.pointer();
        const int scalarStride = entityValues.scalar_stride();
        for (int scalar = 0; scalar < entityValues.num_scalars(); ++scalar) {
          entityPtr[scalar*scalarStride] = ++value*20;
        }
      }
    }

    value = 0;
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::ScalarIdx scalar : constEntityValues.scalars()) {
          EXPECT_EQ(constEntityValues(scalar), ++value*20);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiScalar_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_host_multi_scalar_pointer(get_bulk(),
                                 field.data<stk::mesh::ReadWrite>(),
                                 field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiScalar_pointer_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;

  test_host_multi_scalar_pointer(get_bulk(),
                                 field.data<stk::mesh::ReadWrite>(),
                                 field.data());
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiScalar_pointer_layoutRight)
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

  // Entity indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto entityValues = fieldData.entity_values(entity);
      entityValues() = (fmi.bucket_id*100 + fmi.bucket_ord*10);
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto constEntityValues = constFieldData.entity_values(entity);
      NGP_EXPECT_EQ(constEntityValues(), static_cast<int>(fmi.bucket_id*100 + fmi.bucket_ord*10));
    }
  );

  // FastMeshIndex indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto entityValues = fieldData.entity_values(fmi);
      entityValues() = (fmi.bucket_id*200 + fmi.bucket_ord*20);
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto constEntityValues = constFieldData.entity_values(fmi);
      NGP_EXPECT_EQ(constEntityValues(), static_cast<int>(fmi.bucket_id*200 + fmi.bucket_ord*20));
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_scalar_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_scalar(get_bulk(), field,
                     field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                     field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_scalar_fieldBase)
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

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      entityValues() = (entity.bucket_id*10 + entity.bucket_ord);
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      NGP_EXPECT_EQ(*constEntityPtr, static_cast<int>(entity.bucket_id*10 + entity.bucket_ord));
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_scalar_pointer)
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

  // Entity indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        entityValues(component) = fmi.bucket_id*100 + fmi.bucket_ord*10 + component();
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
        NGP_EXPECT_EQ(constEntityValues(component),
                      static_cast<int>(fmi.bucket_id*100 + fmi.bucket_ord*10 + component()));
      }
    }
  );

  // FastMeshIndex indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto entityValues = fieldData.entity_values(fmi);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        entityValues(component) = fmi.bucket_id*200 + fmi.bucket_ord*20 + component()*2;
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto constEntityValues = constFieldData.entity_values(fmi);
      for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
        NGP_EXPECT_EQ(constEntityValues(component),
                      static_cast<int>(fmi.bucket_id*200 + fmi.bucket_ord*20 + component()*2));
      }
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiComponent_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_component(get_bulk(), field,
                              field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                              field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiComponent_fieldBase)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);

  test_device_multi_component(get_bulk(), fieldBase,
                              fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                              fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiComponent_field_async)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_component(get_bulk(), field,
                              field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(stk::ngp::ExecSpace()),
                              field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(stk::ngp::ExecSpace()));
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiComponent_fieldBase_async)
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

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        entityValues(component) = entity.bucket_id*100 + entity.bucket_ord*10 + component();
      }
    }
  );

  // Read const value from Field<int>
  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      const int componentStride = constEntityValues.component_stride();
      for (int component = 0; component < constEntityValues.num_components(); ++component) {
        NGP_EXPECT_EQ(constEntityPtr[component*componentStride],
                      static_cast<int>(entity.bucket_id*100 + entity.bucket_ord*10 + component));
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiComponent_pointer)
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

  // Entity indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        entityValues(copy) = fmi.bucket_id*100 + fmi.bucket_ord*10 + copy();
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        NGP_EXPECT_EQ(constEntityValues(copy), static_cast<int>(fmi.bucket_id*100 + fmi.bucket_ord*10 + copy()));
      }
    }
  );

  // FastMeshIndex indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto entityValues = fieldData.entity_values(fmi);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        entityValues(copy) = fmi.bucket_id*200 + fmi.bucket_ord*20 + copy()*2;
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto constEntityValues = constFieldData.entity_values(fmi);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        NGP_EXPECT_EQ(constEntityValues(copy), static_cast<int>(fmi.bucket_id*200 + fmi.bucket_ord*20 + copy()*2));
      }
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_copy(get_bulk(), field,
                         field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                         field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_fieldBase)
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

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        entityValues(copy) = entity.bucket_id*100 + entity.bucket_ord*10 + copy();
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      const int copyStride = constEntityValues.copy_stride();
      for (int copy = 0; copy < constEntityValues.num_copies(); ++copy) {
        NGP_EXPECT_EQ(constEntityPtr[copy*copyStride],
                      static_cast<int>(entity.bucket_id*100 + entity.bucket_ord*10 + copy));
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_pointer)
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

  // Entity indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(copy, component) = fmi.bucket_id*1000 + fmi.bucket_ord*100 + copy()*10 + component();
        }
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          NGP_EXPECT_EQ(constEntityValues(copy, component),
                        static_cast<int>(fmi.bucket_id*1000 + fmi.bucket_ord*100 + copy()*10 + component()));
        }
      }
    }
  );

  // FastMeshIndex indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto entityValues = fieldData.entity_values(fmi);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(copy, component) = fmi.bucket_id*2000 + fmi.bucket_ord*200 + copy()*20 + component()*2;
        }
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto constEntityValues = constFieldData.entity_values(fmi);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          NGP_EXPECT_EQ(constEntityValues(copy, component),
                        static_cast<int>(fmi.bucket_id*2000 + fmi.bucket_ord*200 + copy()*20 + component()*2));
        }
      }
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_multiComponent_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_copy_multi_component(get_bulk(), field,
                                         field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                                         field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_multiComponent_fieldBase)
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

  // Write and read values from Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(copy, component) = entity.bucket_id*1000 + entity.bucket_ord*100 + copy()*10 + component();
        }
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      const int copyStride = constEntityValues.copy_stride();
      const int componentStride = constEntityValues.component_stride();
      for (int copy = 0; copy < constEntityValues.num_copies(); ++copy) {
        for (int component = 0; component < constEntityValues.num_components(); ++component) {
          NGP_EXPECT_EQ(constEntityPtr[copy*copyStride + component*componentStride],
                        static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + copy*10 + component));
        }
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_multiComponent_pointer)
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

  // Write and read values from Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < entityValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < entityValues.num_components(); ++component) {
          entityValues(copy, component) = (entity.bucket_id*1000 + entity.bucket_ord*100 + copy()*10 + component());
        }
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < constEntityValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constEntityValues.num_components(); ++component) {
          NGP_EXPECT_EQ(constEntityValues(copy, component),
                        static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + copy()*10 + component()));
        }
      }
    }
  );

  // Write and read values from FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < entityValuesBase.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < entityValuesBase.num_components(); ++component) {
          entityValuesBase(copy, component) = entity.bucket_id*10000 + entity.bucket_ord*1000 + copy()*100 +
                                              component()*10;
        }
      }
    }
  );

  // Read const value from FieldBase
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < constEntityValuesBase.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constEntityValuesBase.num_components(); ++component) {
          NGP_EXPECT_EQ(constEntityValuesBase(copy, component),
                        static_cast<int>(entity.bucket_id*10000 + entity.bucket_ord*1000 + copy()*100 +
                                         component()*10));
        }
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_multiComponent_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_multi_copy_multi_component_traditional_for_loop(get_bulk(), *m_field);
}


//==============================================================================
template <typename FieldType, typename FieldDataType, typename ConstFieldDataType>
void test_device_multi_scalar(const stk::mesh::BulkData& bulk, FieldType& field,
                              const FieldDataType& fieldData, const ConstFieldDataType& constFieldData)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Entity indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar : entityValues.scalars()) {
          entityValues(scalar) = fmi.bucket_id*1000 + fmi.bucket_ord*100 + scalar()*10;
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      const stk::mesh::Entity entity = ngpMesh.get_entity(stk::topology::NODE_RANK, fmi);
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar : constEntityValues.scalars()) {
        NGP_EXPECT_EQ(constEntityValues(scalar),
                      static_cast<int>(fmi.bucket_id*1000 + fmi.bucket_ord*100 + scalar()*10));
      }
    }
  );

  // FastMeshIndex indexing
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto entityValues = fieldData.entity_values(fmi);
      for (stk::mesh::ScalarIdx scalar : entityValues.scalars()) {
        entityValues(scalar) = fmi.bucket_id*2000 + fmi.bucket_ord*200 + scalar()*20;
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
      auto constEntityValues = constFieldData.entity_values(fmi);
      for (stk::mesh::ScalarIdx scalar : constEntityValues.scalars()) {
        NGP_EXPECT_EQ(constEntityValues(scalar),
                      static_cast<int>(fmi.bucket_id*2000 + fmi.bucket_ord*200 + scalar()*20));
      }
    }
  );
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiScalar_field)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  const stk::mesh::Field<int>& field = *m_field;

  test_device_multi_scalar(get_bulk(), field,
                           field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>(),
                           field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>());
}

//------------------------------------------------------------------------------
NGP_TEST_F(FieldDataEntityAccess, device_multiScalar_fieldBase)
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

  // Write and read values from Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar : entityValues.scalars()) {
        entityValues(scalar) = entity.bucket_id*1000 + entity.bucket_ord*100 + scalar()*10;
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      const int scalarStride = constEntityValues.scalar_stride();
      for (int scalar = 0; scalar < constEntityValues.num_scalars(); ++scalar) {
          NGP_EXPECT_EQ(constEntityPtr[scalar*scalarStride],
                        static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + scalar*10));
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiScalar_pointer)
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

  // Write and read values from Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar(0); scalar < entityValues.num_scalars(); ++scalar) {
        entityValues(scalar) = entity.bucket_id*1000 + entity.bucket_ord*100 + scalar()*10;
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar(0); scalar < constEntityValues.num_scalars(); ++scalar) {
        NGP_EXPECT_EQ(constEntityValues(scalar),
                      static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + scalar()*10));
      }
    }
  );

  // Write and read values from FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar(0); scalar < entityValuesBase.num_scalars(); ++scalar) {
        entityValuesBase(scalar) = entity.bucket_id*10000 + entity.bucket_ord*1000 + scalar()*100;
      }
    }
  );

  // Read const value from FieldBase
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::ScalarIdx scalar(0); scalar < constEntityValuesBase.num_scalars(); ++scalar) {
        NGP_EXPECT_EQ(constEntityValuesBase(scalar),
                      static_cast<int>(entity.bucket_id*10000 + entity.bucket_ord*1000 + scalar()*100));
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiScalar_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_multi_scalar_traditional_for_loop(get_bulk(), *m_field);
}


//==============================================================================
template <typename FieldType>
void test_host_to_device_scalar(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto entityValues = fieldData.entity_values(entity);
        entityValues() = nodeId;
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto constEntityValues = constFieldData.entity_values(fmi);
        NGP_EXPECT_EQ(constEntityValues(), nodeId);
      }
    );
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();
  test_host_to_device_scalar(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_scalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_left_field();
  test_host_to_device_scalar(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_scalar_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_right_field();
  test_host_to_device_scalar(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_multi_component(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(component) = nodeId*10 + component;
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto constEntityValues = constFieldData.entity_values(fmi);
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          NGP_EXPECT_EQ(constEntityValues(component), (nodeId*10 + component));
        }
      }
    );
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();
  test_host_to_device_multi_component(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_left_field();
  test_host_to_device_multi_component(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiComponent_layoutRight)
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
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          entityValues(copy) = nodeId*100 + 10*copy;
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto constEntityValues = constFieldData.entity_values(fmi);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          NGP_EXPECT_EQ(constEntityValues(copy), (nodeId*100 + 10*copy));
        }
      }
    );
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiCopy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();
  test_host_to_device_multi_copy(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiCopy_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_left_field();
  test_host_to_device_multi_copy(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiCopy_layoutRight)
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
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            entityValues(copy, component) = nodeId*100 + 10*copy + component();
          }
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto constEntityValues = constFieldData.entity_values(fmi);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
            NGP_EXPECT_EQ(constEntityValues(copy, component), (nodeId*100 + 10*copy + component()));
          }
        }
      }
    );
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiCopy_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_host_to_device_multi_copy_multi_component(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiCopy_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_host_to_device_multi_copy_multi_component(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiCopy_multiComponent_layoutRight)
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
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>(); // Create early so next call is sync instead of update

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::ScalarIdx scalar : entityValues.scalars()) {
          entityValues(scalar) = nodeId*100 + 10*scalar;
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto constEntityValues = constFieldData.entity_values(fmi);
        for (stk::mesh::ScalarIdx scalar : constEntityValues.scalars()) {
          NGP_EXPECT_EQ(constEntityValues(scalar), (nodeId*100 + 10*scalar));
        }
      }
    );
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiScalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_host_to_device_multi_scalar(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiScalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_host_to_device_multi_scalar(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiScalar_layoutRight)
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
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto entityValues = fieldData.entity_values(fmi);
        entityValues() = nodeId;
      }
    );
  }

  {
    auto constFieldData = field.template data<>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto constEntityValues = constFieldData.entity_values(entity);
        NGP_EXPECT_EQ(constEntityValues(), nodeId);
      }
    }
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();
  test_device_to_host_scalar(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_scalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_left_field();
  test_device_to_host_scalar(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_scalar_layoutRight)
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
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto entityValues = fieldData.entity_values(fmi);
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(component) = nodeId*10 + component;
        }
      }
    );
  }

  {
    auto constFieldData = field.template data<>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          NGP_EXPECT_EQ(constEntityValues(component), (nodeId*10 + component));
        }
      }
    }
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();
  test_device_to_host_multi_component(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_left_field();
  test_device_to_host_multi_component(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiComponent_layoutRight)
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
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto entityValues = fieldData.entity_values(fmi);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          entityValues(copy) = nodeId*100 + 10*copy;
        }
      }
    );
  }

  {
    auto constFieldData = field.template data<>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          NGP_EXPECT_EQ(constEntityValues(copy), (nodeId*100 + 10*copy));
        }
      }
    }
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiCopy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();
  test_device_to_host_multi_copy(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiCopy_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_left_field();
  test_device_to_host_multi_copy(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiCopy_layoutRight)
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
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto entityValues = fieldData.entity_values(fmi);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            entityValues(copy, component) = nodeId*100 + 10*copy + component();
          }
        }
      }
    );
  }

  {
    auto constFieldData = field.template data<>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
            NGP_EXPECT_EQ(constEntityValues(copy, component), (nodeId*100 + 10*copy + component()));
          }
        }
      }
    }
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiCopy_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_to_host_multi_copy_multi_component(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiCopy_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_device_to_host_multi_copy_multi_component(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiCopy_multiComponent_layoutRight)
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
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& fmi) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, fmi));
        auto entityValues = fieldData.entity_values(fmi);
        for (stk::mesh::ScalarIdx scalar : entityValues.scalars()) {
          entityValues(scalar) = nodeId*100 + 10*scalar;
        }
      }
    );
  }

  {
    auto constFieldData = field.template data<>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::ScalarIdx scalar : constEntityValues.scalars()) {
          NGP_EXPECT_EQ(constEntityValues(scalar), (nodeId*100 + 10*scalar));
        }
      }
    }
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiScalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_to_host_multi_scalar(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiScalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_device_to_host_multi_scalar(get_bulk(), *m_leftField);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiScalar_layoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_right_field();
  test_device_to_host_multi_scalar(get_bulk(), *m_rightField);
}

//==============================================================================
TEST_F(FieldDataEntityAccess, host_isFieldDefined)
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
  auto elem1Values = fieldData.entity_values(elem1);
  auto elem2Values = fieldData.entity_values(elem2);

  EXPECT_EQ(elem1Values.is_field_defined(), true);
  EXPECT_EQ(elem2Values.is_field_defined(), false);
}

//------------------------------------------------------------------------------
void device_is_field_defined(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field,
                             stk::mesh::Part& part1, stk::mesh::Part& part2)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, part1,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto elemValues = fieldData.entity_values(entity);
      NGP_EXPECT_EQ(elemValues.is_field_defined(), true);
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, part2,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto elemValues = fieldData.entity_values(entity);
      NGP_EXPECT_EQ(elemValues.is_field_defined(), false);
    }
  );
}

TEST_F(FieldDataEntityAccess, device_isFieldDefined)
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
TEST_F(FieldDataEntityAccess, host_consistencyCheck_entity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& elemField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField1");
  stk::mesh::Field<int>& nodeField = get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodeField1");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(nodeField, get_meta().universal_part(), nullptr);
  create_single_element_mesh();
  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);

  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadWrite>().entity_values(node1));                       // Wrong rank entity
  EXPECT_ANY_THROW(elemField.data().entity_values(node1));  // Wrong rank entity

  // Acquire FieldData before opening modification cycle
  {
    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().modification_begin();
    get_bulk().declare_node(100);
    get_bulk().modification_end();

    EXPECT_NO_THROW(elemFieldData.entity_values(elem1));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.entity_values(elem1));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.entity_values(node1));       // Stale FieldData from before mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.entity_values(node1));  // Stale FieldData from before mesh mod
  }

  // Acquire FieldData during modification cycle
  {
    get_bulk().modification_begin();

    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().declare_node(101);

    EXPECT_NO_THROW(elemFieldData.entity_values(elem1));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.entity_values(elem1));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.entity_values(node1));       // Stale FieldData during mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.entity_values(node1));  // Stale FieldData during mesh mod

    get_bulk().modification_end();
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_consistencyCheck_meshIndex)
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
  const stk::mesh::MeshIndex elem1_mi = get_bulk().mesh_index(elem1);
  const stk::mesh::MeshIndex node1_mi = get_bulk().mesh_index(node1);

  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadWrite>().entity_values(node1_mi));                       // Wrong rank entity
  EXPECT_ANY_THROW(elemField.data().entity_values(node1_mi));  // Wrong rank entity

  stk::mesh::MeshIndex elem1_badMi = get_bulk().mesh_index(elem1);
  elem1_badMi.bucket_ordinal = 1;  // Only one element in Bucket
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadWrite>().entity_values(elem1_badMi));                       // Bad Bucket ordinal
  EXPECT_ANY_THROW(elemField.data().entity_values(elem1_badMi));  // Bad Bucket ordinal

  // Acquire FieldData before opening modification cycle
  {
    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().modification_begin();
    get_bulk().declare_node(100);
    get_bulk().modification_end();

    EXPECT_NO_THROW(elemFieldData.entity_values(elem1_mi));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.entity_values(elem1_mi));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.entity_values(node1_mi));       // Stale FieldData from before mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.entity_values(node1_mi));  // Stale FieldData from before mesh mod
  }

  // Acquire FieldData during modification cycle
  {
    get_bulk().modification_begin();

    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().declare_node(101);

    EXPECT_NO_THROW(elemFieldData.entity_values(elem1_mi));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.entity_values(elem1_mi));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.entity_values(node1_mi));       // Stale FieldData during mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.entity_values(node1_mi));  // Stale FieldData during mesh mod

    get_bulk().modification_end();
  }

  auto secondMesh = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_initial_bucket_capacity(1)
                                                          .set_maximum_bucket_capacity(1).create();
  stk::io::fill_mesh("generated:2x1x1", *secondMesh);  // Create two-element mesh
  const stk::mesh::Entity elem2 = secondMesh->get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::MeshIndex elem2_mi = secondMesh->mesh_index(elem2);
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadWrite>().entity_values(elem2_mi));                       // Entity from different mesh
  EXPECT_ANY_THROW(elemField.data().entity_values(elem2_mi));  // Entity from different mesh
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_consistencyCheck_fastMeshIndex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  stk::mesh::Field<int>& elemField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField1");
  stk::mesh::Field<int>& nodeField = get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodeField1");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(nodeField, get_meta().universal_part(), nullptr);
  create_single_element_mesh();

  const stk::mesh::FastMeshIndex elem1_badFmi1{1, 0};
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadWrite>().entity_values(elem1_badFmi1));                       // Bad Bucket ID
  EXPECT_ANY_THROW(elemField.data().entity_values(elem1_badFmi1));  // Bad Bucket ID

  const stk::mesh::FastMeshIndex elem1_badFmi2{0, 1};
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadWrite>().entity_values(elem1_badFmi2));                       // Bad Bucket ordinal
  EXPECT_ANY_THROW(elemField.data().entity_values(elem1_badFmi2));  // Bad Bucket ordinal

  // Acquire FieldData before opening modification cycle
  {
    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().modification_begin();
    get_bulk().declare_node(100);
    get_bulk().modification_end();

    const stk::mesh::FastMeshIndex elem1_fmi{0, 0};
    const stk::mesh::FastMeshIndex node1_fmi{1, 0};
    EXPECT_NO_THROW(elemFieldData.entity_values(elem1_fmi));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.entity_values(elem1_fmi));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.entity_values(node1_fmi));       // Stale FieldData from before mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.entity_values(node1_fmi));  // Stale FieldData from before mesh mod
  }

  // Acquire FieldData during modification cycle
  {
    get_bulk().modification_begin();

    auto elemFieldData = elemField.data<stk::mesh::ReadWrite>();
    auto constElemFieldData = elemField.data();
    auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite>();
    auto constNodeFieldData = nodeField.data();

    get_bulk().declare_node(101);

    const stk::mesh::FastMeshIndex elem1_fmi{0, 0};
    const stk::mesh::FastMeshIndex node1_fmi{1, 0};
    EXPECT_NO_THROW(elemFieldData.entity_values(elem1_fmi));        // Unmodified during mesh mod
    EXPECT_NO_THROW(constElemFieldData.entity_values(elem1_fmi));   // Unmodified during mesh mod
    EXPECT_ANY_THROW(nodeFieldData.entity_values(node1_fmi));       // Stale FieldData during mesh mod
    EXPECT_ANY_THROW(constNodeFieldData.entity_values(node1_fmi));  // Stale FieldData during mesh mod

    get_bulk().modification_end();
  }
}

//------------------------------------------------------------------------------
// This test cannot run normally because it triggers a device-side Kokkos::abort()
// which we can't trap.  Uncomment and run portions of it manually to confirm behavior.
//
// void device_consistency_check_entity(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& elemField,
//                                      stk::mesh::Field<int>& nodeField)
// {
//   {
//     const stk::mesh::Entity node8 = bulk.get_entity(stk::topology::NODE_RANK, 8);
//     auto elemFieldData = elemField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
//     auto constElemFieldData = elemField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
//
//     Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
//         auto entityValues = elemFieldData.entity_values(node8);           // Abort: Bad Bucket ordinal (wrong rank Entity)
//         auto constEntityValues = constElemFieldData.entity_values(node8); // Abort: Bad Bucket ordinal (wrong rank Entity)
//       }
//     );
//   }
//
//   // Acquire FieldData before opening modification cycle
//   {
//     const stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
//     const stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
//     auto elemFieldData = elemField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
//     auto constElemFieldData = elemField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
//     auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
//     auto constNodeFieldData = nodeField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
//
//     bulk.modification_begin();
//     bulk.declare_node(100);
//     bulk.modification_end();
//
//     Kokkos::parallel_for(1,
//       KOKKOS_LAMBDA(int) {
//         elemFieldData.entity_values(elem1);       // Unmodified during mesh mod
//         constElemFieldData.entity_values(elem1);  // Unmodified during mesh mod
//       }
//     );
//
//     Kokkos::parallel_for(1,
//       KOKKOS_LAMBDA(int) {
//         nodeFieldData.entity_values(node1);       // Abort: Stale FieldData from before mesh mod
//         constNodeFieldData.entity_values(node1);  // Abort: Stale FieldData from before mesh mod
//       }
//     );
//   }
//
//   // Acquire FieldData during modification cycle
//   {
//     bulk.modification_begin();
//
//     const stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
//     const stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
//     auto elemFieldData = elemField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
//     auto constElemFieldData = elemField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
//     auto nodeFieldData = nodeField.data<stk::mesh::ReadWrite, stk::ngp::DeviceSpace>();
//     auto constNodeFieldData = nodeField.data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();
//
//     bulk.declare_node(101);
//
//     Kokkos::parallel_for(1,
//       KOKKOS_LAMBDA(int) {
//         elemFieldData.entity_values(elem1);       // Unmodified during mesh mod
//         constElemFieldData.entity_values(elem1);  // Unmodified during mesh mod
//       }
//     );
//
//     Kokkos::parallel_for(1,
//       KOKKOS_LAMBDA(int) {
//         nodeFieldData.entity_values(node1);       // Abort: Stale FieldData during mesh mod
//         constNodeFieldData.entity_values(node1);  // Abort: Stale FieldData during mesh mod
//       }
//     );
//
//     bulk.modification_end();
//   }
// }
//
// TEST_F(FieldDataEntityAccess, device_consistencyCheck_entity)
// {
//   if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
//
//   setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many
//
//   stk::mesh::Field<int>& elemField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField1");
//   stk::mesh::Field<int>& nodeField = get_meta().declare_field<int>(stk::topology::NODE_RANK, "nodeField1");
//   stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
//   stk::mesh::put_field_on_mesh(nodeField, get_meta().universal_part(), nullptr);
//   create_single_element_mesh();
//
//   device_consistency_check_entity(get_bulk(), elemField, nodeField);
// }

//==============================================================================
TEST_F(FieldDataEntityAccess, host_scalarField_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);

  auto entityValues = field.data<stk::mesh::ReadWrite>().entity_values(node1);

  const stk::mesh::CopyIdx goodCopy = 0_copy;  // Only 1 copy
  const stk::mesh::CopyIdx badCopy  = 1_copy;

  const stk::mesh::ComponentIdx goodComponent = 0_comp;  // Only 1 components
  const stk::mesh::ComponentIdx badComponent  = 1_comp;

  EXPECT_NO_THROW(entityValues());                         // Normal scalar access
  EXPECT_NO_THROW(entityValues(goodComponent));            // In-bounds component
  EXPECT_NO_THROW(entityValues(goodCopy));                 // In-bounds copy
  EXPECT_NO_THROW(entityValues(goodCopy, goodComponent));  // In-bounds copy/component

  EXPECT_ANY_THROW(entityValues(badComponent));            // Out-of-bounds component
  EXPECT_ANY_THROW(entityValues(badComponent));            // Out-of-bounds copy
  EXPECT_ANY_THROW(entityValues(badCopy, goodComponent));  // Out-of-bounds copy/component
  EXPECT_ANY_THROW(entityValues(goodCopy, badComponent));  // Out-of-bounds copy/component
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponentField_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), 3, nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);

  auto entityValues = field.data<stk::mesh::ReadWrite>().entity_values(node1);

  const stk::mesh::CopyIdx goodCopy = 0_copy;  // Only 1 copy
  const stk::mesh::CopyIdx badCopy  = 1_copy;

  const stk::mesh::ComponentIdx goodComponent = 2_comp;  // Only 3 components
  const stk::mesh::ComponentIdx badComponent  = 3_comp;

  EXPECT_NO_THROW(entityValues(goodComponent));            // In-bounds component
  EXPECT_NO_THROW(entityValues(goodCopy, goodComponent));  // In-bounds copy/component

  EXPECT_ANY_THROW(entityValues());                        // Mistaken scalar access

  EXPECT_ANY_THROW(entityValues(goodCopy));                // Mistaken copy access

  EXPECT_ANY_THROW(entityValues(badComponent));            // Out-of-bounds component
  EXPECT_ANY_THROW(entityValues(badCopy, goodComponent));  // Out-of-bounds copy/component
  EXPECT_ANY_THROW(entityValues(goodCopy, badComponent));  // Out-of-bounds copy/component
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopyField_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), 1, 8, nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);

  auto entityValues = field.data<stk::mesh::ReadWrite>().entity_values(node1);

  const stk::mesh::CopyIdx goodCopy = 7_copy;  // Only 8 copies
  const stk::mesh::CopyIdx badCopy  = 8_copy;

  const stk::mesh::ComponentIdx goodComponent = 0_comp;  // Only 1 component
  const stk::mesh::ComponentIdx badComponent  = 1_comp;

  EXPECT_NO_THROW(entityValues(goodCopy));                 // In-bounds copy
  EXPECT_NO_THROW(entityValues(goodCopy, goodComponent));  // In-bounds copy/component

  EXPECT_ANY_THROW(entityValues());                        // Mistaken scalar access

  EXPECT_ANY_THROW(entityValues(badCopy));                 // Out-of-bounds copy
  EXPECT_ANY_THROW(entityValues(badCopy, goodComponent));  // Out-of-bounds copy/component
  EXPECT_ANY_THROW(entityValues(goodCopy, badComponent));  // Out-of-bounds copy/component

  EXPECT_ANY_THROW(entityValues(goodComponent));           // Mistaken component access
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopyMultiComponentField_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), 3, 8, nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);

  auto entityValues = field.data<stk::mesh::ReadWrite>().entity_values(node1);

  const stk::mesh::CopyIdx goodCopy = 7_copy;  // Only 8 copies
  const stk::mesh::CopyIdx badCopy  = 8_copy;

  const stk::mesh::ComponentIdx goodComponent = 2_comp;  // Only 3 components
  const stk::mesh::ComponentIdx badComponent  = 3_comp;

  EXPECT_NO_THROW(entityValues(goodCopy, goodComponent));  // In-bounds copy/component

  EXPECT_ANY_THROW(entityValues());                        // Mistaken scalar access

  EXPECT_ANY_THROW(entityValues(goodCopy));                // Mistaken copy access

  EXPECT_ANY_THROW(entityValues(goodComponent));           // Mistaken component access

  EXPECT_ANY_THROW(entityValues(badCopy, goodComponent));  // Out-of-bounds copy/component
  EXPECT_ANY_THROW(entityValues(goodCopy, badComponent));  // Out-of-bounds copy/component
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiScalar_boundsCheck)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), 3, 8, nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);

  auto entityValues = field.data<stk::mesh::ReadWrite>().entity_values(node1);

  const stk::mesh::ScalarIdx goodScalar = 23_scalar;  // Only 24 scalars
  const stk::mesh::ScalarIdx badScalar  = 24_scalar;

  EXPECT_NO_THROW(entityValues(goodScalar));                 // In-bounds scalar
  EXPECT_ANY_THROW(entityValues(badScalar));                 // Out-of-bounds scalar
}

#endif  // STK_FIELD_BOUNDS_CHECK

//==============================================================================

void change_node_val_on_device(stk::mesh::Entity node,
                               const stk::mesh::FieldBase* coordField,
                               double newVal)
{
  auto deviceCoordData = coordField->data<double,stk::mesh::ReadWrite,stk::ngp::DeviceSpace>();
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int& /*idx*/)
  {
    auto deviceCoordVals = deviceCoordData.entity_values(node);
    deviceCoordVals(0_comp) = newVal;
  });
}

TEST(TestFieldData, syncCount_emptyModCycle)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }
  auto mesh = stk::mesh::MeshBuilder(comm).create();
  stk::io::fill_mesh("generated:1x1x1", *mesh);

  auto coordField = mesh->mesh_meta_data().coordinate_field();

  stk::mesh::Entity node1 = mesh->get_entity(stk::topology::NODE_RANK, 1);
  ASSERT_TRUE(mesh->is_valid(node1));
  {
    auto coordFieldData = coordField->data<double>();
    auto coordVals = coordFieldData.entity_values(node1);
    EXPECT_DOUBLE_EQ(0.0, coordVals(0_comp));
  }

  mesh->modification_begin();
  mesh->declare_node(9, stk::mesh::PartVector{});
  mesh->modification_end();

  {
    auto coordFieldData = coordField->data<double,stk::mesh::ReadWrite>();
    auto coordVals = coordFieldData.entity_values(node1);
    coordVals(0_comp) = 99.9;
  }

  const unsigned syncCount = mesh->synchronized_count();

  mesh->modification_begin();
  {
    auto coordFieldData = coordField->data<double>();
    auto coordVals = coordFieldData.entity_values(node1);
    EXPECT_DOUBLE_EQ(99.9, coordVals(0_comp));
  }

  unsigned tmpSyncCount = mesh->synchronized_count();
  EXPECT_EQ((syncCount+1), tmpSyncCount);

  mesh->modification_end();

  EXPECT_EQ(mesh->synchronized_count(), syncCount);

  double newValue = 999;
  {
    stk::mesh::Entity node9 = mesh->get_entity(stk::topology::NODE_RANK,9);
    change_node_val_on_device(node9, coordField, newValue);
  }

  {
    auto coordFieldData = coordField->data<double>();
    stk::mesh::Entity node9 = mesh->get_entity(stk::topology::NODE_RANK,9);
    auto coordVals = coordFieldData.entity_values(node9);
    EXPECT_DOUBLE_EQ(newValue, coordVals(0_comp));
  }
}

} //namespace <anonymous>
