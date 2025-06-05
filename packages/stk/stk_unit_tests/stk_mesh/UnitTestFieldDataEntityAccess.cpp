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
#include <stk_mesh/base/FieldBytes.hpp>
#include <stk_mesh/base/ConstFieldBytes.hpp>
#include <stk_mesh/base/FieldDataBase.hpp>
#include <stk_mesh/base/FieldIndexTypes.hpp>

namespace {

class FieldDataEntityAccess : public FieldDataAccessFixture {};

//==============================================================================
TEST_F(FieldDataEntityAccess, inconsistentTemplateParameters)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::HostMemSpace, stk::mesh::Layout::Left>()));  // Correct
  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>()));  // Correct

  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadOnly,  stk::ngp::HostMemSpace, stk::mesh::Layout::Left>()));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>()));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::HostMemSpace, stk::mesh::Layout::Right>()));  // Wrong layout
  EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostMemSpace, stk::mesh::Layout::Right>()));  // Wrong layout

  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::MemSpace, stk::mesh::Layout::Left>()));  // Correct
  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace, stk::mesh::Layout::Left>()));  // Correct

  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadOnly,  stk::ngp::MemSpace, stk::mesh::Layout::Left>()));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::MemSpace, stk::mesh::Layout::Left>()));  // Wrong datatype
  //EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::MemSpace, stk::mesh::Layout::Right>()));  // Trapped by static_assert()
  //EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace, stk::mesh::Layout::Right>()));  // Trapped by static_assert()
}

TEST_F(FieldDataEntityAccess, inconsistentTemplateParameters_async)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);

  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::HostMemSpace,
                                  stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Correct
  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostMemSpace,
                                  stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Correct

  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadOnly, stk::ngp::HostMemSpace,
                                   stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::HostMemSpace,
                                   stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::HostMemSpace,
                                   stk::mesh::Layout::Right>(stk::ngp::ExecSpace())));  // Wrong layout
  EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostMemSpace,
                                   stk::mesh::Layout::Right>(stk::ngp::ExecSpace())));  // Wrong layout

  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::MemSpace,
                                  stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Correct
  EXPECT_NO_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace,
                                  stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Correct


  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadOnly,  stk::ngp::MemSpace,
                                   stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Wrong datatype
  EXPECT_ANY_THROW((fieldBase.data<double, stk::mesh::ReadWrite, stk::ngp::MemSpace,
                                   stk::mesh::Layout::Left>(stk::ngp::ExecSpace())));  // Wrong datatype
  // EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadOnly,  stk::ngp::MemSpace,
  //                                  stk::mesh::Layout::Right>(stk::ngp::ExecSpace())));  // Trapped by static_assert()
  // EXPECT_ANY_THROW((fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace,
  //                                  stk::mesh::Layout::Right>(stk::ngp::ExecSpace())));  // Trapped by static_assert()
}

//==============================================================================
TEST_F(FieldDataEntityAccess, host_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      entityValues() = ++value;
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      EXPECT_EQ(constEntityValues(), ++value);
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      entityValuesBase() = ++value*10;
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      EXPECT_EQ(constEntityValuesBase(), ++value*10);
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int, stk::mesh::Layout::Left>
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      entityValues() = ++value;
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      EXPECT_EQ(constEntityValues(), ++value);
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      entityValuesBase() = ++value*10;
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      EXPECT_EQ(constEntityValuesBase(), ++value*10);
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_scalar_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      entityValues() = ++value;
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      EXPECT_EQ(*constEntityPtr, ++value);
    }
  }
}


//==============================================================================
TEST_F(FieldDataEntityAccess, host_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  stk::mesh::Field<int>& field = *m_field;
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        entityValues(component) = ++value;
      }
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
        EXPECT_EQ(constEntityValues(component), ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValuesBase.components()) {
        entityValuesBase(component) = ++value*10;
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::ComponentIdx component : constEntityValuesBase.components()) {
        EXPECT_EQ(constEntityValuesBase(component), ++value*10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int, stk::mesh::Layout::Left>
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        entityValues(component) = ++value;
      }
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
        EXPECT_EQ(constEntityValues(component), ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValuesBase.components()) {
        entityValuesBase(component) = ++value*10;
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::ComponentIdx component : constEntityValuesBase.components()) {
        EXPECT_EQ(constEntityValuesBase(component), ++value*10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiComponent_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        entityValues(component) = ++value;
      }
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      const int componentStride = constEntityValues.component_stride();
      for (int component = 0; component < constEntityValues.num_components(); ++component) {
        EXPECT_EQ(constEntityPtr[component*componentStride], ++value);
      }
    }
  }
}


//==============================================================================
TEST_F(FieldDataEntityAccess, host_multiCopy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();

  stk::mesh::Field<int>& field = *m_field;
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        entityValues(copy) = ++value;
      }
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        EXPECT_EQ(constEntityValues(copy), ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValuesBase.copies()) {
        entityValuesBase(copy) = ++value*10;
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValuesBase.copies()) {
        EXPECT_EQ(constEntityValuesBase(copy), ++value*10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int, stk::mesh::Layout::Left>
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        entityValues(copy) = ++value;
      }
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        EXPECT_EQ(constEntityValues(copy), ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValuesBase.copies()) {
        entityValuesBase(copy) = ++value*10;
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValuesBase.copies()) {
        EXPECT_EQ(constEntityValuesBase(copy), ++value*10);
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();

  stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        entityValues(copy) = ++value;
      }
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      const int copyStride = constEntityValues.copy_stride();
      for (int copy = 0; copy < constEntityValues.num_copies(); ++copy) {
        EXPECT_EQ(constEntityPtr[copy*copyStride], ++value);
      }
    }
  }
}


//==============================================================================
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  stk::mesh::Field<int>& field = *m_field;
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(copy, component) = ++value;
        }
      }
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          EXPECT_EQ(constEntityValues(copy, component), ++value);
        }
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValuesBase.copies()) {
        for (stk::mesh::ComponentIdx component : entityValuesBase.components()) {
          entityValuesBase(copy, component) = ++value*10;
        }
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValuesBase.copies()) {
        for (stk::mesh::ComponentIdx component : constEntityValuesBase.components()) {
          EXPECT_EQ(constEntityValuesBase(copy, component), ++value*10);
        }
      }
    }
  }
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
  auto fieldData = field.data();
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
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
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
  auto fieldDataBase = fieldBase.data<int>();
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
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly>();
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
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();

  const stk::mesh::Field<int, stk::mesh::Layout::Left>& field = *m_leftField;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_leftField);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int, stk::mesh::Layout::Left>
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(copy, component) = ++value;
        }
      }
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          EXPECT_EQ(constEntityValues(copy, component), ++value);
        }
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValuesBase.copies()) {
        for (stk::mesh::ComponentIdx component : entityValuesBase.components()) {
          entityValuesBase(copy, component) = ++value*10;
        }
      }
    }
  }

  value = 0;
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::HostMemSpace, stk::mesh::Layout::Left>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValuesBase.copies()) {
        for (stk::mesh::ComponentIdx component : constEntityValuesBase.components()) {
          EXPECT_EQ(constEntityValuesBase(copy, component), ++value*10);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_multiCopy_multiComponent_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  int value = 0;
  auto fieldData = field.data();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(copy, component) = ++value;
        }
      }
    }
  }

  value = 0;
  auto constFieldData = field.data<stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      const int copyStride = constEntityValues.copy_stride();
      const int componentStride = constEntityValues.component_stride();
      for (int copy = 0; copy < constEntityValues.num_copies(); ++copy) {
        for (int component = 0; component < constEntityValues.num_components(); ++component) {
          EXPECT_EQ(constEntityPtr[copy*copyStride + component*componentStride], ++value);
        }
      }
    }
  }
}


//==============================================================================
void test_device_scalar_entity_values(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(field);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Write and read values from Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      entityValues() = (entity.bucket_id*10 + entity.bucket_ord);
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      NGP_EXPECT_EQ(constEntityValues(), static_cast<int>(entity.bucket_id*10 + entity.bucket_ord));
    }
  );

  // Write and read values from FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      entityValuesBase() = (entity.bucket_id*10 + entity.bucket_ord);
    }
  );

  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      NGP_EXPECT_EQ(constEntityValuesBase(), static_cast<int>(entity.bucket_id*10 + entity.bucket_ord));
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();
  test_device_scalar_entity_values(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_scalar_entity_values_pointer(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      entityValues() = (entity.bucket_id*10 + entity.bucket_ord);
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
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
  test_device_scalar_entity_values_pointer(get_bulk(), *m_field);
}


//==============================================================================
void test_device_multi_component_entity_values(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(field);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Write value to Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        entityValues(component) = (entity.bucket_id*100 + entity.bucket_ord*10 + component);
      }
    }
  );

  // Read const value from Field<int>
  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
        NGP_EXPECT_EQ(constEntityValues(component),
                      static_cast<int>(entity.bucket_id*100 + entity.bucket_ord*10 + component));
      }
    }
  );

  // Write value to FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValuesBase.components()) {
        entityValuesBase(component) = (entity.bucket_id*1000 + entity.bucket_ord*100 + component*10);
      }
    }
  );

  // Read const value from FieldBase
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::ComponentIdx component : constEntityValuesBase.components()) {
        NGP_EXPECT_EQ(constEntityValuesBase(component),
                      static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + component*10));
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();
  test_device_multi_component_entity_values(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_multi_component_entity_values_async(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(field);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Use the default device exec space to trigger the async code path.  No proper streams exercised here.

  // Write value to Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>(stk::ngp::ExecSpace());
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        entityValues(component) = (entity.bucket_id*100 + entity.bucket_ord*10 + component);
      }
    }
  );

  // Read const value from Field<int>
  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>(stk::ngp::ExecSpace());
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
        NGP_EXPECT_EQ(constEntityValues(component),
                      static_cast<int>(entity.bucket_id*100 + entity.bucket_ord*10 + component));
      }
    }
  );

  // Write value to FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace>(stk::ngp::ExecSpace());
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValuesBase.components()) {
        entityValuesBase(component) = (entity.bucket_id*1000 + entity.bucket_ord*100 + component*10);
      }
    }
  );

  // Read const value from FieldBase
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::MemSpace>(stk::ngp::ExecSpace());
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::ComponentIdx component : constEntityValuesBase.components()) {
        NGP_EXPECT_EQ(constEntityValuesBase(component),
                      static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + component*10));
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiComponent_async)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();
  test_device_multi_component_entity_values_async(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_multi_component_entity_values_pointer(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::ComponentIdx component : entityValues.components()) {
        entityValues(component) = (entity.bucket_id*100 + entity.bucket_ord*10 + component);
      }
    }
  );

  // Read const value from Field<int>
  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
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
  test_device_multi_component_entity_values_pointer(get_bulk(), *m_field);
}


//==============================================================================
void test_device_multi_copy_entity_values(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(field);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Write value to Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        entityValues(copy) = (entity.bucket_id*100 + entity.bucket_ord*10 + copy);
      }
    }
  );

  // Read const value from Field<int>
  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        NGP_EXPECT_EQ(constEntityValues(copy),
                      static_cast<int>(entity.bucket_id*100 + entity.bucket_ord*10 + copy));
      }
    }
  );

  // Write value to FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValuesBase.copies()) {
        entityValuesBase(copy) = (entity.bucket_id*1000 + entity.bucket_ord*100 + copy*10);
      }
    }
  );

  // Read const value from FieldBase
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValuesBase.copies()) {
        NGP_EXPECT_EQ(constEntityValuesBase(copy),
                      static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + copy*10));
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiCopy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();
  test_device_multi_copy_entity_values(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_multi_copy_entity_values_pointer(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Write the values normally and read them through a raw pointer to make sure
  // indexing is consistent between the two APIs
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        entityValues(copy) = (entity.bucket_id*100 + entity.bucket_ord*10 + copy);
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
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
  test_device_multi_copy_entity_values_pointer(get_bulk(), *m_field);
}


//==============================================================================
void test_device_multi_copy_multi_component_entity_values(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(field);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Write and read values from Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(copy, component) = (entity.bucket_id*1000 + entity.bucket_ord*100 + static_cast<int>(copy)*10 +
                                           static_cast<int>(component));
        }
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
        for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
          NGP_EXPECT_EQ(constEntityValues(copy, component),
                        static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + static_cast<int>(copy)*10 +
                                         static_cast<int>(component)));
        }
      }
    }
  );

  // Write and read values from FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValuesBase.copies()) {
        for (stk::mesh::ComponentIdx component : entityValuesBase.components()) {
          entityValuesBase(copy, component) = (entity.bucket_id*10000 + entity.bucket_ord*1000 + static_cast<int>(copy)*100 +
                                               static_cast<int>(component)*10);
        }
      }
    }
  );

  // Read const value from FieldBase
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy : constEntityValuesBase.copies()) {
        for (stk::mesh::ComponentIdx component : constEntityValuesBase.components()) {
          NGP_EXPECT_EQ(constEntityValuesBase(copy, component),
                        static_cast<int>(entity.bucket_id*10000 + entity.bucket_ord*1000 + static_cast<int>(copy)*100 +
                                         static_cast<int>(component)*10));
        }
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_multi_copy_multi_component_entity_values(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_multi_copy_multi_component_entity_values_pointer(stk::mesh::BulkData& bulk,
                                                                  stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Write and read values from Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(copy, component) = (entity.bucket_id*1000 + entity.bucket_ord*100 + static_cast<int>(copy)*10 +
                                           static_cast<int>(component));
        }
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      const int* constEntityPtr = constEntityValues.pointer();
      const int copyStride = constEntityValues.copy_stride();
      const int componentStride = constEntityValues.component_stride();
      for (int copy = 0; copy < constEntityValues.num_copies(); ++copy) {
        for (int component = 0; component < constEntityValues.num_components(); ++component) {
          NGP_EXPECT_EQ(constEntityPtr[copy*copyStride + component*componentStride],
                        static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + static_cast<int>(copy)*10 +
                                         static_cast<int>(component)));
        }
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_multiComponent_pointer)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_multi_copy_multi_component_entity_values_pointer(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_multi_copy_multi_component_entity_values_traditional_for_loop(stk::mesh::BulkData& bulk,
                                                                               stk::mesh::Field<int>& field)
{
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(field);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  // Write and read values from Field<int>
  auto fieldData = field.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValues = fieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < entityValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < entityValues.num_components(); ++component) {
          entityValues(copy, component) = (entity.bucket_id*1000 + entity.bucket_ord*100 + static_cast<int>(copy)*10 +
                                           static_cast<int>(component));
        }
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValues = constFieldData.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < constEntityValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constEntityValues.num_components(); ++component) {
          NGP_EXPECT_EQ(constEntityValues(copy, component),
                        static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + static_cast<int>(copy)*10 +
                                         static_cast<int>(component)));
        }
      }
    }
  );

  // Write and read values from FieldBase
  auto fieldDataBase = fieldBase.data<int, stk::mesh::ReadWrite, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto entityValuesBase = fieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < entityValuesBase.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < entityValuesBase.num_components(); ++component) {
          entityValuesBase(copy, component) = (entity.bucket_id*10000 + entity.bucket_ord*1000 + static_cast<int>(copy)*100 +
                                               static_cast<int>(component)*10);
        }
      }
    }
  );

  // Read const value from FieldBase
  auto constFieldDataBase = fieldBase.data<int, stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto constEntityValuesBase = constFieldDataBase.entity_values(entity);
      for (stk::mesh::CopyIdx copy(0); copy < constEntityValuesBase.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constEntityValuesBase.num_components(); ++component) {
          NGP_EXPECT_EQ(constEntityValuesBase(copy, component),
                        static_cast<int>(entity.bucket_id*10000 + entity.bucket_ord*1000 + static_cast<int>(copy)*100 +
                                         static_cast<int>(component)*10));
        }
      }
    }
  );
}

NGP_TEST_F(FieldDataEntityAccess, device_multiCopy_multiComponent_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_device_multi_copy_multi_component_entity_values_traditional_for_loop(get_bulk(), *m_field);
}


//==============================================================================
template <typename FieldType>
void test_host_to_device_scalar_entity_values(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>(); // Create early so next call is sync instead of update

  {
    auto fieldData = field.data();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto entityValues = fieldData.entity_values(entity);
        entityValues() = nodeId;
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, entity));
        auto constEntityValues = constFieldData.entity_values(entity);
        NGP_EXPECT_EQ(constEntityValues(), nodeId);
      }
    );
  }
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_field();
  test_host_to_device_scalar_entity_values(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_scalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_left_field();
  test_host_to_device_scalar_entity_values(get_bulk(), *m_leftField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_multi_component_entity_values(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>(); // Create early so next call is sync instead of update

  {
    auto fieldData = field.data();
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
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, entity));
        auto constEntityValues = constFieldData.entity_values(entity);
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
  test_host_to_device_multi_component_entity_values(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_left_field();
  test_host_to_device_multi_component_entity_values(get_bulk(), *m_leftField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_multi_copy_entity_values(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>(); // Create early so next call is sync instead of update

  {
    auto fieldData = field.data();
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
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, entity));
        auto constEntityValues = constFieldData.entity_values(entity);
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
  test_host_to_device_multi_copy_entity_values(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiCopy_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_left_field();
  test_host_to_device_multi_copy_entity_values(get_bulk(), *m_leftField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_multi_copy_multi_component_entity_values(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>(); // Create early so next call is sync instead of update

  {
    auto fieldData = field.data();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            entityValues(copy, component) = nodeId*100 + 10*copy + static_cast<int>(component);
          }
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, entity));
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
            NGP_EXPECT_EQ(constEntityValues(copy, component), (nodeId*100 + 10*copy + static_cast<int>(component)));
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
  test_host_to_device_multi_copy_multi_component_entity_values(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedHostToDevice_multiCopy_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_host_to_device_multi_copy_multi_component_entity_values(get_bulk(), *m_leftField);
}


//==============================================================================
template <typename FieldType>
void test_device_to_host_scalar_entity_values(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, entity));
        auto entityValues = fieldData.entity_values(entity);
        entityValues() = nodeId;
      }
    );
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly>();
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
  test_device_to_host_scalar_entity_values(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_scalar_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_scalar_left_field();
  test_device_to_host_scalar_entity_values(get_bulk(), *m_leftField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_multi_component_entity_values(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, entity));
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          entityValues(component) = nodeId*10 + component;
        }
      }
    );
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly>();
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
  test_device_to_host_multi_component_entity_values(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_left_field();
  test_device_to_host_multi_component_entity_values(get_bulk(), *m_leftField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_multi_copy_entity_values(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, entity));
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          entityValues(copy) = nodeId*100 + 10*copy;
        }
      }
    );
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly>();
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
  test_device_to_host_multi_copy_entity_values(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiCopy_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_left_field();
  test_device_to_host_multi_copy_entity_values(get_bulk(), *m_leftField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_multi_copy_multi_component_entity_values(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);

  {
    auto fieldData = field.template data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
        int nodeId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::NODE_RANK, entity));
        auto entityValues = fieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            entityValues(copy, component) = nodeId*100 + 10*copy + static_cast<int>(component);
          }
        }
      }
    );
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly>();
    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity entity : *bucket) {
        int nodeId = bulk.identifier(entity);
        auto constEntityValues = constFieldData.entity_values(entity);
        for (stk::mesh::CopyIdx copy : constEntityValues.copies()) {
          for (stk::mesh::ComponentIdx component : constEntityValues.components()) {
            NGP_EXPECT_EQ(constEntityValues(copy, component), (nodeId*100 + 10*copy + static_cast<int>(component)));
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
  test_device_to_host_multi_copy_multi_component_entity_values(get_bulk(), *m_field);
}

NGP_TEST_F(FieldDataEntityAccess, mixedDeviceToHost_multiCopy_multiComponent_layoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_left_field();
  test_device_to_host_multi_copy_multi_component_entity_values(get_bulk(), *m_leftField);
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

  auto fieldData = field.data<stk::mesh::ReadOnly>();
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
  auto fieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();

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
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  stk::mesh::Field<int>& elemField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField1");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  create_single_element_mesh();
  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);

  EXPECT_ANY_THROW(elemField.data().entity_values(node1));                       // Wrong rank entity
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadOnly>().entity_values(node1));  // Wrong rank entity

  auto fieldData = elemField.data();
  auto constFieldData = elemField.data<stk::mesh::ReadOnly>();
  create_node(100);
  EXPECT_ANY_THROW(fieldData.entity_values(elem1));       // Stale FieldData from before mesh mod
  EXPECT_ANY_THROW(constFieldData.entity_values(elem1));  // Stale FieldData from before mesh mod
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_consistencyCheck_meshIndex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  stk::mesh::Field<int>& elemField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField1");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  create_single_element_mesh();
  const stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);

  const stk::mesh::MeshIndex node1_mi = get_bulk().mesh_index(node1);
  EXPECT_ANY_THROW(elemField.data().entity_values(node1_mi));                       // Wrong rank entity
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadOnly>().entity_values(node1_mi));  // Wrong rank entity

  stk::mesh::MeshIndex elem1_badMi = get_bulk().mesh_index(elem1);
  elem1_badMi.bucket_ordinal = 1;  // Only one element in Bucket
  EXPECT_ANY_THROW(elemField.data().entity_values(elem1_badMi));                       // Bad Bucket ordinal
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadOnly>().entity_values(elem1_badMi));  // Bad Bucket ordinal

  auto fieldData = elemField.data();
  auto constFieldData = elemField.data<stk::mesh::ReadOnly>();
  create_node(100);
  const stk::mesh::MeshIndex elem1_mi = get_bulk().mesh_index(elem1);
  EXPECT_ANY_THROW(fieldData.entity_values(elem1_mi));       // Stale FieldData from before mesh mod
  EXPECT_ANY_THROW(constFieldData.entity_values(elem1_mi));  // Stale FieldData from before mesh mod

  auto secondMesh = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_initial_bucket_capacity(1)
                                                          .set_maximum_bucket_capacity(1).create();
  stk::io::fill_mesh("generated:2x1x1", *secondMesh);  // Create two-element mesh
  const stk::mesh::Entity elem2 = secondMesh->get_entity(stk::topology::ELEM_RANK, 2);
  const stk::mesh::MeshIndex elem2_mi = secondMesh->mesh_index(elem2);
  EXPECT_ANY_THROW(elemField.data().entity_values(elem2_mi));                       // Entity from different mesh
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadOnly>().entity_values(elem2_mi));  // Entity from different mesh
}

//------------------------------------------------------------------------------
TEST_F(FieldDataEntityAccess, host_consistencyCheck_fastMeshIndex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  stk::mesh::Field<int>& elemField = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField1");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  create_single_element_mesh();

  const stk::mesh::FastMeshIndex elem1_badFmi1{1, 0};
  EXPECT_ANY_THROW(elemField.data().entity_values(elem1_badFmi1));                       // Bad Bucket ID
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadOnly>().entity_values(elem1_badFmi1));  // Bad Bucket ID

  const stk::mesh::FastMeshIndex elem1_badFmi2{0, 1};
  EXPECT_ANY_THROW(elemField.data().entity_values(elem1_badFmi2));                       // Bad Bucket ordinal
  EXPECT_ANY_THROW(elemField.data<stk::mesh::ReadOnly>().entity_values(elem1_badFmi2));  // Bad Bucket ordinal

  auto fieldData = elemField.data();
  auto constFieldData = elemField.data<stk::mesh::ReadOnly>();
  create_node(100);
  const stk::mesh::MeshIndex elem1_fmi{0, 0};
  EXPECT_ANY_THROW(fieldData.entity_values(elem1_fmi));       // Stale FieldData from before mesh mod
  EXPECT_ANY_THROW(constFieldData.entity_values(elem1_fmi));  // Stale FieldData from before mesh mod
}

//------------------------------------------------------------------------------
// This test cannot run normally because it triggers a device-side Kokkos::abort()
// which we can't trap.  Uncomment and run manually to confirm behavior.
//
// void device_consistency_check_entity(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& elemField)
// {
//   {
//     const stk::mesh::Entity node8 = bulk.get_entity(stk::topology::NODE_RANK, 8);
//     auto fieldData = elemField.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
//     auto constFieldData = elemField.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
//
//     Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
//         auto entityValues = fieldData.entity_values(node8);           // Abort: Bad Bucket ordinal (wrong rank Entity)
//         auto constEntityValues = constFieldData.entity_values(node8); // Abort: Bad Bucket ordinal (wrong rank Entity)
//       }
//     );
//   }
//
//   {
//     const stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
//     auto fieldData = elemField.data<stk::mesh::ReadWrite, stk::ngp::MemSpace>();
//     auto constFieldData = elemField.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
//
//     bulk.modification_begin();
//     bulk.declare_node(100);
//     bulk.modification_end();
//
//     Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
//         auto entityValues = fieldData.entity_values(elem1);           // Abort: Stale FieldData from before mesh mod
//         auto constEntityValues = constFieldData.entity_values(elem1); // Abort: Stale FieldData from before mesh mod
//       }
//     );
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
//   stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
//   create_single_element_mesh();
//
//   device_consistency_check_entity(get_bulk(), elemField);
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

  auto entityValues = field.data().entity_values(node1);

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

  auto entityValues = field.data().entity_values(node1);

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

  auto entityValues = field.data().entity_values(node1);

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
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::Field<int>& field = get_meta().declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, get_meta().universal_part(), 3, 8, nullptr);
  const stk::mesh::Entity node1 = create_node(1);
  create_node(2);

  auto entityValues = field.data().entity_values(node1);

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

#endif  // STK_FIELD_BOUNDS_CHECK

//==============================================================================

} //namespace <anonymous>
