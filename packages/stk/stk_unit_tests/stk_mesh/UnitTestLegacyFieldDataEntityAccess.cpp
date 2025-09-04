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
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>

namespace {

#ifndef STK_UNIFIED_MEMORY

class LegacyFieldDataEntityAccess : public FieldDataAccessFixture {};

//==============================================================================
TEST_F(LegacyFieldDataEntityAccess, host_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();

  const stk::mesh::Field<int>& field = *m_field;
  const stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);

  // Write and read values from Field<int>
  int value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      int* entityValues = stk::mesh::field_data(field, entity);
      *entityValues = ++value;
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      const int* constEntityValues = stk::mesh::field_data(field, entity);
      EXPECT_EQ(*constEntityValues, ++value);
    }
  }

  // Write and read values from FieldBase
  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      int* entityValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, entity));
      *entityValuesBase = ++value*10;
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      const int* constEntityValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, entity));
      EXPECT_EQ(*constEntityValuesBase, ++value*10);
    }
  }
}

//==============================================================================
TEST_F(LegacyFieldDataEntityAccess, host_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();

  stk::mesh::Field<int>& field = *m_field;
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);
  const int numComponents = stk::mesh::field_extent0_per_entity(field, *buckets.front());

  // Write and read values from Field<int>
  int value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      int* entityValues = stk::mesh::field_data(field, entity);
      for (int component = 0; component < numComponents; ++component) {
        entityValues[component] = ++value;
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      const int* constEntityValues = stk::mesh::field_data(field, entity);
      for (int component = 0; component < numComponents; ++component) {
        EXPECT_EQ(constEntityValues[component], ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      int* entityValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, entity));
      for (int component = 0; component < numComponents; ++component) {
        entityValuesBase[component] = ++value*10;
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      const int* constEntityValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, entity));
      for (int component = 0; component < numComponents; ++component) {
        EXPECT_EQ(constEntityValuesBase[component], ++value*10);
      }
    }
  }
}

//==============================================================================
TEST_F(LegacyFieldDataEntityAccess, host_multiCopy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();

  stk::mesh::Field<int>& field = *m_field;
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);
  const int numCopies = stk::mesh::field_extent1_per_entity(field, *buckets.front());

  // Write and read values from Field<int>
  int value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      int* entityValues = stk::mesh::field_data(field, entity);
      for (int copy = 0; copy < numCopies; ++copy) {
        entityValues[copy] = ++value;
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      const int* constEntityValues = stk::mesh::field_data(field, entity);
      for (int copy = 0; copy < numCopies; ++copy) {
        EXPECT_EQ(constEntityValues[copy], ++value);
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      int* entityValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, entity));
      for (int copy = 0; copy < numCopies; ++copy) {
        entityValuesBase[copy] = ++value*10;
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      const int* constEntityValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, entity));
      for (int copy = 0; copy < numCopies; ++copy) {
        EXPECT_EQ(constEntityValuesBase[copy], ++value*10);
      }
    }
  }
}

//==============================================================================
TEST_F(LegacyFieldDataEntityAccess, host_multiCopy_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();

  stk::mesh::Field<int>& field = *m_field;
  stk::mesh::FieldBase& fieldBase = static_cast<stk::mesh::FieldBase&>(*m_field);
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::NODE_RANK);
  const int numCopies = stk::mesh::field_extent1_per_entity(field, *buckets.front());
  const int numComponents = stk::mesh::field_extent0_per_entity(field, *buckets.front());

  // Write and read values from Field<int>
  int value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      int* entityValues = stk::mesh::field_data(field, entity);
      for (int copy = 0; copy < numCopies; ++copy) {
        for (int component = 0; component < numComponents; ++component) {
          entityValues[copy*numComponents + component] = ++value;
        }
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      const int* constEntityValues = stk::mesh::field_data(field, entity);
      for (int copy = 0; copy < numCopies; ++copy) {
        for (int component = 0; component < numComponents; ++component) {
          EXPECT_EQ(constEntityValues[copy*numComponents + component], ++value);
        }
      }
    }
  }

  // Write and read values from FieldBase
  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      int* entityValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, entity));
      for (int copy = 0; copy < numCopies; ++copy) {
        for (int component = 0; component < numComponents; ++component) {
          entityValuesBase[copy*numComponents + component] = ++value*10;
        }
      }
    }
  }

  value = 0;
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity entity : *bucket) {
      const int* constEntityValuesBase = static_cast<int*>(stk::mesh::field_data(fieldBase, entity));
      for (int copy = 0; copy < numCopies; ++copy) {
        for (int component = 0; component < numComponents; ++component) {
          EXPECT_EQ(constEntityValuesBase[copy*numComponents + component], ++value*10);
        }
      }
    }
  }
}


//==============================================================================
void test_legacy_device_scalar_entity_values(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(field);

  // Write value to Field<int>
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      ngpField(entity, 0) = (entity.bucket_id*10 + entity.bucket_ord);
    }
  );

  // Read value from Field<int>
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      NGP_EXPECT_EQ(ngpField(entity, 0), static_cast<int>(entity.bucket_id*10 + entity.bucket_ord));
    }
  );
}

NGP_TEST_F(LegacyFieldDataEntityAccess, device_scalar)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();
  build_mesh_with_scalar_field();
  test_legacy_device_scalar_entity_values(get_bulk(), *m_field);
}


//==============================================================================
void test_legacy_device_multi_component_entity_values(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(field);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  const int numComponents = stk::mesh::field_extent0_per_entity(field, *buckets.front());

  // Write and read values from Field<int>
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      for (int component = 0; component < numComponents; ++component) {
        ngpField(entity, component) = (entity.bucket_id*100 + entity.bucket_ord*10 + component);
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      for (int component = 0; component < numComponents; ++component) {
        NGP_EXPECT_EQ(ngpField(entity, component), static_cast<int>(entity.bucket_id*100 + entity.bucket_ord*10 + component));
      }
    }
  );
}

NGP_TEST_F(LegacyFieldDataEntityAccess, device_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_component_field();
  test_legacy_device_multi_component_entity_values(get_bulk(), *m_field);
}


//==============================================================================
void test_legacy_device_multi_copy_entity_values(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(field);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  const int numCopies = stk::mesh::field_extent1_per_entity(field, *buckets.front());

 // Write and read values from Field<int>
 stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
     for (int copy = 0; copy < numCopies; ++copy) {
       ngpField(entity, copy) = (entity.bucket_id*100 + entity.bucket_ord*10 + copy);
     }
   }
 );

 // Read value from Field<int>
 stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
   KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
     for (int copy = 0; copy < numCopies; ++copy) {
       NGP_EXPECT_EQ(ngpField(entity, copy), static_cast<int>(entity.bucket_id*100 + entity.bucket_ord*10 + copy));
     }
   }
 );
}

NGP_TEST_F(LegacyFieldDataEntityAccess, device_multiCopy)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_field();
  test_legacy_device_multi_copy_entity_values(get_bulk(), *m_field);
}


//==============================================================================
void test_legacy_device_multi_copy_multi_component_entity_values(stk::mesh::BulkData& bulk,
                                                                 stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(field);
  const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::NODE_RANK);
  const int numCopies = stk::mesh::field_extent1_per_entity(field, *buckets.front());
  const int numComponents = stk::mesh::field_extent0_per_entity(field, *buckets.front());

  // Write value to Field<int>
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      for (int copy = 0; copy < numCopies; ++copy) {
        for (int component = 0; component < numComponents; ++component) {
          ngpField(entity, copy*numComponents + component) =
              (entity.bucket_id*1000 + entity.bucket_ord*100 + copy*10 + component);
        }
      }
    }
  );

  // Read value from Field<int>
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      for (int copy = 0; copy < numCopies; ++copy) {
        for (int component = 0; component < numComponents; ++component) {
          NGP_EXPECT_EQ(ngpField(entity, copy*numComponents + component),
                        static_cast<int>(entity.bucket_id*1000 + entity.bucket_ord*100 + copy*10 + component));
        }
      }
    }
  );
}

NGP_TEST_F(LegacyFieldDataEntityAccess, device_multiCopy_multiComponent)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  build_mesh_with_multi_copy_multi_component_field();
  test_legacy_device_multi_copy_multi_component_entity_values(get_bulk(), *m_field);
}

#endif

} //namespace <anonymous>
