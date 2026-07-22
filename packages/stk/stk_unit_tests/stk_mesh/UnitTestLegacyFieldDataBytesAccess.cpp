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
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
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

//==============================================================================
class LegacyFieldDataBytesAccess : public stk::unit_test_util::MeshFixture
{
public:
  LegacyFieldDataBytesAccess()
    : m_field(nullptr)
  {}

  void build_two_element_mesh() {
    stk::mesh::fixtures::HexFixture::fill_mesh(2, 1, 1, get_bulk());
    m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();  // Trigger creation of default-initialized device data object
  }

  void build_two_bucket_mesh() {
    stk::mesh::fixtures::HexFixture::fill_mesh(stk::mesh::get_default_maximum_bucket_capacity(), 2, 1, get_bulk());
    m_field->data<stk::mesh::ReadOnly, stk::ngp::DeviceSpace>();  // Trigger creation of default-initialized device data object
  }

protected:
  stk::mesh::Field<int>* m_field;
};

//==============================================================================
TEST_F(LegacyFieldDataBytesAccess, host_multiCopyMultiComponent_entityBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh();

  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);
  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      const int elemId = get_bulk().identifier(elem);
      std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                   2*elemId, 20*elemId, 200*elemId};
      int* data = stk::mesh::field_data(*m_field, elem);

      std::memcpy(data, setValue.data(), 6*sizeof(int));
    }
  }

  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      const int elemId = get_bulk().identifier(elem);
      std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                                        2*elemId, 20*elemId, 200*elemId};

      const int* data = stk::mesh::field_data(*m_field, elem);
      for (int i = 0; i < 6; ++i) {
        EXPECT_EQ(data[i], expectedValue[i]);
      }

      EXPECT_EQ(std::memcmp(data, expectedValue.data(), 6*sizeof(int)), 0);
    }
  }
}

TEST_F(LegacyFieldDataBytesAccess, host_multiCopyMultiComponent_bucketBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh();

  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);

  for (stk::mesh::Bucket* bucket : buckets) {
    const int bucketCapacity = bucket->capacity();
    const int scalarsPerEntity = stk::mesh::field_scalars_per_entity(*m_field, *bucket);
    const int bytesPerEntity = stk::mesh::field_bytes_per_entity(*m_field, *bucket);
    int* data = stk::mesh::field_data(*m_field, *bucket);

    for (int bucketOrd = 0; bucketOrd < bucketCapacity; ++bucketOrd) {
      std::array<int, 6> setValue {1*bucketOrd, 10*bucketOrd, 100*bucketOrd,
                                   2*bucketOrd, 20*bucketOrd, 200*bucketOrd};

      std::memcpy(data + bucketOrd*scalarsPerEntity, setValue.data(), bytesPerEntity);
    }
  }

  for (stk::mesh::Bucket* bucket : buckets) {
    const int bucketCapacity = bucket->capacity();
    const int scalarsPerEntity = stk::mesh::field_scalars_per_entity(*m_field, *bucket);
    const int bytesPerEntity = stk::mesh::field_bytes_per_entity(*m_field, *bucket);
    const int* data = stk::mesh::field_data(*m_field, *bucket);

    for (int bucketOrd = 0; bucketOrd < bucketCapacity; ++bucketOrd) {
      std::array<int, 6> expectedValue {1*bucketOrd, 10*bucketOrd, 100*bucketOrd,
                                        2*bucketOrd, 20*bucketOrd, 200*bucketOrd};

      for (int i = 0; i < scalarsPerEntity; ++i) {
        EXPECT_EQ(data[bucketOrd*scalarsPerEntity + i], expectedValue[i]);
      }

      EXPECT_EQ(std::memcmp(data + bucketOrd*scalarsPerEntity, expectedValue.data(), bytesPerEntity), 0);
    }
  }
}


//------------------------------------------------------------------------------
void test_legacy_device_entity_bytes(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::mesh::NgpField<int>& ngpField = stk::mesh::get_updated_ngp_field<int>(field);
  const int bytesPerScalar = sizeof(int);
  const int bytesPerEntity = field.max_size() * bytesPerScalar;

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
      const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
      const int setValue[6] = {1*elemId, 10*elemId, 100*elemId,
                               2*elemId, 20*elemId, 200*elemId};
      const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue);
#ifdef STK_USE_DEVICE_MESH
      const auto& deviceBucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, elem.bucket_id);
      int scalarStride = deviceBucket.capacity() * bytesPerScalar;
#else
      int scalarStride = bytesPerScalar;
#endif

      std::byte* entityBytes = reinterpret_cast<std::byte*>(&ngpField(elem, 0));
      for (int byte = 0; byte < bytesPerEntity; ++byte) {
        const int scalar = byte / bytesPerScalar;
        const int byteInScalar = byte % bytesPerScalar;
        entityBytes[scalar*scalarStride + byteInScalar] = setValueBytes[byte];
      }
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
      const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
      const int expectedValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                    2*elemId, 20*elemId, 200*elemId};
      const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue);
#ifdef STK_USE_DEVICE_MESH
      const auto& deviceBucket = ngpMesh.get_bucket(stk::topology::ELEM_RANK, elem.bucket_id);
      int scalarStride = deviceBucket.capacity() * bytesPerScalar;
#else
      int scalarStride = bytesPerScalar;
#endif

      for (int component = 0; component < 6; ++component) {
        NGP_EXPECT_EQ(ngpField(elem, component), expectedValue[component]);
      }

      std::byte* entityBytes = reinterpret_cast<std::byte*>(&ngpField(elem, 0));
      for (int byte = 0; byte < bytesPerEntity; ++byte) {
        const int scalar = byte / bytesPerScalar;
        const int byteInScalar = byte % bytesPerScalar;
        NGP_EXPECT_EQ(entityBytes[scalar*scalarStride + byteInScalar], expectedValueBytes[byte]);
      }
    }
  );
}

TEST_F(LegacyFieldDataBytesAccess, device_multiCopyMultiComponent_entityBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh();

  test_legacy_device_entity_bytes(get_bulk(), *m_field);
}


//------------------------------------------------------------------------------
TEST_F(LegacyFieldDataBytesAccess, device_multiCopyMultiComponent_bucketBytes)
{
  // We don't have a legacy Bucket-based device-side byte access method
  // that is equivalent to the new API.
}

//==============================================================================

#endif

}
