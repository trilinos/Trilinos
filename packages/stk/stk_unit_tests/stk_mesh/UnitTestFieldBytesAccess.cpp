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

TEST(FieldLayout, ostreamOperator)
{
  {
    std::ostringstream os;
    os<<stk::mesh::Layout::Left;
    EXPECT_EQ(std::string("Layout::Left"), os.str());
  }
  {
    std::ostringstream os;
    os<<stk::mesh::Layout::Right;
    EXPECT_EQ(std::string("Layout::Right"), os.str());
  }
  {
    std::ostringstream os;
    os<<static_cast<stk::mesh::Layout>(42);
    EXPECT_EQ(std::string("Unknown Layout"), os.str());
  }
}

//==============================================================================
class FieldBytesAccess : public stk::unit_test_util::MeshFixture
{
public:
  FieldBytesAccess()
    : m_field(nullptr),
      m_leftField(nullptr)
  {}

  stk::mesh::Entity create_node(stk::mesh::EntityId nodeId) {
    get_bulk().modification_begin();
    stk::mesh::Entity node = get_bulk().declare_node(nodeId);
    get_bulk().modification_end();

    return node;
  }

  template <typename FieldType>
  void build_two_element_mesh(FieldType* field) {
    stk::mesh::fixtures::HexFixture::fill_mesh(2, 1, 1, get_bulk());
    field->template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();  // Trigger creation of default-initialized device data object
  }

  template <typename FieldType>
  void build_two_bucket_mesh(FieldType* field) {
    stk::mesh::fixtures::HexFixture::fill_mesh(stk::mesh::get_default_maximum_bucket_capacity(), 2, 1, get_bulk());
    field->template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();  // Trigger creation of default-initialized device data object
  }

protected:
  stk::mesh::Field<int>* m_field;
  stk::mesh::Field<int, stk::mesh::Layout::Left>* m_leftField;
  stk::mesh::Field<int, stk::mesh::Layout::Right>* m_rightField;
};


//==============================================================================
TEST_F(FieldBytesAccess, host_multiCopyMultiComponent_entityBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  auto fieldBytes = m_field->bytes();
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);

  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      const int elemId = get_bulk().identifier(elem);
      std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                   2*elemId, 20*elemId, 200*elemId};
      const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

      auto entityBytes = fieldBytes.entity_bytes(elem);
      for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
        entityBytes(byte) = setValueBytes[byte];
      }
    }
  }

  auto fieldData = m_field->data<stk::mesh::ReadOnly>();
  auto constFieldBytes = m_field->const_bytes();

  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      const int elemId = get_bulk().identifier(elem);
      std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                                        2*elemId, 20*elemId, 200*elemId};
      const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

      auto entityValues = fieldData.entity_values(elem);
      for (stk::mesh::CopyIdx copy : entityValues.copies()) {
        for (stk::mesh::ComponentIdx component : entityValues.components()) {
          EXPECT_EQ(entityValues(copy, component), expectedValue[3*static_cast<int>(copy) + component]);
        }
      }

      auto constEntityBytes = constFieldBytes.entity_bytes(elem);
      for (stk::mesh::ByteIdx byte : constEntityBytes.bytes()) {
        EXPECT_EQ(constEntityBytes(byte), expectedValueBytes[byte]);
      }
    }
  }
}

TEST_F(FieldBytesAccess, host_multiCopyMultiComponent_entityBytes_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  auto fieldBytes = m_field->bytes();
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);

  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      const int elemId = get_bulk().identifier(elem);
      std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                   2*elemId, 20*elemId, 200*elemId};
      const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

      auto entityBytes = fieldBytes.entity_bytes(elem);
      for (stk::mesh::ByteIdx byte(0); byte < entityBytes.num_bytes(); ++byte) {
        entityBytes(byte) = setValueBytes[byte];
      }
    }
  }

  auto fieldData = m_field->data<stk::mesh::ReadOnly>();
  auto constFieldBytes = m_field->const_bytes();

  for (stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      const int elemId = get_bulk().identifier(elem);
      std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                                        2*elemId, 20*elemId, 200*elemId};
      const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

      auto entityValues = fieldData.entity_values(elem);
      for (stk::mesh::CopyIdx copy(0); copy < entityValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < entityValues.num_components(); ++component) {
          EXPECT_EQ(entityValues(copy, component), expectedValue[3*static_cast<int>(copy) + component]);
        }
      }

      auto constEntityBytes = constFieldBytes.entity_bytes(elem);
      for (stk::mesh::ByteIdx byte : constEntityBytes.bytes()) {
        EXPECT_EQ(constEntityBytes(byte), expectedValueBytes[byte]);
      }
    }
  }
}

TEST_F(FieldBytesAccess, host_multiCopyMultiComponent_bucketBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  auto& fieldBytes = m_field->bytes();
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);

  for (stk::mesh::Bucket* bucket : buckets) {
    auto bucketBytes = fieldBytes.bucket_bytes(*bucket);
    for (stk::mesh::EntityIdx elem : bucket->entities()) {
      const int elemId = get_bulk().identifier((*bucket)[elem]);
      std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                   2*elemId, 20*elemId, 200*elemId};
      const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

      for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
        bucketBytes(elem, byte) = setValueBytes[byte];
      }
    }
  }

  auto constFieldData = m_field->data<stk::mesh::ReadOnly>();
  auto constFieldBytes = m_field->const_bytes();

  for (stk::mesh::Bucket* bucket : buckets) {
    auto constBucketValues = constFieldData.bucket_values(*bucket);
    auto constBucketBytes = constFieldBytes.bucket_bytes(*bucket);

    for (stk::mesh::EntityIdx elem : bucket->entities()) {
      const int elemId = get_bulk().identifier((*bucket)[elem]);
      std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                                        2*elemId, 20*elemId, 200*elemId};
      const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

      for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
        for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
          EXPECT_EQ(constBucketValues(elem, copy, component), expectedValue[3*static_cast<int>(copy) + component]);
        }
      }

      for (stk::mesh::ByteIdx byte : constBucketBytes.bytes()) {
        EXPECT_EQ(constBucketBytes(elem, byte), expectedValueBytes[byte]);
      }
    }
  }
}

TEST_F(FieldBytesAccess, host_multiCopyMultiComponent_bucketBytes_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  auto& fieldBytes = m_field->bytes();
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);

  for (stk::mesh::Bucket* bucket : buckets) {
    auto bucketBytes = fieldBytes.bucket_bytes(*bucket);
    for (stk::mesh::EntityIdx elem(0); elem < bucket->num_entities(); ++elem) {
      const int elemId = get_bulk().identifier((*bucket)[elem]);
      std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                   2*elemId, 20*elemId, 200*elemId};
      const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

      for (stk::mesh::ByteIdx byte(0); byte < bucketBytes.num_bytes(); ++byte) {
        bucketBytes(elem, byte) = setValueBytes[byte];
      }
    }
  }

  auto constFieldData = m_field->data<stk::mesh::ReadOnly>();
  auto constFieldBytes = m_field->const_bytes();

  for (stk::mesh::Bucket* bucket : buckets) {
    auto constBucketValues = constFieldData.bucket_values(*bucket);
    auto constBucketBytes = constFieldBytes.bucket_bytes(*bucket);

    for (stk::mesh::EntityIdx elem(0); elem < bucket->num_entities(); ++elem) {
      const int elemId = get_bulk().identifier((*bucket)[elem]);
      std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                                        2*elemId, 20*elemId, 200*elemId};
      const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

      for (stk::mesh::CopyIdx copy(0); copy < constBucketValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constBucketValues.num_components(); ++component) {
          EXPECT_EQ(constBucketValues(elem, copy, component), expectedValue[3*static_cast<int>(copy) + component]);
        }
      }

      for (stk::mesh::ByteIdx byte(0); byte < constBucketBytes.num_bytes(); ++byte) {
        EXPECT_EQ(constBucketBytes(elem, byte), expectedValueBytes[byte]);
      }
    }
  }
}

TEST_F(FieldBytesAccess, host_multiCopyMultiComponent_bucketBytes_fieldBytesIndexing)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  auto& fieldBytes = m_field->bytes();
  const stk::mesh::BucketVector& buckets = get_bulk().buckets(stk::topology::ELEM_RANK);

  for (stk::mesh::Bucket* bucket : buckets) {
    auto bucketBytes = fieldBytes.bucket_bytes(*bucket);
    for (stk::mesh::EntityIdx elem : bucketBytes.entities()) {
      const int elemId = get_bulk().identifier((*bucket)[elem]);
      std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                   2*elemId, 20*elemId, 200*elemId};
      const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

      for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
        bucketBytes(elem, byte) = setValueBytes[byte];
      }
    }
  }

  auto constFieldData = m_field->data<stk::mesh::ReadOnly>();
  auto constFieldBytes = m_field->const_bytes();

  for (stk::mesh::Bucket* bucket : buckets) {
    auto constBucketValues = constFieldData.bucket_values(*bucket);
    auto constBucketBytes = constFieldBytes.bucket_bytes(*bucket);

    for (stk::mesh::EntityIdx elem : constBucketBytes.entities()) {
      const int elemId = get_bulk().identifier((*bucket)[elem]);
      std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                                        2*elemId, 20*elemId, 200*elemId};
      const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

      for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
        for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
          EXPECT_EQ(constBucketValues(elem, copy, component), expectedValue[3*static_cast<int>(copy) + component]);
        }
      }

      for (stk::mesh::ByteIdx byte : constBucketBytes.bytes()) {
        EXPECT_EQ(constBucketBytes(elem, byte), expectedValueBytes[byte]);
      }
    }
  }
}


//==============================================================================
void test_device_entity_bytes(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  auto& fieldBytes = field.bytes<stk::ngp::MemSpace>();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
      const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
      const int setValue[6] = {1*elemId, 10*elemId, 100*elemId,
                               2*elemId, 20*elemId, 200*elemId};
      const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue);

      auto entityBytes = fieldBytes.entity_bytes(elem);
      for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
        entityBytes(byte) = setValueBytes[byte];
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  auto& constFieldBytes = field.const_bytes<stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
      const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
      const int expectedValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                    2*elemId, 20*elemId, 200*elemId};
      const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue);

      auto constFieldValues = constFieldData.entity_values(elem);
      const int numComponents = constFieldValues.num_components();
      for (stk::mesh::CopyIdx copy : constFieldValues.copies()) {
        for (stk::mesh::ComponentIdx component : constFieldValues.components()) {
          NGP_EXPECT_EQ(constFieldValues(copy, component),
                        expectedValue[static_cast<int>(copy)*numComponents + component]);
        }
      }

      auto constEntityBytes = constFieldBytes.entity_bytes(elem);
      for (stk::mesh::ByteIdx byte : constEntityBytes.bytes()) {
        NGP_EXPECT_EQ(constEntityBytes(byte), expectedValueBytes[byte]);
      }
    }
  );
}

TEST_F(FieldBytesAccess, device_multiCopyMultiComponent_entityBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  test_device_entity_bytes(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_entity_bytes_traditional_for_loop(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  auto& fieldBytes = field.bytes<stk::ngp::MemSpace>();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
      const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
      const int setValue[6] = {1*elemId, 10*elemId, 100*elemId,
                               2*elemId, 20*elemId, 200*elemId};
      const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue);

      auto entityBytes = fieldBytes.entity_bytes(elem);
      for (stk::mesh::ByteIdx byte(0); byte <entityBytes.num_bytes(); ++byte) {
        entityBytes(byte) = setValueBytes[byte];
      }
    }
  );

  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  auto& constFieldBytes = field.const_bytes<stk::ngp::MemSpace>();
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
      const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
      const int expectedValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                    2*elemId, 20*elemId, 200*elemId};
      const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue);

      auto constFieldValues = constFieldData.entity_values(elem);
      const int numComponents = constFieldValues.num_components();
      for (stk::mesh::CopyIdx copy(0); copy < constFieldValues.num_copies(); ++copy) {
        for (stk::mesh::ComponentIdx component(0); component < constFieldValues.num_components(); ++component) {
          NGP_EXPECT_EQ(constFieldValues(copy, component),
                        expectedValue[static_cast<int>(copy)*numComponents + component]);
        }
      }

      auto constEntityBytes = constFieldBytes.entity_bytes(elem);
      for (stk::mesh::ByteIdx byte(0); byte < constEntityBytes.num_bytes(); ++byte) {
        NGP_EXPECT_EQ(constEntityBytes(byte), expectedValueBytes[byte]);
      }
    }
  );
}

TEST_F(FieldBytesAccess, device_multiCopyMultiComponent_entityBytes_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  test_device_entity_bytes_traditional_for_loop(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_bucket_bytes(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  auto& fieldBytes = field.bytes<stk::ngp::MemSpace>();

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::ELEM_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketBytes = fieldBytes.bucket_bytes(bucketId);

      const stk::mesh::EntityIdx numElems = bucketBytes.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numElems),
        [&](stk::mesh::EntityIdx elem) {
          const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK,
                                                                   stk::mesh::FastMeshIndex{static_cast<unsigned>(bucketId),
                                                                                            static_cast<unsigned>(elem)}));
          const int setValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                   2*elemId, 20*elemId, 200*elemId};
          const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue);

          for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
            bucketBytes(elem, byte) = setValueBytes[byte];
          }
        }
      );
    }
  );

  auto& constFieldBytes = field.const_bytes<stk::ngp::MemSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = constFieldData.bucket_values(bucketId);
      auto bucketBytes = constFieldBytes.bucket_bytes(bucketId);

      const stk::mesh::EntityIdx numElems = bucketBytes.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numElems),
        [&](stk::mesh::EntityIdx elem) {
          const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK,
                                                                   stk::mesh::FastMeshIndex{static_cast<unsigned>(bucketId),
                                                                                            static_cast<unsigned>(elem)}));
          const int expectedValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                        2*elemId, 20*elemId, 200*elemId};
          const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue);

          const int numComponents = bucketValues.num_components();
          for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
            for (stk::mesh::ComponentIdx component : bucketValues.components()) {
              NGP_EXPECT_EQ(bucketValues(elem, copy, component),
                            expectedValue[static_cast<int>(copy)*numComponents + component]);
            }
          }

          for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
            NGP_EXPECT_EQ(bucketBytes(elem, byte), expectedValueBytes[byte]);
          }
        }
      );
    }
  );
}

TEST_F(FieldBytesAccess, device_multiCopyMultiComponent_bucketBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  test_device_bucket_bytes(get_bulk(), *m_field);
}

//------------------------------------------------------------------------------
void test_device_bucket_bytes_traditional_for_loop(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  auto constFieldData = field.data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
  auto& fieldBytes = field.bytes<stk::ngp::MemSpace>();

  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::ELEM_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketBytes = fieldBytes.bucket_bytes(bucketId);

      const stk::mesh::EntityIdx numElems = bucketBytes.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numElems),
        [&](stk::mesh::EntityIdx elem) {
          const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK,
                                                                   stk::mesh::FastMeshIndex{static_cast<unsigned>(bucketId),
                                                                                            static_cast<unsigned>(elem)}));
          const int setValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                   2*elemId, 20*elemId, 200*elemId};
          const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue);

          for (stk::mesh::ByteIdx byte(0); byte < bucketBytes.num_bytes(); ++byte) {
            bucketBytes(elem, byte) = setValueBytes[byte];
          }
        }
      );
    }
  );

  auto& constFieldBytes = field.const_bytes<stk::ngp::MemSpace>();
  Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
    KOKKOS_LAMBDA(const TeamHandleType& team) {
      const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
      auto bucketValues = constFieldData.bucket_values(bucketId);
      auto bucketBytes = constFieldBytes.bucket_bytes(bucketId);

      const stk::mesh::EntityIdx numElems = bucketBytes.num_entities();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numElems),
        [&](stk::mesh::EntityIdx elem) {
          const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK,
                                                                   stk::mesh::FastMeshIndex{static_cast<unsigned>(bucketId),
                                                                                            static_cast<unsigned>(elem)}));
          const int expectedValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                        2*elemId, 20*elemId, 200*elemId};
          const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue);

          const int numComponents = bucketValues.num_components();
          for (stk::mesh::CopyIdx copy(0); copy < bucketValues.num_copies(); ++copy) {
            for (stk::mesh::ComponentIdx component(0); component < bucketValues.num_components(); ++component) {
              NGP_EXPECT_EQ(bucketValues(elem, copy, component),
                            expectedValue[static_cast<int>(copy)*numComponents + component]);
            }
          }

          for (stk::mesh::ByteIdx byte(0); byte < bucketBytes.num_bytes(); ++byte) {
            NGP_EXPECT_EQ(bucketBytes(elem, byte), expectedValueBytes[byte]);
          }
        }
      );
    }
  );
}

TEST_F(FieldBytesAccess, device_multiCopyMultiComponent_bucketBytes_traditionalForLoop)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  test_device_bucket_bytes_traditional_for_loop(get_bulk(), *m_field);
}


//==============================================================================
template <typename FieldType>
void test_host_to_device_entity_bytes(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  {
    field.template synchronize<stk::mesh::ReadWrite>();  // Mark it as modified so that we will sync on the other side
    auto fieldBytes = field.bytes();
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity elem : *bucket) {
        const int elemId = bulk.identifier(elem);
        std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                     2*elemId, 20*elemId, 200*elemId};
        const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

        auto entityBytes = fieldBytes.entity_bytes(elem);
        for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
          entityBytes(byte) = setValueBytes[byte];
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
    auto& constFieldBytes = field.template const_bytes<stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
        const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
        const int expectedValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                      2*elemId, 20*elemId, 200*elemId};
        const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue);

        auto constFieldValues = constFieldData.entity_values(elem);
        const int numComponents = constFieldValues.num_components();
        for (stk::mesh::CopyIdx copy : constFieldValues.copies()) {
          for (stk::mesh::ComponentIdx component : constFieldValues.components()) {
            NGP_EXPECT_EQ(constFieldValues(copy, component),
                          expectedValue[static_cast<int>(copy)*numComponents + component]);
          }
        }

        auto constEntityBytes = constFieldBytes.entity_bytes(elem);
        for (stk::mesh::ByteIdx byte : constEntityBytes.bytes()) {
          NGP_EXPECT_EQ(constEntityBytes(byte), expectedValueBytes[byte]);
        }
      }
    );
  }
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_entityBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  test_host_to_device_entity_bytes(get_bulk(), *m_field);
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_entityBytes_automaticLayoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "leftField1");
  stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_leftField);

  test_host_to_device_entity_bytes(get_bulk(), *m_leftField);
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_entityBytes_automaticLayoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "rightField1");
  stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_rightField);

  test_host_to_device_entity_bytes(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_entity_bytes_forced_layout(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  {
    field.template synchronize<stk::mesh::ReadWrite>();  // Mark it as modified so that we will sync on the other side
    auto fieldBytes = field.bytes();
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity elem : *bucket) {
        const int elemId = bulk.identifier(elem);
        std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                     2*elemId, 20*elemId, 200*elemId};
        const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

        if (fieldBytes.data_layout() == stk::mesh::Layout::Left) {
          auto entityBytes = fieldBytes.entity_bytes_left(elem);
          for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
            entityBytes(byte) = setValueBytes[byte];
          }
        }
        else {
          auto entityBytes = fieldBytes.entity_bytes_right(elem);
          for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
            entityBytes(byte) = setValueBytes[byte];
          }
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
    auto& constFieldBytes = field.template const_bytes<stk::ngp::MemSpace>();
    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
        const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
        const int expectedValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                      2*elemId, 20*elemId, 200*elemId};
        const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue);

        auto constFieldValues = constFieldData.entity_values(elem);
        const int numComponents = constFieldValues.num_components();
        for (stk::mesh::CopyIdx copy : constFieldValues.copies()) {
          for (stk::mesh::ComponentIdx component : constFieldValues.components()) {
            NGP_EXPECT_EQ(constFieldValues(copy, component),
                          expectedValue[static_cast<int>(copy)*numComponents + component]);
          }
        }

        auto constEntityBytes = constFieldBytes.entity_bytes(elem);
        for (stk::mesh::ByteIdx byte : constEntityBytes.bytes()) {
          NGP_EXPECT_EQ(constEntityBytes(byte), expectedValueBytes[byte]);
        }
      }
    );
  }
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_entityBytes_forcedLayoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "leftField1");
  stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_leftField);

  test_host_to_device_entity_bytes_forced_layout(get_bulk(), *m_leftField);
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_entityBytes_forcedLayoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "rightField1");
  stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_rightField);

  test_host_to_device_entity_bytes_forced_layout(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_bucket_bytes(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::ELEM_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    field.template synchronize<stk::mesh::ReadWrite>();  // Mark it as modified so that we will sync on the other side
    auto& fieldBytes = field.bytes();
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketBytes = fieldBytes.bucket_bytes(*bucket);
      for (stk::mesh::EntityIdx elem : bucket->entities()) {
        const int elemId = bulk.identifier((*bucket)[elem]);
        std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                     2*elemId, 20*elemId, 200*elemId};
        const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

        for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
          bucketBytes(elem, byte) = setValueBytes[byte];
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
    auto& constFieldBytes = field.template const_bytes<stk::ngp::MemSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto bucketValues = constFieldData.bucket_values(bucketId);
        auto bucketBytes = constFieldBytes.bucket_bytes(bucketId);

        const stk::mesh::EntityIdx numElems = bucketBytes.num_entities();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numElems),
          [&](stk::mesh::EntityIdx elem) {
            const int elemId = ngpMesh.identifier(
                  ngpMesh.get_entity(stk::topology::ELEM_RANK,
                  stk::mesh::FastMeshIndex{static_cast<unsigned>(bucketId),
                                           static_cast<unsigned>(elem)}));
            const int expectedValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                          2*elemId, 20*elemId, 200*elemId};
            const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue);

            const int numComponents = bucketValues.num_components();
            for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
              for (stk::mesh::ComponentIdx component : bucketValues.components()) {
                NGP_EXPECT_EQ(bucketValues(elem, copy, component),
                              expectedValue[static_cast<int>(copy)*numComponents + component]);
              }
            }

            for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
              NGP_EXPECT_EQ(bucketBytes(elem, byte), expectedValueBytes[byte]);
            }
          }
        );
      }
    );
  }
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_bucketBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_field);

  test_host_to_device_bucket_bytes(get_bulk(), *m_field);
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_bucketBytes_automaticLayoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "leftField1");
  stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_leftField);

  test_host_to_device_bucket_bytes(get_bulk(), *m_leftField);
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_bucketBytes_automaticLayoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "rightField1");
  stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_rightField);

  test_host_to_device_bucket_bytes(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_host_to_device_bucket_bytes_forced_layout(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::ELEM_RANK, field);
  unsigned numBuckets = bucketIds.size();
  using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

  {
    field.template synchronize<stk::mesh::ReadWrite>();  // Mark it as modified so that we will sync on the other side
    auto& fieldBytes = field.bytes();
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      if (fieldBytes.data_layout() == stk::mesh::Layout::Left) {
        auto bucketBytes = fieldBytes.bucket_bytes_left(*bucket);
        for (stk::mesh::EntityIdx elem : bucket->entities()) {
          const int elemId = bulk.identifier((*bucket)[elem]);
          std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                       2*elemId, 20*elemId, 200*elemId};
          const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

          for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
            bucketBytes(elem, byte) = setValueBytes[byte];
          }
        }
      }
      else {
        auto bucketBytes = fieldBytes.bucket_bytes_right(*bucket);
        for (stk::mesh::EntityIdx elem : bucket->entities()) {
          const int elemId = bulk.identifier((*bucket)[elem]);
          std::array<int, 6> setValue {1*elemId, 10*elemId, 100*elemId,
                                       2*elemId, 20*elemId, 200*elemId};
          const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue.data());

          for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
            bucketBytes(elem, byte) = setValueBytes[byte];
          }
        }
      }
    }
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly, stk::ngp::MemSpace>();
    auto& constFieldBytes = field.template const_bytes<stk::ngp::MemSpace>();
    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto bucketValues = constFieldData.bucket_values(bucketId);
        auto bucketBytes = constFieldBytes.bucket_bytes(bucketId);

        const stk::mesh::EntityIdx numElems = bucketBytes.num_entities();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numElems),
          [&](stk::mesh::EntityIdx elem) {
            const int elemId = ngpMesh.identifier(
                  ngpMesh.get_entity(stk::topology::ELEM_RANK,
                  stk::mesh::FastMeshIndex{static_cast<unsigned>(bucketId),
                                           static_cast<unsigned>(elem)}));
            const int expectedValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                          2*elemId, 20*elemId, 200*elemId};
            const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue);

            const int numComponents = bucketValues.num_components();
            for (stk::mesh::CopyIdx copy : bucketValues.copies()) {
              for (stk::mesh::ComponentIdx component : bucketValues.components()) {
                NGP_EXPECT_EQ(bucketValues(elem, copy, component),
                              expectedValue[static_cast<int>(copy)*numComponents + component]);
              }
            }

            for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
              NGP_EXPECT_EQ(bucketBytes(elem, byte), expectedValueBytes[byte]);
            }
          }
        );
      }
    );
  }
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_bucketBytes_forcedLayoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "leftField1");
  stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_leftField);

  test_host_to_device_bucket_bytes_forced_layout(get_bulk(), *m_leftField);
}

TEST_F(FieldBytesAccess, mixedHostToDevice_multiCopyMultiComponent_bucketBytes_forcedLayoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "rightField1");
  stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_rightField);

  test_host_to_device_bucket_bytes_forced_layout(get_bulk(), *m_rightField);
}



//==============================================================================
template <typename FieldType>
void test_device_to_host_entity_bytes(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  {
    field.template synchronize<stk::mesh::OverwriteAll, stk::ngp::MemSpace>();  // Mark it as modified so that we will sync on the other side
    auto fieldBytes = field.template bytes<stk::ngp::MemSpace>();

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
        const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
        const int setValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                 2*elemId, 20*elemId, 200*elemId};
        const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue);

        auto entityBytes = fieldBytes.entity_bytes(elem);
        for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
          entityBytes(byte) = setValueBytes[byte];
        }
      }
    );
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly>();
    auto constFieldBytes = field.const_bytes();
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity elem : *bucket) {
        const int elemId = bulk.identifier(elem);
        std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                                          2*elemId, 20*elemId, 200*elemId};
        const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

        auto entityValues = constFieldData.entity_values(elem);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            EXPECT_EQ(entityValues(copy, component), expectedValue[3*static_cast<int>(copy) + component]);
          }
        }

        auto constEntityBytes = constFieldBytes.entity_bytes(elem);
        for (stk::mesh::ByteIdx byte : constEntityBytes.bytes()) {
          EXPECT_EQ(constEntityBytes(byte), expectedValueBytes[byte]);
        }
      }
    }
  }
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_entityBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_field);

  test_device_to_host_entity_bytes(get_bulk(), *m_field);
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_entityBytes_automaticLayoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "leftField1");
  stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_leftField);

  test_device_to_host_entity_bytes(get_bulk(), *m_leftField);
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_entityBytes_automaticLayoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "rightField1");
  stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_rightField);

  test_device_to_host_entity_bytes(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_entity_bytes_forced_layout(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  {
    field.template synchronize<stk::mesh::OverwriteAll, stk::ngp::MemSpace>();  // Mark it as modified so that we will sync on the other side
    auto fieldBytes = field.template bytes<stk::ngp::MemSpace>();

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, field,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem) {
        const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK, elem));
        const int setValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                 2*elemId, 20*elemId, 200*elemId};
        const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue);

        auto entityBytes = fieldBytes.entity_bytes(elem);
        for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
          entityBytes(byte) = setValueBytes[byte];
        }
      }
    );
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly>();
    auto constFieldBytes = field.const_bytes();
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      for (stk::mesh::Entity elem : *bucket) {
        const int elemId = bulk.identifier(elem);
        std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                                          2*elemId, 20*elemId, 200*elemId};
        const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

        auto entityValues = constFieldData.entity_values(elem);
        for (stk::mesh::CopyIdx copy : entityValues.copies()) {
          for (stk::mesh::ComponentIdx component : entityValues.components()) {
            EXPECT_EQ(entityValues(copy, component), expectedValue[3*static_cast<int>(copy) + component]);
          }
        }

        if (constFieldBytes.data_layout() == stk::mesh::Layout::Left) {
          auto constEntityBytes = constFieldBytes.entity_bytes_left(elem);
          for (stk::mesh::ByteIdx byte : constEntityBytes.bytes()) {
            EXPECT_EQ(constEntityBytes(byte), expectedValueBytes[byte]);
          }
        }
        else {
          auto constEntityBytes = constFieldBytes.entity_bytes_right(elem);
          for (stk::mesh::ByteIdx byte : constEntityBytes.bytes()) {
            EXPECT_EQ(constEntityBytes(byte), expectedValueBytes[byte]);
          }
        }
      }
    }
  }
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_entityBytes_forcedLayoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "leftField1");
  stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_leftField);

  test_device_to_host_entity_bytes_forced_layout(get_bulk(), *m_leftField);
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_entityBytes_forcedLayoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 4, 4);
  m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "rightField1");
  stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_element_mesh(m_rightField);

  test_device_to_host_entity_bytes_forced_layout(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_bucket_bytes(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  {
    field.template synchronize<stk::mesh::ReadWrite, stk::ngp::MemSpace>();  // Mark it as modified so that we will sync on the other side
    auto fieldBytes = field.template bytes<stk::ngp::MemSpace>();

    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::ELEM_RANK, field);
    unsigned numBuckets = bucketIds.size();
    using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto bucketBytes = fieldBytes.bucket_bytes(bucketId);

        const stk::mesh::EntityIdx numElems = bucketBytes.num_entities();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numElems),
          [&](stk::mesh::EntityIdx elem) {
            const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK,
                                                                     stk::mesh::FastMeshIndex{static_cast<unsigned>(bucketId),
                                                                                              static_cast<unsigned>(elem)}));
            const int setValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                     2*elemId, 20*elemId, 200*elemId};
            const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue);

            for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
              bucketBytes(elem, byte) = setValueBytes[byte];
            }
          }
        );
      }
    );
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly>();
    auto constFieldBytes = field.const_bytes();
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      auto constBucketValues = constFieldData.bucket_values(*bucket);
      auto constBucketBytes = constFieldBytes.bucket_bytes(*bucket);

      for (stk::mesh::EntityIdx elem : bucket->entities()) {
        const int elemId = bulk.identifier((*bucket)[elem]);
        std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                                          2*elemId, 20*elemId, 200*elemId};
        const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

        for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
          for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
            EXPECT_EQ(constBucketValues(elem, copy, component), expectedValue[3*static_cast<int>(copy) + component]);
          }
        }

        for (stk::mesh::ByteIdx byte : constBucketBytes.bytes()) {
          EXPECT_EQ(constBucketBytes(elem, byte), expectedValueBytes[byte]);
        }
      }
    }
  }
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_bucketBytes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_field = &get_meta().declare_field<int>(stk::topology::ELEM_RANK, "field1");
  stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_field);

  test_device_to_host_bucket_bytes(get_bulk(), *m_field);
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_bucketBytes_automaticLayoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "leftField1");
  stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_leftField);

  test_device_to_host_bucket_bytes(get_bulk(), *m_leftField);
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_bucketBytes_automaticLayoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "rightField1");
  stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_rightField);

  test_device_to_host_bucket_bytes(get_bulk(), *m_rightField);
}

//------------------------------------------------------------------------------
template <typename FieldType>
void test_device_to_host_bucket_bytes_forced_layout(stk::mesh::BulkData& bulk, FieldType& field)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  {
    field.template synchronize<stk::mesh::ReadWrite, stk::ngp::MemSpace>();  // Mark it as modified so that we will sync on the other side
    auto fieldBytes = field.template bytes<stk::ngp::MemSpace>();

    stk::NgpVector<unsigned> bucketIds = ngpMesh.get_bucket_ids(stk::topology::ELEM_RANK, field);
    unsigned numBuckets = bucketIds.size();
    using TeamHandleType = typename stk::ngp::TeamPolicy<stk::ngp::ExecSpace>::member_type;

    Kokkos::parallel_for(stk::ngp::TeamPolicy<stk::ngp::ExecSpace>(numBuckets, Kokkos::AUTO),
      KOKKOS_LAMBDA(const TeamHandleType& team) {
        const int bucketId = bucketIds.get<stk::ngp::ExecSpace>(team.league_rank());
        auto bucketBytes = fieldBytes.bucket_bytes(bucketId);

        const stk::mesh::EntityIdx numElems = bucketBytes.num_entities();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0_entity, numElems),
          [&](stk::mesh::EntityIdx elem) {
            const int elemId = ngpMesh.identifier(ngpMesh.get_entity(stk::topology::ELEM_RANK,
                                                                     stk::mesh::FastMeshIndex{static_cast<unsigned>(bucketId),
                                                                                              static_cast<unsigned>(elem)}));
            const int setValue[6] = {1*elemId, 10*elemId, 100*elemId,
                                     2*elemId, 20*elemId, 200*elemId};
            const std::byte* setValueBytes = reinterpret_cast<const std::byte*>(setValue);

            for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
              bucketBytes(elem, byte) = setValueBytes[byte];
            }
          }
        );
      }
    );
  }

  {
    auto constFieldData = field.template data<stk::mesh::ReadOnly>();
    auto constFieldBytes = field.const_bytes();
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      if (constFieldBytes.data_layout() == stk::mesh::Layout::Left) {
        auto constBucketValues = constFieldData.bucket_values(*bucket);
        auto constBucketBytes = constFieldBytes.bucket_bytes_left(*bucket);

        for (stk::mesh::EntityIdx elem : bucket->entities()) {
          const int elemId = bulk.identifier((*bucket)[elem]);
          std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                2*elemId, 20*elemId, 200*elemId};
          const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

          for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
            for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
              EXPECT_EQ(constBucketValues(elem, copy, component), expectedValue[3*static_cast<int>(copy) + component]);
            }
          }

          for (stk::mesh::ByteIdx byte : constBucketBytes.bytes()) {
            EXPECT_EQ(constBucketBytes(elem, byte), expectedValueBytes[byte]);
          }
        }
      }
      else {
        auto constBucketValues = constFieldData.bucket_values(*bucket);
        auto constBucketBytes = constFieldBytes.bucket_bytes_right(*bucket);

        for (stk::mesh::EntityIdx elem : bucket->entities()) {
          const int elemId = bulk.identifier((*bucket)[elem]);
          std::array<int, 6> expectedValue {1*elemId, 10*elemId, 100*elemId,
                2*elemId, 20*elemId, 200*elemId};
          const std::byte* expectedValueBytes = reinterpret_cast<const std::byte*>(expectedValue.data());

          for (stk::mesh::CopyIdx copy : constBucketValues.copies()) {
            for (stk::mesh::ComponentIdx component : constBucketValues.components()) {
              EXPECT_EQ(constBucketValues(elem, copy, component), expectedValue[3*static_cast<int>(copy) + component]);
            }
          }

          for (stk::mesh::ByteIdx byte : constBucketBytes.bytes()) {
            EXPECT_EQ(constBucketBytes(elem, byte), expectedValueBytes[byte]);
          }
        }
      }
    }
  }
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_bucketBytes_forcedLayoutLeft)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_leftField = &get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "leftField1");
  stk::mesh::put_field_on_mesh(*m_leftField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_leftField);

  test_device_to_host_bucket_bytes_forced_layout(get_bulk(), *m_leftField);
}

TEST_F(FieldBytesAccess, mixedDeviceToHost_multiCopyMultiComponent_bucketBytes_forcedLayoutRight)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) GTEST_SKIP();

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  m_rightField = &get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "rightField1");
  stk::mesh::put_field_on_mesh(*m_rightField, get_meta().universal_part(), 3, 2, nullptr);
  build_two_bucket_mesh(m_rightField);

  test_device_to_host_bucket_bytes_forced_layout(get_bulk(), *m_rightField);
}


//==============================================================================
TEST_F(FieldBytesAccess, host_isFieldDefined)
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

  auto fieldBytes = field.const_bytes();
  auto elem1Bytes = fieldBytes.entity_bytes(elem1);
  auto elem2Bytes = fieldBytes.entity_bytes(elem2);

  EXPECT_EQ(elem1Bytes.is_field_defined(), true);
  EXPECT_EQ(elem2Bytes.is_field_defined(), false);

  const stk::mesh::BucketVector& buckets1 = get_bulk().get_buckets(stk::topology::ELEM_RANK, part1);
  const stk::mesh::BucketVector& buckets2 = get_bulk().get_buckets(stk::topology::ELEM_RANK, part2);
  auto bucket1Bytes = fieldBytes.bucket_bytes(*buckets1.front());
  auto bucket2Bytes = fieldBytes.bucket_bytes(*buckets2.front());

  EXPECT_EQ(bucket1Bytes.is_field_defined(), true);
  EXPECT_EQ(bucket2Bytes.is_field_defined(), false);
}

//------------------------------------------------------------------------------
void device_is_field_defined(stk::mesh::BulkData& bulk, stk::mesh::Field<int>& field,
                             stk::mesh::Part& part1, stk::mesh::Part& part2)
{
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  auto fieldBytes = field.const_bytes();

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, part1,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto elemBytes = fieldBytes.entity_bytes(entity);
      NGP_EXPECT_EQ(elemBytes.is_field_defined(), true);

      auto bucketBytes = fieldBytes.bucket_bytes(entity.bucket_id);
      NGP_EXPECT_EQ(bucketBytes.is_field_defined(), true);
    }
  );

  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, part2,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entity) {
      auto elemBytes = fieldBytes.entity_bytes(entity);
      NGP_EXPECT_EQ(elemBytes.is_field_defined(), false);

      auto bucketBytes = fieldBytes.bucket_bytes(entity.bucket_id);
      NGP_EXPECT_EQ(bucketBytes.is_field_defined(), false);
    }
  );
}

TEST_F(FieldBytesAccess, device_isFieldDefined)
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

// Note: We can only test host-side consistency and bounds-checking because the
// device-side checking issues a Kokkos::abort() which cannot be trapped in a test.

#ifdef STK_FIELD_BOUNDS_CHECK

//==============================================================================
TEST_F(FieldBytesAccess, host_consistencyCheck_entity)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  auto& elemField      = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField");
  auto& elemFieldLeft  = get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "elemFieldLeft");
  auto& elemFieldRight = get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "elemFieldRight");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldLeft, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldRight, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x1", get_bulk());
  const stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);

  EXPECT_ANY_THROW(elemField.bytes().entity_bytes(node1));                   // Wrong rank entity
  EXPECT_ANY_THROW(elemField.const_bytes().entity_bytes(node1));             // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldLeft.bytes().entity_bytes_left(node1));          // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldLeft.const_bytes().entity_bytes_left(node1));    // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldRight.bytes().entity_bytes_right(node1));        // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldRight.const_bytes().entity_bytes_right(node1));  // Wrong rank entity
}

//==============================================================================
TEST_F(FieldBytesAccess, host_consistencyCheck_meshIndex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  auto& elemField      = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField");
  auto& elemFieldLeft  = get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "elemFieldLeft");
  auto& elemFieldRight = get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "elemFieldRight");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldLeft, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldRight, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x1", get_bulk());
  const stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);

  const stk::mesh::MeshIndex node1_mi = get_bulk().mesh_index(node1);
  EXPECT_ANY_THROW(elemField.bytes().entity_bytes(node1_mi));                   // Wrong rank entity
  EXPECT_ANY_THROW(elemField.const_bytes().entity_bytes(node1_mi));             // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldLeft.bytes().entity_bytes_left(node1_mi));          // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldLeft.const_bytes().entity_bytes_left(node1_mi));    // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldRight.bytes().entity_bytes_right(node1_mi));        // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldRight.const_bytes().entity_bytes_right(node1_mi));  // Wrong rank entity

  auto secondMesh = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_initial_bucket_capacity(1)
                                                          .set_maximum_bucket_capacity(1).create();
  stk::io::fill_mesh("generated:1x1x1", *secondMesh);  // Create two-element mesh
  const stk::mesh::Entity elem1 = secondMesh->get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::MeshIndex elem1_mi = secondMesh->mesh_index(elem1);
  EXPECT_ANY_THROW(elemField.bytes().entity_bytes(elem1_mi));                   // Entity from different mesh
  EXPECT_ANY_THROW(elemField.const_bytes().entity_bytes(elem1_mi));             // Entity from different mesh
  EXPECT_ANY_THROW(elemFieldLeft.bytes().entity_bytes_left(elem1_mi));          // Entity from different mesh
  EXPECT_ANY_THROW(elemFieldLeft.const_bytes().entity_bytes_left(elem1_mi));    // Entity from different mesh
  EXPECT_ANY_THROW(elemFieldRight.bytes().entity_bytes_right(elem1_mi));        // Entity from different mesh
  EXPECT_ANY_THROW(elemFieldRight.const_bytes().entity_bytes_right(elem1_mi));  // Entity from different mesh
}

//------------------------------------------------------------------------------
TEST_F(FieldBytesAccess, host_consistencyCheck_fastMeshIndex)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  auto& elemField      = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField");
  auto& elemFieldLeft  = get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "elemFieldLeft");
  auto& elemFieldRight = get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "elemFieldRight");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldLeft, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldRight, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  const stk::mesh::FastMeshIndex elem1_badFmi{1, 0};
  EXPECT_ANY_THROW(elemField.bytes().entity_bytes(elem1_badFmi));                   // Bad Bucket Id
  EXPECT_ANY_THROW(elemField.const_bytes().entity_bytes(elem1_badFmi));             // Bad Bucket Id
  EXPECT_ANY_THROW(elemFieldLeft.bytes().entity_bytes_left(elem1_badFmi));          // Bad Bucket Id
  EXPECT_ANY_THROW(elemFieldLeft.const_bytes().entity_bytes_left(elem1_badFmi));    // Bad Bucket Id
  EXPECT_ANY_THROW(elemFieldRight.bytes().entity_bytes_right(elem1_badFmi));        // Bad Bucket Id
  EXPECT_ANY_THROW(elemFieldRight.const_bytes().entity_bytes_right(elem1_badFmi));  // Bad Bucket Id
}

//------------------------------------------------------------------------------
TEST_F(FieldBytesAccess, host_consistencyCheck_bucket)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  auto& elemField      = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField");
  auto& elemFieldLeft  = get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "elemFieldLeft");
  auto& elemFieldRight = get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "elemFieldRight");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldLeft, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldRight, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x1", get_bulk());
  const stk::mesh::Entity node1 = get_bulk().get_entity(stk::topology::NODE_RANK, 1);

  const stk::mesh::Bucket& bucket_node1 = get_bulk().bucket(node1);
  EXPECT_ANY_THROW(elemField.bytes().bucket_bytes(bucket_node1));                   // Wrong rank entity
  EXPECT_ANY_THROW(elemField.const_bytes().bucket_bytes(bucket_node1));             // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldLeft.bytes().bucket_bytes_left(bucket_node1));          // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldLeft.const_bytes().bucket_bytes_left(bucket_node1));    // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldRight.bytes().bucket_bytes_right(bucket_node1));        // Wrong rank entity
  EXPECT_ANY_THROW(elemFieldRight.const_bytes().bucket_bytes_right(bucket_node1));  // Wrong rank entity

  auto secondMesh = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_initial_bucket_capacity(1)
                                                          .set_maximum_bucket_capacity(1).create();
  stk::io::fill_mesh("generated:1x1x1", *secondMesh);  // Create two-element mesh
  const stk::mesh::Entity elem1 = secondMesh->get_entity(stk::topology::ELEM_RANK, 1);
  const stk::mesh::Bucket& bucket_elem1 = secondMesh->bucket(elem1);
  EXPECT_ANY_THROW(elemField.bytes().bucket_bytes(bucket_elem1));                   // Bucket from different mesh
  EXPECT_ANY_THROW(elemField.const_bytes().bucket_bytes(bucket_elem1));             // Bucket from different mesh
  EXPECT_ANY_THROW(elemFieldLeft.bytes().bucket_bytes_left(bucket_elem1));          // Bucket from different mesh
  EXPECT_ANY_THROW(elemFieldLeft.const_bytes().bucket_bytes_left(bucket_elem1));    // Bucket from different mesh
  EXPECT_ANY_THROW(elemFieldRight.bytes().bucket_bytes_right(bucket_elem1));        // Bucket from different mesh
  EXPECT_ANY_THROW(elemFieldRight.const_bytes().bucket_bytes_right(bucket_elem1));  // Bucket from different mesh
}

//------------------------------------------------------------------------------
TEST_F(FieldBytesAccess, host_consistencyCheck_bucketId)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, 1, 1);  // Small Buckets to force creation of many

  auto& elemField      = get_meta().declare_field<int>(stk::topology::ELEM_RANK, "elemField");
  auto& elemFieldLeft  = get_meta().declare_field<int, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK, "elemFieldLeft");
  auto& elemFieldRight = get_meta().declare_field<int, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK, "elemFieldRight");
  stk::mesh::put_field_on_mesh(elemField, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldLeft, get_meta().universal_part(), nullptr);
  stk::mesh::put_field_on_mesh(elemFieldRight, get_meta().universal_part(), nullptr);
  stk::io::fill_mesh("generated:1x1x1", get_bulk());

  EXPECT_ANY_THROW(elemField.bytes().bucket_bytes(1));                   // Bad Bucket Id
  EXPECT_ANY_THROW(elemField.const_bytes().bucket_bytes(1));             // Bad Bucket Id
  EXPECT_ANY_THROW(elemFieldLeft.bytes().bucket_bytes_left(1));          // Bad Bucket Id
  EXPECT_ANY_THROW(elemFieldLeft.const_bytes().bucket_bytes_left(1));    // Bad Bucket Id
  EXPECT_ANY_THROW(elemFieldRight.bytes().bucket_bytes_right(1));        // Bad Bucket Id
  EXPECT_ANY_THROW(elemFieldRight.const_bytes().bucket_bytes_right(1));  // Bad Bucket Id
}

#endif

//==============================================================================

}
