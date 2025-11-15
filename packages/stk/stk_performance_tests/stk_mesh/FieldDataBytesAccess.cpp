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
#include <stk_util/stk_config.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

namespace {

class FieldDataBytesAccess : public stk::unit_test_util::MeshFixture
{
public:
  FieldDataBytesAccess()
    : batchTimer(get_comm())
  { }

protected:
  void declare_element_field()
  {
    m_elementField = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "elementField");
    stk::mesh::put_field_on_mesh(*m_elementField, get_meta().universal_part(), 3, nullptr);
  }

  void declare_element_field_left()
  {
    m_elementFieldLeft = &get_meta().declare_field<double, stk::mesh::Layout::Left>(stk::topology::ELEM_RANK,
                                                                                    "elementFieldLeft");
    stk::mesh::put_field_on_mesh(*m_elementFieldLeft, get_meta().universal_part(), 3, nullptr);
  }

  void declare_element_field_right()
  {
    m_elementFieldRight = &get_meta().declare_field<double, stk::mesh::Layout::Right>(stk::topology::ELEM_RANK,
                                                                                      "elementFieldRight");
    stk::mesh::put_field_on_mesh(*m_elementFieldRight, get_meta().universal_part(), 3, nullptr);
  }

  void fill_multi_block_mesh(unsigned numElemsPerDim)
  {
    stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(numElemsPerDim), get_bulk());
    stk::performance_tests::move_elements_to_other_blocks(get_bulk(), numElemsPerDim);
  }

  template <typename BYTES_FUNCTOR, typename VERIFIER_FUNCTOR>
  void run_test(int numIters, const BYTES_FUNCTOR& bytesFunctor, const VERIFIER_FUNCTOR& verifierFunctor)
  {
    const unsigned NUM_RUNS = 5;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    declare_element_field();
    stk::io::fill_mesh("generated:512x100x100", get_bulk());

    stk::mesh::Selector selector(get_meta().locally_owned_part());

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        bytesFunctor(selector, *m_elementField);
      }
      batchTimer.stop_batch_timer();
    }

    verifierFunctor(selector, *m_elementField);

    batchTimer.print_batch_timing(numIters);
  }

  template <typename BYTES_FUNCTOR, typename VERIFIER_FUNCTOR>
  void run_test_left(int numIters, const BYTES_FUNCTOR& bytesFunctor, const VERIFIER_FUNCTOR& verifierFunctor)
  {
    const unsigned NUM_RUNS = 5;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    declare_element_field_left();
    stk::io::fill_mesh("generated:512x100x100", get_bulk());

    stk::mesh::Selector selector(get_meta().locally_owned_part());

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        bytesFunctor(selector, *m_elementFieldLeft);
      }
      batchTimer.stop_batch_timer();
    }

    verifierFunctor(selector, *m_elementFieldLeft);

    batchTimer.print_batch_timing(numIters);
  }

  template <typename BYTES_FUNCTOR, typename VERIFIER_FUNCTOR>
  void run_test_right(int numIters, const BYTES_FUNCTOR& bytesFunctor, const VERIFIER_FUNCTOR& verifierFunctor)
  {
    const unsigned NUM_RUNS = 5;

    batchTimer.initialize_batch_timer();

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    declare_element_field_right();
    stk::io::fill_mesh("generated:512x100x100", get_bulk());

    stk::mesh::Selector selector(get_meta().locally_owned_part());

    for (unsigned run = 0; run < NUM_RUNS; ++run) {
      batchTimer.start_batch_timer();
      for (int iter = 0; iter < numIters; ++iter) {
        bytesFunctor(selector, *m_elementFieldRight);
      }
      batchTimer.stop_batch_timer();
    }

    verifierFunctor(selector, *m_elementFieldRight);

    batchTimer.print_batch_timing(numIters);
  }

  stk::unit_test_util::BatchTimer batchTimer;
  stk::mesh::Field<double> *m_elementField;
  stk::mesh::Field<double, stk::mesh::Layout::Left> *m_elementFieldLeft;
  stk::mesh::Field<double, stk::mesh::Layout::Right> *m_elementFieldRight;
};

auto verify_initialized_field = [](const stk::mesh::Selector& selector, auto& elementField)
{
  const stk::mesh::BulkData& bulk = elementField.get_mesh();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  std::array<double, 3> initValue { 1.0, 2.0, 3.0 };
  auto fieldData = elementField.template data<>();

  for (const stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      auto entityValues = fieldData.entity_values(elem);
      EXPECT_EQ(entityValues(0_comp), initValue[0]);
      EXPECT_EQ(entityValues(1_comp), initValue[1]);
      EXPECT_EQ(entityValues(2_comp), initValue[2]);
    }
  }
};


#ifndef STK_UNIFIED_MEMORY
class LegacyFieldBytesAccess : public FieldDataBytesAccess {};

auto legacy_initialize_field_entity_bytes = [](const stk::mesh::Selector& selector,
                                               stk::mesh::Field<double>& elementField)
{
  const stk::mesh::BulkData& bulk = elementField.get_mesh();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  std::array<double, 3> initValue { 1.0, 2.0, 3.0 };
  const std::byte* initValueBytes = reinterpret_cast<std::byte*>(initValue.data());

  for (const stk::mesh::Bucket* bucket : buckets) {
    const int bytesPerEntity = stk::mesh::field_bytes_per_entity(elementField, *bucket);
    for (stk::mesh::Entity elem : *bucket) {
      std::byte* entityBytes = reinterpret_cast<std::byte*>(stk::mesh::field_data(elementField, elem));
      std::memcpy(entityBytes, initValueBytes, bytesPerEntity);
    }
  }
};


auto legacy_initialize_field_bucket_bytes = [](const stk::mesh::Selector& selector,
                                               stk::mesh::Field<double>& elementField)
{
  const stk::mesh::BulkData& bulk = elementField.get_mesh();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  std::array<double, 3> initValue { 1.0, 2.0, 3.0 };
  const std::byte* initValueBytes = reinterpret_cast<std::byte*>(initValue.data());

  for (const stk::mesh::Bucket* bucket : buckets) {
    auto& fieldMetaData = elementField.get_meta_data_for_field()[bucket->bucket_id()];
    for (int bucketOrd = 0; bucketOrd < static_cast<int>(bucket->size()); ++bucketOrd) {
      std::memcpy(fieldMetaData.m_data + bucketOrd*fieldMetaData.m_bytesPerEntity, initValueBytes,
                  fieldMetaData.m_bytesPerEntity);
    }
  }
};

//------------------------------------------------------------------------------
TEST_F(LegacyFieldBytesAccess, initializeFieldEntityBytes)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_test(100, legacy_initialize_field_entity_bytes, verify_initialized_field);
}

TEST_F(LegacyFieldBytesAccess, initializeFieldBucketBytes)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_test(100, legacy_initialize_field_bucket_bytes, verify_initialized_field);
}
#endif  // For not STK_UNIFIED_MEMORY


auto initialize_field_entity_bytes = [](const stk::mesh::Selector& selector,
                                        stk::mesh::Field<double>& elementField)
{
  const stk::mesh::BulkData& bulk = elementField.get_mesh();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  std::array<double, 3> initValue { 1.0, 2.0, 3.0 };
  const std::byte* initValueBytes = reinterpret_cast<std::byte*>(initValue.data());

  auto& fieldDataBytes = elementField.data_bytes<std::byte>();

  for (const stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      auto entityBytes = fieldDataBytes.entity_bytes(elem);
      for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
        entityBytes(byte) = initValueBytes[byte];
      }
    }
  }
};

auto initialize_field_entity_bytes_left = [](const stk::mesh::Selector& selector,
                                             stk::mesh::Field<double, stk::mesh::Layout::Left>& elementField)
{
  const stk::mesh::BulkData& bulk = elementField.get_mesh();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  std::array<double, 3> initValue { 1.0, 2.0, 3.0 };
  const std::byte* initValueBytes = reinterpret_cast<std::byte*>(initValue.data());

  auto& fieldDataBytes = elementField.data_bytes<std::byte>();

  for (const stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      auto entityBytes = fieldDataBytes.entity_bytes<stk::mesh::Layout::Left>(elem);
      for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
        entityBytes(byte) = initValueBytes[byte];
      }
    }
  }
};

auto initialize_field_entity_bytes_right = [](const stk::mesh::Selector& selector,
                                              stk::mesh::Field<double, stk::mesh::Layout::Right>& elementField)
{
  const stk::mesh::BulkData& bulk = elementField.get_mesh();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  std::array<double, 3> initValue { 1.0, 2.0, 3.0 };
  const std::byte* initValueBytes = reinterpret_cast<std::byte*>(initValue.data());

  auto& fieldDataBytes = elementField.data_bytes<std::byte>();

  for (const stk::mesh::Bucket* bucket : buckets) {
    for (stk::mesh::Entity elem : *bucket) {
      auto entityBytes = fieldDataBytes.entity_bytes<stk::mesh::Layout::Right>(elem);
      for (stk::mesh::ByteIdx byte : entityBytes.bytes()) {
        entityBytes(byte) = initValueBytes[byte];
      }
    }
  }
};

auto initialize_field_bucket_bytes = [](const stk::mesh::Selector& selector,
                                        stk::mesh::Field<double>& elementField)
{
  const stk::mesh::BulkData& bulk = elementField.get_mesh();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  std::array<double, 3> initValue { 1.0, 2.0, 3.0 };
  const std::byte* initValueBytes = reinterpret_cast<std::byte*>(initValue.data());

  auto& fieldDataBytes = elementField.data_bytes<std::byte>();

  for (const stk::mesh::Bucket* bucket : buckets) {
    auto bucketBytes = fieldDataBytes.bucket_bytes(*bucket);
    for (stk::mesh::EntityIdx elem : bucket->entities()) {
      for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
        bucketBytes(elem, byte) = initValueBytes[byte];
      }
    }
  }
};

auto initialize_field_bucket_bytes_left = [](const stk::mesh::Selector& selector,
                                             stk::mesh::Field<double, stk::mesh::Layout::Left>& elementField)
{
  const stk::mesh::BulkData& bulk = elementField.get_mesh();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  std::array<double, 3> initValue { 1.0, 2.0, 3.0 };
  const std::byte* initValueBytes = reinterpret_cast<std::byte*>(initValue.data());

  auto& fieldDataBytes = elementField.data_bytes<std::byte>();

  for (const stk::mesh::Bucket* bucket : buckets) {
    auto bucketBytes = fieldDataBytes.bucket_bytes<stk::mesh::Layout::Left>(*bucket);
    for (stk::mesh::EntityIdx elem : bucket->entities()) {
      for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
        bucketBytes(elem, byte) = initValueBytes[byte];
      }
    }
  }
};

auto initialize_field_bucket_bytes_right = [](const stk::mesh::Selector& selector,
                                              stk::mesh::Field<double, stk::mesh::Layout::Right>& elementField)
{
  const stk::mesh::BulkData& bulk = elementField.get_mesh();
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  std::array<double, 3> initValue { 1.0, 2.0, 3.0 };
  const std::byte* initValueBytes = reinterpret_cast<std::byte*>(initValue.data());

  auto& fieldDataBytes = elementField.data_bytes<std::byte>();

  for (const stk::mesh::Bucket* bucket : buckets) {
    auto bucketBytes = fieldDataBytes.bucket_bytes<stk::mesh::Layout::Right>(*bucket);
    for (stk::mesh::EntityIdx elem : bucket->entities()) {
      for (stk::mesh::ByteIdx byte : bucketBytes.bytes()) {
        bucketBytes(elem, byte) = initValueBytes[byte];
      }
    }
  }
};

//------------------------------------------------------------------------------
TEST_F(FieldDataBytesAccess, initializeFieldEntityBytes)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_test(100, initialize_field_entity_bytes, verify_initialized_field);
}

TEST_F(FieldDataBytesAccess, initializeFieldEntityBytesLeft)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_test_left(100, initialize_field_entity_bytes_left, verify_initialized_field);
}

TEST_F(FieldDataBytesAccess, initializeFieldEntityBytesRight)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_test_right(100, initialize_field_entity_bytes_right, verify_initialized_field);
}


//------------------------------------------------------------------------------
TEST_F(FieldDataBytesAccess, initializeFieldBucketBytes)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_test(100, initialize_field_bucket_bytes, verify_initialized_field);
}

TEST_F(FieldDataBytesAccess, initializeFieldBucketBytesLeft)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_test_left(100, initialize_field_bucket_bytes_left, verify_initialized_field);
}

TEST_F(FieldDataBytesAccess, initializeFieldBucketBytesRight)
{
  if (get_parallel_size() != 1) GTEST_SKIP();

  run_test_right(100, initialize_field_bucket_bytes_right, verify_initialized_field);
}

//------------------------------------------------------------------------------
}
