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

#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for BucketVector, PartVector, etc
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <string>                       // for string, operator==, etc

namespace {

class FieldQueryFunctions : public stk::unit_test_util::MeshFixture
{
public:
  FieldQueryFunctions()
    : MeshFixture(),
      m_doubleField(nullptr),
      m_intField(nullptr),
      m_block1Part(nullptr),
      m_block2Part(nullptr)
  {}

  void setup_mesh_scalar_fields_one_block() {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_block1Part = &get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    m_doubleField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "double_field");
    m_intField = &get_meta().declare_field<int>(stk::topology::ELEMENT_RANK, "int_field");
    stk::mesh::put_field_on_mesh(*m_doubleField, *m_block1Part, nullptr);
    stk::mesh::put_field_on_mesh(*m_intField, *m_block1Part, nullptr);

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n";

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_mesh_vector_fields_one_block() {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_block1Part = &get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    m_doubleField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "double_field");
    m_intField = &get_meta().declare_field<int>(stk::topology::ELEMENT_RANK, "int_field");
    stk::mesh::put_field_on_mesh(*m_doubleField, *m_block1Part, 3, nullptr);
    stk::mesh::put_field_on_mesh(*m_intField, *m_block1Part, 3, nullptr);

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n";

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_mesh_vector_fields_four_copies_one_block() {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_block1Part = &get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    m_doubleField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "double_field");
    m_intField = &get_meta().declare_field<int>(stk::topology::ELEMENT_RANK, "int_field");
    stk::mesh::put_field_on_mesh(*m_doubleField, *m_block1Part, 3, 4, nullptr);
    stk::mesh::put_field_on_mesh(*m_intField, *m_block1Part, 3, 4, nullptr);

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "0,2,HEX_8,5,6,7,8,9,10,11,12,block_1\n";

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_mesh_vector_fields_not_on_all_blocks() {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_block1Part = &get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    m_block2Part = &get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    m_doubleField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "double_field");
    m_intField = &get_meta().declare_field<int>(stk::topology::ELEMENT_RANK, "int_field");
    stk::mesh::put_field_on_mesh(*m_doubleField, *m_block1Part, 3, nullptr);
    stk::mesh::put_field_on_mesh(*m_intField, *m_block2Part, 3, nullptr);

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2\n";

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  void setup_mesh_vector_fields_variable_num_copies() {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

    m_block1Part = &get_meta().declare_part_with_topology("block_1", stk::topology::HEX_8);
    m_block2Part = &get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    m_doubleField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "double_field");
    m_intField = &get_meta().declare_field<int>(stk::topology::ELEMENT_RANK, "int_field");
    stk::mesh::put_field_on_mesh(*m_doubleField, *m_block1Part, 3, 1, nullptr);
    stk::mesh::put_field_on_mesh(*m_doubleField, *m_block2Part, 3, 2, nullptr);
    stk::mesh::put_field_on_mesh(*m_intField, *m_block1Part, 3, 4, nullptr);
    stk::mesh::put_field_on_mesh(*m_intField, *m_block2Part, 3, 1, nullptr);

    const std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                                 "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2\n";

    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);
  }

  stk::mesh::Field<double> * m_doubleField;
  stk::mesh::Field<int> * m_intField;
  stk::mesh::Part * m_block1Part;
  stk::mesh::Part * m_block2Part;
};


TEST_F(FieldQueryFunctions, fieldBytesPerEntity_scalarFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_scalar_fields_one_block();

  const unsigned expectedDoubleBytes = sizeof(double);
  const unsigned expectedIntBytes    = sizeof(int);

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, *bucket), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, *bucket), expectedIntBytes);

    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, bucket->bucket_id()), expectedIntBytes);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, element), expectedDoubleBytes);
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, element), expectedIntBytes);
    }
  }
}

TEST_F(FieldQueryFunctions, fieldBytesPerEntity_vectorFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_one_block();

  const unsigned expectedDoubleBytes = 3 * sizeof(double);
  const unsigned expectedIntBytes    = 3 * sizeof(int);

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, *bucket), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, *bucket), expectedIntBytes);

    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, bucket->bucket_id()), expectedIntBytes);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, element), expectedDoubleBytes);
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, element), expectedIntBytes);
    }
  }
}

TEST_F(FieldQueryFunctions, fieldBytesPerEntity_vectorFieldsFourCopies_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_four_copies_one_block();

  const unsigned expectedDoubleBytes = 3 * 4 * sizeof(double);
  const unsigned expectedIntBytes    = 3 * 4 * sizeof(int);

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, *bucket), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, *bucket), expectedIntBytes);

    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, bucket->bucket_id()), expectedIntBytes);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, element), expectedDoubleBytes);
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, element), expectedIntBytes);
    }
  }
}

TEST_F(FieldQueryFunctions, fieldBytesPerEntity_vectorFields_notOnAllBlocks)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_not_on_all_blocks();

  unsigned expectedDoubleBytes = 0;
  unsigned expectedIntBytes    = 0;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    if (bucket->member(*m_block1Part)) {
      expectedDoubleBytes = 3 * sizeof(double);
      expectedIntBytes    = 0;
    }
    else if (bucket->member(*m_block2Part)) {
      expectedDoubleBytes = 0;
      expectedIntBytes    = 3 * sizeof(int);
    }

    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, *bucket), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, *bucket), expectedIntBytes);

    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, bucket->bucket_id()), expectedIntBytes);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, element), expectedDoubleBytes);
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, element), expectedIntBytes);
    }
  }
}

TEST_F(FieldQueryFunctions, fieldBytesPerEntity_vectorFields_variableNumCopies)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_variable_num_copies();

  unsigned expectedDoubleBytes = 0;
  unsigned expectedIntBytes    = 0;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    if (bucket->member(*m_block1Part)) {
      expectedDoubleBytes = 3 * 1 * sizeof(double);
      expectedIntBytes    = 3 * 4 * sizeof(int);
    }
    else if (bucket->member(*m_block2Part)) {
      expectedDoubleBytes = 3 * 2 * sizeof(double);
      expectedIntBytes    = 3 * 1 * sizeof(int);
    }

    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, *bucket), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, *bucket), expectedIntBytes);

    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleBytes);
    EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, bucket->bucket_id()), expectedIntBytes);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_doubleField, element), expectedDoubleBytes);
      EXPECT_EQ(stk::mesh::field_bytes_per_entity(*m_intField, element), expectedIntBytes);
    }
  }
}


TEST_F(FieldQueryFunctions, isMatchingRank_scalarFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_scalar_fields_one_block();

  EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, stk::topology::ELEMENT_RANK), true);
  EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, stk::topology::ELEMENT_RANK), true);

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, *bucket), true);
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, *bucket), true);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, element), true);
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, element), true);
    }
  }
}

TEST_F(FieldQueryFunctions, isMatchingRank_vectorFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_one_block();

  EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, stk::topology::ELEMENT_RANK), true);
  EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, stk::topology::ELEMENT_RANK), true);

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, *bucket), true);
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, *bucket), true);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, element), true);
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, element), true);
    }
  }
}

TEST_F(FieldQueryFunctions, isMatchingRank_vectorFieldsFourCopies_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_four_copies_one_block();

  EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, stk::topology::ELEMENT_RANK), true);
  EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, stk::topology::ELEMENT_RANK), true);

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, *bucket), true);
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, *bucket), true);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, element), true);
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, element), true);
    }
  }
}

TEST_F(FieldQueryFunctions, isMatchingRank_vectorFields_notOnAllBlocks)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_not_on_all_blocks();

  EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, stk::topology::ELEMENT_RANK), true);
  EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, stk::topology::ELEMENT_RANK), true);

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, *bucket), true);
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, *bucket), true);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, element), true);
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, element), true);
    }
  }
}

TEST_F(FieldQueryFunctions, isMatchingRank_vectorFields_variableNumCopies)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_variable_num_copies();

  EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, stk::topology::ELEMENT_RANK), true);
  EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, stk::topology::ELEMENT_RANK), true);

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, *bucket), true);
    EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, *bucket), true);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_doubleField, element), true);
      EXPECT_EQ(stk::mesh::is_matching_rank(*m_intField, element), true);
    }
  }
}


TEST_F(FieldQueryFunctions, fieldScalarsPerEntity_scalarFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_scalar_fields_one_block();

  const unsigned expectedSize = 1;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, *bucket), expectedSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, *bucket), expectedSize);

    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, bucket->bucket_id()), expectedSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, bucket->bucket_id()), expectedSize);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, element), expectedSize);
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, element), expectedSize);
    }
  }
}

TEST_F(FieldQueryFunctions, fieldScalarsPerEntity_vectorFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_one_block();

  const unsigned expectedSize = 3;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, *bucket), expectedSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, *bucket), expectedSize);

    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, bucket->bucket_id()), expectedSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, bucket->bucket_id()), expectedSize);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, element), expectedSize);
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, element), expectedSize);
    }
  }
}

TEST_F(FieldQueryFunctions, fieldScalarsPerEntity_vectorFieldsFourCopies_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_four_copies_one_block();

  const unsigned expectedSize = 3 * 4;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, *bucket), expectedSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, *bucket), expectedSize);

    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, bucket->bucket_id()), expectedSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, bucket->bucket_id()), expectedSize);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, element), expectedSize);
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, element), expectedSize);
    }
  }
}

TEST_F(FieldQueryFunctions, fieldScalarsPerEntity_vectorFields_notOnAllBlocks)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_not_on_all_blocks();

  unsigned expectedDoubleSize = 0;
  unsigned expectedIntSize    = 0;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    if (bucket->member(*m_block1Part)) {
      expectedDoubleSize = 3;
      expectedIntSize    = 0;
    }
    else if (bucket->member(*m_block2Part)) {
      expectedDoubleSize = 0;
      expectedIntSize    = 3;
    }

    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, *bucket), expectedDoubleSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, *bucket), expectedIntSize);

    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, bucket->bucket_id()), expectedIntSize);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, element), expectedDoubleSize);
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, element), expectedIntSize);
    }
  }
}

TEST_F(FieldQueryFunctions, fieldScalarsPerEntity_vectorFields_variableNumCopies)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_variable_num_copies();

  unsigned expectedDoubleSize = 0;
  unsigned expectedIntSize    = 0;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    if (bucket->member(*m_block1Part)) {
      expectedDoubleSize = 3 * 1;
      expectedIntSize    = 3 * 4;
    }
    else if (bucket->member(*m_block2Part)) {
      expectedDoubleSize = 3 * 2;
      expectedIntSize    = 3 * 1;
    }

    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, *bucket), expectedDoubleSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, *bucket), expectedIntSize);

    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleSize);
    EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, bucket->bucket_id()), expectedIntSize);

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_doubleField, element), expectedDoubleSize);
      EXPECT_EQ(stk::mesh::field_scalars_per_entity(*m_intField, element), expectedIntSize);
    }
  }
}


TEST_F(FieldQueryFunctions, fieldExtentPerEntity_scalarFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_scalar_fields_one_block();

  const unsigned expectedSize = 1;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, *bucket), expectedSize);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, *bucket), expectedSize);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, *bucket), expectedSize);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, *bucket), expectedSize);

      for (unsigned i = 0; i < 5; ++i) {
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, *bucket), expectedSize);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, *bucket), expectedSize);
      }
    }

    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, bucket->bucket_id()), expectedSize);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, bucket->bucket_id()), expectedSize);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, bucket->bucket_id()), expectedSize);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, bucket->bucket_id()), expectedSize);

      for (unsigned i = 0; i < 5; ++i) {
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, bucket->bucket_id()), expectedSize);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, bucket->bucket_id()), expectedSize);
      }
    }

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, element), expectedSize);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, element), expectedSize);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, element), expectedSize);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, element), expectedSize);

      for (unsigned i = 0; i < 5; ++i) {
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, element), expectedSize);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, element), expectedSize);
      }
    }
  }
}

TEST_F(FieldQueryFunctions, fieldExtentPerEntity_vectorFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_one_block();

  const unsigned expectedFirstExtent      = 3;
  const unsigned expectedRemainingExtents = 1;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, *bucket), expectedFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, *bucket), expectedFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, *bucket), expectedRemainingExtents);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, *bucket), expectedRemainingExtents);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expected = (i == 0) ? expectedFirstExtent : expectedRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, *bucket), expected);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, *bucket), expected);
      }
    }

    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, bucket->bucket_id()), expectedFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, bucket->bucket_id()), expectedFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, bucket->bucket_id()), expectedRemainingExtents);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, bucket->bucket_id()), expectedRemainingExtents);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expected = (i == 0) ? expectedFirstExtent : expectedRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, bucket->bucket_id()), expected);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, bucket->bucket_id()), expected);
      }
    }

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, element), expectedFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, element), expectedFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, element), expectedRemainingExtents);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, element), expectedRemainingExtents);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expected = (i == 0) ? expectedFirstExtent : expectedRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, element), expected);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, element), expected);
      }
    }
  }
}

TEST_F(FieldQueryFunctions, fieldExtentPerEntity_vectorFieldsFourCopies_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_four_copies_one_block();

  const unsigned expectedFirstExtent      = 3;
  const unsigned expectedSecondExtent     = 4;
  const unsigned expectedRemainingExtents = 1;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, *bucket), expectedFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, *bucket), expectedFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, *bucket), expectedSecondExtent);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, *bucket), expectedSecondExtent);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expected = (i == 0) ? expectedFirstExtent
                                           : (i == 1) ? expectedSecondExtent : expectedRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, *bucket), expected);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, *bucket), expected);
      }
    }

    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, bucket->bucket_id()), expectedFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, bucket->bucket_id()), expectedFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, bucket->bucket_id()), expectedSecondExtent);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, bucket->bucket_id()), expectedSecondExtent);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expected = (i == 0) ? expectedFirstExtent
                                           : (i == 1) ? expectedSecondExtent : expectedRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, bucket->bucket_id()), expected);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, bucket->bucket_id()), expected);
      }
    }

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, element), expectedFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, element), expectedFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, element), expectedSecondExtent);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, element), expectedSecondExtent);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expected = (i == 0) ? expectedFirstExtent
                                           : (i == 1) ? expectedSecondExtent : expectedRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, element), expected);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, element), expected);
      }
    }
  }
}

TEST_F(FieldQueryFunctions, fieldExtentPerEntity_vectorFields_notOnAllBlocks)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_not_on_all_blocks();

  unsigned expectedDoubleFirstExtent      = 0;
  unsigned expectedIntFirstExtent         = 0;
  unsigned expectedDoubleSecondExtent     = 0;
  unsigned expectedIntSecondExtent        = 0;
  unsigned expectedDoubleRemainingExtents = 0;
  unsigned expectedIntRemainingExtents    = 0;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    if (bucket->member(*m_block1Part)) {
      expectedDoubleFirstExtent      = 3;
      expectedIntFirstExtent         = 0;
      expectedDoubleSecondExtent     = 1;
      expectedIntSecondExtent        = 0;
      expectedDoubleRemainingExtents = 1;
      expectedIntRemainingExtents    = 0;
    }
    else if (bucket->member(*m_block2Part)) {
      expectedDoubleFirstExtent      = 0;
      expectedIntFirstExtent         = 3;
      expectedDoubleSecondExtent     = 0;
      expectedIntSecondExtent        = 1;
      expectedDoubleRemainingExtents = 0;
      expectedIntRemainingExtents    = 1;
    }

    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, *bucket), expectedDoubleFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, *bucket), expectedIntFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, *bucket), expectedDoubleSecondExtent);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, *bucket), expectedIntSecondExtent);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expectedDouble = (i == 0) ? expectedDoubleFirstExtent
                                                 : (i == 1) ? expectedDoubleSecondExtent : expectedDoubleRemainingExtents;
        const unsigned expectedInt = (i == 0) ? expectedIntFirstExtent
                                              : (i == 1) ? expectedIntSecondExtent : expectedIntRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, *bucket), expectedDouble);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, *bucket), expectedInt);
      }
    }

    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, bucket->bucket_id()), expectedIntFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleSecondExtent);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, bucket->bucket_id()), expectedIntSecondExtent);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expectedDouble = (i == 0) ? expectedDoubleFirstExtent
                                                 : (i == 1) ? expectedDoubleSecondExtent : expectedDoubleRemainingExtents;
        const unsigned expectedInt = (i == 0) ? expectedIntFirstExtent
                                              : (i == 1) ? expectedIntSecondExtent : expectedIntRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, bucket->bucket_id()), expectedDouble);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, bucket->bucket_id()), expectedInt);
      }
    }

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, element), expectedDoubleFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, element), expectedIntFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, element), expectedDoubleSecondExtent);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, element), expectedIntSecondExtent);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expectedDouble = (i == 0) ? expectedDoubleFirstExtent
                                                 : (i == 1) ? expectedDoubleSecondExtent : expectedDoubleRemainingExtents;
        const unsigned expectedInt = (i == 0) ? expectedIntFirstExtent
                                              : (i == 1) ? expectedIntSecondExtent : expectedIntRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, element), expectedDouble);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, element), expectedInt);
      }
    }
  }
}

TEST_F(FieldQueryFunctions, fieldExtentPerEntity_vectorFields_variableNumCopies)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_variable_num_copies();

  unsigned expectedDoubleFirstExtent      = 0;
  unsigned expectedIntFirstExtent         = 0;
  unsigned expectedDoubleSecondExtent     = 0;
  unsigned expectedIntSecondExtent        = 0;
  unsigned expectedDoubleRemainingExtents = 0;
  unsigned expectedIntRemainingExtents    = 0;

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK))
  {
    if (bucket->member(*m_block1Part)) {
      expectedDoubleFirstExtent      = 3;
      expectedIntFirstExtent         = 3;
      expectedDoubleSecondExtent     = 1;
      expectedIntSecondExtent        = 4;
      expectedDoubleRemainingExtents = 1;
      expectedIntRemainingExtents    = 1;
    }
    else if (bucket->member(*m_block2Part)) {
      expectedDoubleFirstExtent      = 3;
      expectedIntFirstExtent         = 3;
      expectedDoubleSecondExtent     = 2;
      expectedIntSecondExtent        = 1;
      expectedDoubleRemainingExtents = 1;
      expectedIntRemainingExtents    = 1;
    }

    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, *bucket), expectedDoubleFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, *bucket), expectedIntFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, *bucket), expectedDoubleSecondExtent);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, *bucket), expectedIntSecondExtent);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expectedDouble = (i == 0) ? expectedDoubleFirstExtent
                                                 : (i == 1) ? expectedDoubleSecondExtent : expectedDoubleRemainingExtents;
        const unsigned expectedInt = (i == 0) ? expectedIntFirstExtent
                                              : (i == 1) ? expectedIntSecondExtent : expectedIntRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, *bucket), expectedDouble);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, *bucket), expectedInt);
      }
    }

    {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, bucket->bucket_id()), expectedIntFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, bucket->bucket_id()), expectedDoubleSecondExtent);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, bucket->bucket_id()), expectedIntSecondExtent);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expectedDouble = (i == 0) ? expectedDoubleFirstExtent
                                                 : (i == 1) ? expectedDoubleSecondExtent : expectedDoubleRemainingExtents;
        const unsigned expectedInt = (i == 0) ? expectedIntFirstExtent
                                              : (i == 1) ? expectedIntSecondExtent : expectedIntRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, bucket->bucket_id()), expectedDouble);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, bucket->bucket_id()), expectedInt);
      }
    }

    for (const stk::mesh::Entity element : *bucket) {
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_doubleField, element), expectedDoubleFirstExtent);
      EXPECT_EQ(stk::mesh::field_extent0_per_entity(*m_intField, element), expectedIntFirstExtent);

      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_doubleField, element), expectedDoubleSecondExtent);
      EXPECT_EQ(stk::mesh::field_extent1_per_entity(*m_intField, element), expectedIntSecondExtent);

      for (unsigned i = 0; i < 5; ++i) {
        const unsigned expectedDouble = (i == 0) ? expectedDoubleFirstExtent
                                                 : (i == 1) ? expectedDoubleSecondExtent : expectedDoubleRemainingExtents;
        const unsigned expectedInt = (i == 0) ? expectedIntFirstExtent
                                              : (i == 1) ? expectedIntSecondExtent : expectedIntRemainingExtents;
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_doubleField, i, element), expectedDouble);
        EXPECT_EQ(stk::mesh::field_extent_per_entity(*m_intField, i, element), expectedInt);
      }
    }
  }
}


TEST_F(FieldQueryFunctions, fieldIsAllocatedForBucket_scalarFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_scalar_fields_one_block();

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK)) {
    EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_doubleField, *bucket), true);
    EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_intField, *bucket), true);
  }
}

TEST_F(FieldQueryFunctions, fieldIsAllocatedForBucket_vectorFields_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_one_block();

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK)) {
    EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_doubleField, *bucket), true);
    EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_intField, *bucket), true);
  }
}

TEST_F(FieldQueryFunctions, fieldIsAllocatedForBucket_vectorFieldsFourCopies_singleBlock)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_four_copies_one_block();

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK)) {
    EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_doubleField, *bucket), true);
    EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_intField, *bucket), true);
  }
}

TEST_F(FieldQueryFunctions, fieldIsAllocatedForBucket_vectorFields_notOnAllBlocks)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_not_on_all_blocks();

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK)) {
    if (bucket->member(*m_block1Part)) {
      EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_doubleField, *bucket), true);
      EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_intField, *bucket), false);
    }
    else if (bucket->member(*m_block2Part)) {
      EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_doubleField, *bucket), false);
      EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_intField, *bucket), true);
    }
  }
}

TEST_F(FieldQueryFunctions, fieldIsAllocatedForBucket_vectorFields_variableNumCopies)
{
  if (get_parallel_size() != 1) return;

  setup_mesh_vector_fields_variable_num_copies();

  for (const stk::mesh::Bucket * bucket : get_bulk().buckets(stk::topology::ELEMENT_RANK)) {
    EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_doubleField, *bucket), true);
    EXPECT_EQ(stk::mesh::field_is_allocated_for_bucket(*m_intField, *bucket), true);
  }
}

}
