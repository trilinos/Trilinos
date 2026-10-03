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

#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_topology/topology.hpp>
#include <algorithm>
#include <cstdint>
#include <vector>

namespace {

class GetBucketsFieldDataModCycle : public stk::unit_test_util::MeshFixture {};

// On the cached get_buckets() fast-path with a rank dirty from a pending
// sync_from_partitions, a caller iterating the returned list and reading field data can
// trip FieldBase's in-modification sync guard-rail, which reorganizes the cached list
// mid-iteration and yields duplicate/missing entities. get_buckets() must sync up front.
TEST_F(GetBucketsFieldDataModCycle, cachedGetBucketsSyncsBeforeIteration)
{
  if (stk::parallel_machine_size(get_comm()) != 1) { GTEST_SKIP(); }

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);

  stk::mesh::MetaData& meta = get_meta();

  stk::mesh::Field<uint64_t>& parentIdField =
      meta.declare_field<uint64_t>(stk::topology::ELEMENT_RANK, "parentId");
  stk::mesh::put_field_on_mesh(parentIdField, meta.universal_part(), nullptr);

  // Three element parts so the universal-part selector's cached bucket list has > 1 bucket.
  stk::mesh::Part& partA = meta.declare_part_with_topology("partA", stk::topology::HEX_8);
  stk::mesh::Part& partB = meta.declare_part_with_topology("partB", stk::topology::HEX_8);
  stk::mesh::Part& partC = meta.declare_part_with_topology("partC", stk::topology::HEX_8);

  stk::io::fill_mesh("generated:1x1x6", get_bulk());
  stk::mesh::BulkData& bulk = get_bulk();

  // Spread the 6 elements across the three parts -> three distinct ELEMENT_RANK buckets.
  bulk.modification_begin();
  for (stk::mesh::EntityId id = 1; id <= 6; ++id) {
    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEMENT_RANK, id);
    stk::mesh::Part& target = (id <= 2) ? partA : ((id <= 4) ? partB : partC);
    bulk.change_entity_parts(elem, stk::mesh::PartVector{&target});
  }
  bulk.modification_end();

  // Seed parentId with each element's own id so the field read below returns a known value.
  {
    auto fieldData = parentIdField.data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket* b : bulk.buckets(stk::topology::ELEMENT_RANK)) {
      for (stk::mesh::Entity e : *b) {
        fieldData.entity_values(e)() = static_cast<uint64_t>(bulk.identifier(e));
      }
    }
  }

  stk::mesh::Selector sel = meta.universal_part();

  bulk.modification_begin();

  // Prime the get_buckets cache so the critical call below hits the cached fast-path.
  const stk::mesh::BucketVector& primed = bulk.get_buckets(stk::topology::ELEMENT_RANK, sel);
  ASSERT_EQ(3u, primed.size());

  // Empty a non-last bucket mid-modification so a bucket destroy/renumber is left pending.
  for (stk::mesh::EntityId id = 3; id <= 4; ++id) {
    stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEMENT_RANK, id);
    bulk.change_entity_parts(elem, stk::mesh::PartVector{&partA}, stk::mesh::PartVector{&partB});
  }

  // Ground truth: get_entities reads no field data, so it does not trip the guard-rail sync.
  stk::mesh::EntityVector groundTruthEntities;
  bulk.get_entities(stk::topology::ELEMENT_RANK, sel, groundTruthEntities);
  std::vector<stk::mesh::EntityId> groundTruth;
  for (stk::mesh::Entity e : groundTruthEntities) {
    groundTruth.push_back(bulk.identifier(e));
  }
  std::sort(groundTruth.begin(), groundTruth.end());
  ASSERT_EQ(6u, groundTruth.size());

  // Critical call: cached fast-path while the rank is still dirty.
  const stk::mesh::BucketVector& bkts = bulk.get_buckets(stk::topology::ELEMENT_RANK, sel);

  std::vector<stk::mesh::Bucket*> before(bkts.begin(), bkts.end());

  // Read field data per entity while iterating the cached list; the first read trips the
  // guard-rail sync that (pre-fix) reorganizes 'bkts' out from under this loop.
  std::vector<stk::mesh::EntityId> visited;
  for (stk::mesh::Bucket* b : bkts) {
    for (stk::mesh::Entity e : *b) {
      const uint64_t value = parentIdField.data<stk::mesh::ReadOnly>().entity_values(e)();
      visited.push_back(static_cast<stk::mesh::EntityId>(value));
    }
  }

  std::vector<stk::mesh::Bucket*> after(bkts.begin(), bkts.end());

  EXPECT_EQ(before, after)
      << "get_buckets() cached list was reorganized while a caller was iterating it";

  std::sort(visited.begin(), visited.end());
  EXPECT_EQ(groundTruth, visited)
      << "iteration over the cached bucket list visited duplicate or missing elements";

  bulk.modification_end();
}

}  // namespace
