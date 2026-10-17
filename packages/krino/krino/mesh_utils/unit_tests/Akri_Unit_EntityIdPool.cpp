// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_EntityIdPool.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_StkMeshFixture.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>

#include <limits>
#include <vector>

namespace krino
{

// A committed single-tet mesh is enough to give EntityIdPool a BulkData to
// generate ids from.
class EntityIdPoolFixture : public StkMeshTetFixture
{
protected:
  EntityIdPoolFixture()
  {
    RegularTet meshSpec;
    set_valid_proc_sizes_for_test({1});
    if (is_valid_proc_size_for_test())
      StkMeshTetFixture::build_mesh(meshSpec.nodeLocs, meshSpec.allElementConn, {1});
  }
};

//--------------------------------------------------------------------------
// reserve() fills the pool and get_EntityId() hands ids out from the back;
// because reserve reverses the generated (ascending) ids, the ids come out
// in ascending order.
//--------------------------------------------------------------------------
TEST_F(EntityIdPoolFixture, reserveThenGet_returnsRequestedCountInAscendingOrder)
{
  if (!is_valid_proc_size_for_test())
    return;

  EntityIdPool pool(get_mesh().mesh_meta_data());

  const size_t count = 5;
  pool.reserve(stk::topology::NODE_RANK, count, false, false);

  stk::mesh::EntityId previous = 0;
  for (size_t i = 0; i < count; ++i)
  {
    const stk::mesh::EntityId id = pool.get_EntityId(stk::topology::NODE_RANK);
    EXPECT_GT(id, 0u);
    if (i > 0)
    {
      EXPECT_GT(id, previous) << "ids should be handed out in ascending order";
    }
    previous = id;
  }
}

//--------------------------------------------------------------------------
// Drawing more ids than were reserved must throw rather than return garbage.
//--------------------------------------------------------------------------
TEST_F(EntityIdPoolFixture, getBeyondReserved_throws)
{
  if (!is_valid_proc_size_for_test())
    return;

  EntityIdPool pool(get_mesh().mesh_meta_data());
  pool.reserve(stk::topology::NODE_RANK, 1, false, false);

  EXPECT_NO_THROW(pool.get_EntityId(stk::topology::NODE_RANK));
  EXPECT_ANY_THROW(pool.get_EntityId(stk::topology::NODE_RANK));
}

//--------------------------------------------------------------------------
// The static generate_new_ids honors the make_64bit_ids flag: 32-bit ids are
// pushed above the 32-bit limit.
//--------------------------------------------------------------------------
TEST_F(EntityIdPoolFixture, generateNewIds_make64bit_pushesIdsAbove32BitLimit)
{
  if (!is_valid_proc_size_for_test())
    return;

  const uint64_t max32 = std::numeric_limits<uint32_t>::max();

  std::vector<stk::mesh::EntityId> ids32;
  EntityIdPool::generate_new_ids(get_mesh(), stk::topology::NODE_RANK, 4, ids32, false, false);
  ASSERT_FALSE(ids32.empty());
  for (const auto id : ids32)
    EXPECT_LE(id, max32) << "freshly generated ids should start in the 32-bit range";

  std::vector<stk::mesh::EntityId> ids64;
  EntityIdPool::generate_new_ids(get_mesh(), stk::topology::NODE_RANK, 4, ids64, false, true);
  ASSERT_FALSE(ids64.empty());
  for (const auto id : ids64)
    EXPECT_GT(id, max32) << "make_64bit_ids should push ids above the 32-bit limit";
}

}
