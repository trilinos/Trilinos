// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_EntityIdPool.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino{

EntityIdPool::EntityIdPool(stk::mesh::MetaData & meta_data)
  : my_meta_data(meta_data),
    my_entity_id_pool(meta_data.entity_rank_count())
{
}

stk::mesh::EntityId
EntityIdPool::get_EntityId(stk::mesh::EntityRank rank)
{
  STK_ThrowAssert(static_cast<unsigned>(rank) < my_entity_id_pool.size());
  STK_ThrowRequireMsg(!my_entity_id_pool[rank].empty(), "EntityIdPool is empty for " << rank << ".");

  stk::mesh::EntityId entity_id = my_entity_id_pool[rank].back();
  my_entity_id_pool[rank].pop_back();
  return entity_id;
}

void
EntityIdPool::reserve(stk::mesh::EntityRank rank, size_t count, bool assert_32bit_ids, bool make_64bit_ids)
{
  auto & ids = my_entity_id_pool[rank];
  generate_new_ids(my_meta_data.mesh_bulk_data(), rank, count, ids, assert_32bit_ids, make_64bit_ids);

  // Reverse order because these are used from the back
  const size_t num = ids.size();
  for (size_t i=0; i<num/2; ++i)
    std::swap(ids[i], ids[num-i-1]);
}

void
EntityIdPool::generate_new_ids(stk::mesh::BulkData & mesh, stk::mesh::EntityRank rank, size_t count, std::vector<stk::mesh::EntityId> & ids, bool assert_32bit_ids, bool make_64bit_ids)
{
  mesh.generate_new_ids( rank, count, ids );
  STK_ThrowAssert(!make_64bit_ids || !assert_32bit_ids);
  if (make_64bit_ids)
  {
    push_ids_to_64_bit(ids);
  }
  if (assert_32bit_ids && (rank == stk::topology::NODE_RANK || rank == stk::topology::ELEMENT_RANK)) // Only worry about node and elements since they are output
  {
    check_ids_are_32_bit(rank, ids);
  }
}

void
EntityIdPool::push_ids_to_64_bit(std::vector<stk::mesh::EntityId> & ids)
{
  const uint64_t max_32_bit_id = std::numeric_limits<uint32_t>::max();
  const bool ids_are_32_bit = !ids.empty() && ids[0] <= max_32_bit_id;
  if (ids_are_32_bit)
  {
    for (auto && id : ids)
    {
      STK_ThrowRequireMsg(id <= max_32_bit_id, "Mixture of ids above and below 32 bit limit not allowed in push_ids_to_64_bit.");
      id += max_32_bit_id;
    }
  }
}

void
EntityIdPool::check_ids_are_32_bit(stk::mesh::EntityRank rank, std::vector<stk::mesh::EntityId> & ids)
{
  static const uint64_t max_allowed_id = std::numeric_limits<uint32_t>::max();
  bool have_bad_id = false;
  for (auto && id : ids)
  {
    if (id > max_allowed_id)
    {
      have_bad_id = true;
      break;
    }
  }
  STK_ThrowRequireMsg(!have_bad_id, "Exhausted valid 32 bit ids for rank " << rank << "!");
}

} // namespace krino
