// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_EntityIdPool_h
#define Akri_EntityIdPool_h

#include <stk_mesh/base/Types.hpp>
#include <vector>

namespace stk { namespace mesh { class BulkData; } }

namespace krino {

class EntityIdPool {
public:
  EntityIdPool(stk::mesh::MetaData & meta_data);
  ~EntityIdPool() {}

  // Like stk::mesh::generate_new_ids, this method has lowest ids at front of vector
  static void generate_new_ids(stk::mesh::BulkData & mesh, stk::mesh::EntityRank rank, size_t count, std::vector<stk::mesh::EntityId> & ids, bool assert_32bit_ids, bool make_64bit_ids);

  void reserve(stk::mesh::EntityRank rank, size_t count, bool assert_32bit_ids, bool make_64bit_ids);
  stk::mesh::EntityId get_EntityId(stk::mesh::EntityRank rank);
private:
  static void push_ids_to_64_bit(std::vector<stk::mesh::EntityId> & ids);
  static void check_ids_are_32_bit(stk::mesh::EntityRank rank, std::vector<stk::mesh::EntityId> & ids);

protected:
  stk::mesh::MetaData & my_meta_data;
  std::vector< std::vector<stk::mesh::EntityId> > my_entity_id_pool;
};

} // namespace krino

#endif // Akri_EntityIdPool_h
