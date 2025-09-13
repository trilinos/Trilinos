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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_ENTITYKEY_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_ENTITYKEY_HPP_

#include <iosfwd>                       // for ostream
#include <ostream>
#include <memory>
#include <vector>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace search {
namespace spmd {

struct EntityKeyPair
{
  KOKKOS_FUNCTION
  EntityKeyPair()
  {}

  KOKKOS_FUNCTION
  EntityKeyPair( stk::mesh::Entity entity, stk::mesh::EntityKey key )
    : m_entity(entity)
    , m_key(key)
  {}

  KOKKOS_FUNCTION
  stk::mesh::EntityId   id() const   { return m_key.id(); }

  KOKKOS_FUNCTION
  stk::mesh::EntityRank rank() const { return m_key.rank(); }

  KOKKOS_FUNCTION
  bool is_valid() const { return m_key.is_valid(); }

  operator stk::mesh::Entity() const { return m_entity; }
  operator stk::mesh::EntityKey() const { return m_key; }

  stk::mesh::Entity m_entity;
  stk::mesh::EntityKey m_key;
};

std::ostream & operator << ( std::ostream &os , const EntityKeyPair &key );

inline bool operator < (stk::mesh::EntityKey const& l, EntityKeyPair const& r)
{
  return l < r.m_key;
}
inline bool operator < (EntityKeyPair const& l, stk::mesh::EntityKey const& r)
{
  return l.m_key < r;
}
inline bool operator < (EntityKeyPair const& l, EntityKeyPair const& r)
{
  return l.m_key < r.m_key;
}
inline bool operator == (EntityKeyPair const& l, EntityKeyPair const& r)
{
  return (l.m_key == r.m_key) && (l.m_entity == r.m_entity);
}

EntityKeyPair make_entity_key_pair(const stk::mesh::BulkData& bulk, stk::mesh::EntityKey key);
EntityKeyPair make_entity_key_pair(const stk::mesh::BulkData& bulk, stk::mesh::Entity entity);

EntityKeyPair make_entity_key_pair(const stk::mesh::BulkData* bulk, stk::mesh::EntityKey key);
EntityKeyPair make_entity_key_pair(const stk::mesh::BulkData* bulk, stk::mesh::Entity entity);

EntityKeyPair make_entity_key_pair(const std::shared_ptr<stk::mesh::BulkData>& bulk, stk::mesh::EntityKey key);
EntityKeyPair make_entity_key_pair(const std::shared_ptr<stk::mesh::BulkData>& bulk, stk::mesh::Entity entity);

stk::mesh::EntityVector convert_keys_to_entities(const stk::mesh::BulkData &bulk, const std::vector<EntityKeyPair>& node_keys);

} // namespace spmd
} // namespace search
} // namespace stk




#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_SPMD_ENTITYKEY_HPP_ */
