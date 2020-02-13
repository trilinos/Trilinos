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

#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stddef.h>                     // for NULL
#include <sstream>                      // for operator<<, basic_ostream, etc
#include <vector>
#include "stk_mesh/base/Entity.hpp"     // for Entity, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc


namespace stk {
namespace mesh {
namespace impl {

EntityRepository::EntityRepository()
 : m_entities(stk::topology::NUM_RANKS)
{
}

EntityRepository::~EntityRepository()
{
}

template<typename VecType>
size_t capacity_in_bytes(const VecType& v)
{
  return sizeof(typename VecType::value_type)*v.capacity();
}

size_t EntityRepository::heap_memory_in_bytes() const
{
    size_t bytes = 0;
    for(auto entityMap : m_entities) { bytes += entityMap.size()*sizeof(EntityKeyEntity); }
    return bytes;
}

std::pair<stk::mesh::entity_iterator ,bool>
EntityRepository::internal_create_entity( const EntityKey & key)
{
  EntityKeyEntityMap& entMap = m_entities[key.rank()];
  return entMap.insert(EntityKeyEntity(key,Entity()));
}

Entity EntityRepository::get_entity(const EntityKey &key) const
{
  const EntityRank rank = key.rank();

  entity_iterator ent = m_entities[rank].find(key);
  if (ent != m_entities[rank].end()) {
    return ent->second;
  }

  return Entity();
}

void EntityRepository::update_entity_key(EntityKey new_key, EntityKey old_key, Entity entity)
{
  const EntityRank rank = new_key.rank();

  m_entities[rank].erase(old_key);
  m_entities[rank].insert(EntityKeyEntity(new_key, entity));
}

void EntityRepository::destroy_entity(EntityKey key)
{ 
  m_entities[key.rank()].erase(key);
} 

} // namespace impl
} // namespace mesh
} // namespace stk

