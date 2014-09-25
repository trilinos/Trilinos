// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef stk_mesh_baseImpl_EntityRepository_hpp
#define stk_mesh_baseImpl_EntityRepository_hpp

#include <stddef.h>                     // for size_t
#include <map>                          // for map, map<>::value_compare
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <utility>                      // for pair
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Types.hpp"      // for EntityRank
namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace mesh {
namespace impl {

class EntityRepository {

public:

  typedef std::map<EntityKey,Entity> EntityMap;

    typedef EntityMap::const_iterator const_iterator;
    typedef EntityMap::iterator iterator;

    EntityRepository(BulkData &mesh)
      : m_mesh(mesh), m_entities() {}

    ~EntityRepository();

    Entity get_entity( const EntityKey &key ) const;

    const_iterator begin() const { return m_entities.begin(); }
    const_iterator end() const { return m_entities.end(); }

    const_iterator begin_rank(EntityRank ent_rank) const { return m_entities.lower_bound(EntityKey(ent_rank, 0)); }
    const_iterator end_rank(EntityRank ent_rank) const { return m_entities.upper_bound(EntityKey(static_cast<EntityRank>(ent_rank+1), 0)); }

    // Return a pair: the relevant entity, and whether it had to be created
    // or not. If there was already an active entity with the specified key, the second item in the pair
    // will be false; otherwise it will be true (even if the Entity was present
    // but marked as destroyed).
    std::pair<Entity ,bool>
    internal_create_entity( const EntityKey & key, size_t preferred_offset = 0 );

    void change_entity_bucket( Bucket & b, Entity e, unsigned ordinal);

    void update_entity_key(EntityKey new_key, EntityKey old_key, Entity entity);

    void destroy_entity(EntityKey key, Entity entity);

  private:
    void internal_expunge_entity( EntityMap::iterator i);

    // Entity internal_allocate_entity(EntityKey entity_key);
    Entity allocate_entity();

    BulkData &m_mesh;
    EntityMap m_entities;

    //disable copy constructor and assignment operator
    EntityRepository(const EntityRepository &);
    EntityRepository & operator =(const EntityRepository &);
};

inline
void EntityRepository::destroy_entity(EntityKey key, Entity entity)
{
  m_entities.erase(key);
}

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_EntityRepository_hpp
