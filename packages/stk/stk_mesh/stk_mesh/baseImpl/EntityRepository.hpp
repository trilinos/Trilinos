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

#ifndef stk_mesh_baseImpl_EntityRepository_hpp
#define stk_mesh_baseImpl_EntityRepository_hpp

#include <stddef.h>                     // for size_t
#include <map>                          // for map, map<>::value_compare
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <utility>                      // for pair
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Types.hpp"      // for EntityRank

namespace stk {
namespace mesh {
namespace impl {

class EntityRepository {

public:

    typedef std::pair<EntityKey,Entity> EntityKeyEntity;
    typedef std::vector<EntityKeyEntity> EntityKeyEntityVector;
    typedef std::map<EntityKey,Entity> EntityKeyEntityMap;
    typedef EntityKeyEntityMap::const_iterator const_iterator;
    typedef EntityKeyEntityMap::iterator iterator;

    EntityRepository();

    ~EntityRepository();

    Entity get_entity( const EntityKey &key ) const;

    const_entity_iterator begin_rank(EntityRank ent_rank) const
    {
      return m_entities[ent_rank].begin();
    }
    const_entity_iterator end_rank(EntityRank ent_rank) const
    {
      return m_entities[ent_rank].end();
    }

    // Return a pair: the relevant entity, and whether it had to be created
    // or not. If there was already an active entity with the specified key,
    // the second item in the pair will be false;
    // otherwise it will be true (even if the Entity was present
    // but marked as destroyed).
    std::pair<entity_iterator, bool>
    internal_create_entity( const EntityKey & key);

    void update_entity_key(EntityKey new_key, EntityKey old_key, Entity entity);

    void destroy_entity(EntityKey key);

    void update_num_ranks(unsigned numRanks) {
        m_entities.resize(numRanks);
    }

    size_t heap_memory_in_bytes() const;

  private:

    mutable std::vector<EntityKeyEntityMap> m_entities;

    //disable copy constructor and assignment operator
    EntityRepository(const EntityRepository &);
    EntityRepository & operator =(const EntityRepository &);
};

} // namespace impl

} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_EntityRepository_hpp

