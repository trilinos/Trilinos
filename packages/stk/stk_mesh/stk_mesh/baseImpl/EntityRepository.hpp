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

namespace stk {
namespace mesh {
namespace impl {

class EntityRepository {

public:

    typedef std::vector<std::pair<EntityKey,Entity> > EntityKeyEntityVector;
    typedef EntityKeyEntityVector::const_iterator const_iterator;
    typedef EntityKeyEntityVector::iterator iterator;

    EntityRepository();

    ~EntityRepository();

    Entity get_entity( const EntityKey &key ) const;

    const_entity_iterator begin_rank(EntityRank ent_rank) const
    {
      clear_cache(ent_rank);
      return m_entities[ent_rank].begin();
    }
    const_entity_iterator end_rank(EntityRank ent_rank) const
    {
      clear_cache(ent_rank);
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

    void destroy_entity(EntityKey key, Entity entity);

    void update_num_ranks(unsigned numRanks) {
        m_entities.resize(numRanks);
        m_create_cache.resize(numRanks);
        m_update_cache.resize(numRanks);
        m_destroy_cache.resize(numRanks);
    }

    void clear_all_cache();
    void clear_cache(EntityRank rank) const;

  private:
    void clear_destroyed_entity_cache(EntityRank rank) const;
    void clear_updated_entity_cache(EntityRank rank) const;
    void clear_created_entity_cache(EntityRank rank) const;
    entity_iterator get_from_cache(const EntityKey& key) const;
    std::pair<entity_iterator,bool> add_to_cache(const EntityKey& key);

    mutable std::vector<EntityKeyEntityVector> m_entities;

    mutable std::vector<EntityKeyEntityVector> m_create_cache;
    mutable std::vector<std::vector<std::pair<EntityKey,EntityKey> > > m_update_cache;
    mutable std::vector<std::vector<EntityKey> > m_destroy_cache;
    mutable unsigned m_maxCreateCacheSize;
    unsigned m_maxUpdateCacheSize;

    //disable copy constructor and assignment operator
    EntityRepository(const EntityRepository &);
    EntityRepository & operator =(const EntityRepository &);
};

} // namespace impl

} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_EntityRepository_hpp

