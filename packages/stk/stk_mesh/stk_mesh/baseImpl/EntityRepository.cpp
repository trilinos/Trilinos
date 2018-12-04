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

#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stddef.h>                     // for NULL
#include <sstream>                      // for operator<<, basic_ostream, etc
#include <stk_mesh/base/Bucket.hpp>     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc


namespace stk {
namespace mesh {
namespace impl {

struct EntityKeyEntityLess {
inline bool operator()(const std::pair<EntityKey,Entity>& key_ent_pair, const EntityKey& key) const
{
  return key_ent_pair.first.m_value < key.m_value;
}
};

struct match_EntityKey {
  match_EntityKey(const EntityKey& key)
  : m_key(key)
  {}
  bool operator()(const std::pair<EntityKey,Entity>& item) const
  { return item.first == m_key; }

  bool operator==(const std::pair<EntityKey,Entity>& item) const
  { return item.first == m_key; }

private:
  const EntityKey& m_key;
};

struct ToBeDestroyed {
  ToBeDestroyed(const std::vector<EntityKey>& destroyKeys)
  : m_keys(destroyKeys)
  {}

  bool operator()(const std::pair<EntityKey,Entity>& item) const
  { return std::binary_search(m_keys.begin(), m_keys.end(), item.first); }

private:
  const std::vector<EntityKey>& m_keys;
};

EntityRepository::EntityRepository()
 : m_entities(stk::topology::NUM_RANKS),
   m_create_cache(stk::topology::NUM_RANKS),
   m_update_cache(stk::topology::NUM_RANKS),
   m_destroy_cache(stk::topology::NUM_RANKS),
   m_maxCreateCacheSize(512),
   m_maxUpdateCacheSize(4096)
{
}

EntityRepository::~EntityRepository()
{
}

void EntityRepository::clear_all_cache()
{
  for(EntityRank rank=stk::topology::BEGIN_RANK; rank<m_create_cache.size(); ++rank) {
    clear_cache(rank);
  }
}

void EntityRepository::clear_destroyed_entity_cache(EntityRank rank) const
{
  if (!m_destroy_cache[rank].empty()) {
    std::vector<EntityKey>& destroy = m_destroy_cache[rank];
    std::sort(destroy.begin(), destroy.end());
    EntityKeyEntityVector& entities = m_entities[rank];
    EntityKeyEntityVector::iterator newEnd = std::remove_if(entities.begin(), entities.end(), ToBeDestroyed(destroy));
    entities.resize(newEnd-entities.begin());
    
    destroy.clear();
  }
}

void EntityRepository::clear_updated_entity_cache(EntityRank rank) const
{
  if (!m_update_cache[rank].empty()) {
    std::vector<std::pair<EntityKey,EntityKey> >& update = m_update_cache[rank];
    std::sort(update.begin(), update.end());
    EntityKeyEntityVector::iterator iter = m_entities[rank].begin();
    EntityKeyEntityVector::iterator end = m_entities[rank].end();
    for(const std::pair<EntityKey,EntityKey>& oldnew : update) {
      EntityKeyEntityVector::iterator thisIter = std::lower_bound(iter, end, oldnew.first, EntityKeyEntityLess());
      if (thisIter != end && thisIter->first == oldnew.first) {
        thisIter->first = oldnew.second;
        iter = thisIter+1;
      }
    }
    m_update_cache[rank].clear();
    std::sort(m_entities[rank].begin(), m_entities[rank].end());
  }
}

void EntityRepository::clear_created_entity_cache(EntityRank rank) const
{
  if (!m_create_cache[rank].empty()) {
    std::sort(m_create_cache[rank].begin(), m_create_cache[rank].end());
    unsigned numOld = m_entities[rank].size();
    m_entities[rank].insert(m_entities[rank].end(), m_create_cache[rank].begin(), m_create_cache[rank].end());
    if (numOld > 0) {
      const EntityKey& firstNewKey = m_create_cache[rank][0].first;
      const EntityKey& lastOldKey = m_entities[rank][numOld-1].first;
      if (firstNewKey < lastOldKey) {
        EntityKeyEntityVector::iterator oldEnd = m_entities[rank].begin()+numOld;
        EntityKeyEntityVector::iterator loc = std::lower_bound(m_entities[rank].begin(), oldEnd, firstNewKey, EntityKeyEntityLess());
        std::inplace_merge(loc, oldEnd, m_entities[rank].end());
      }
    }
    m_create_cache[rank].clear();
    size_t num = std::max(m_entities[stk::topology::NODE_RANK].size(),
                          m_entities[stk::topology::ELEM_RANK].size());
    unsigned possibleCacheSize = num/1000;
    m_maxCreateCacheSize = std::max(m_maxCreateCacheSize, possibleCacheSize);
  }
}

void EntityRepository::clear_cache(EntityRank rank) const
{
  clear_created_entity_cache(rank);

  clear_updated_entity_cache(rank);

  clear_destroyed_entity_cache(rank);
}

std::pair<stk::mesh::entity_iterator,bool>
EntityRepository::add_to_cache(const EntityKey& key)
{
    bool inserted_new_entity = false;
    EntityRank rank = key.rank();
    EntityKeyEntityVector& cache = m_create_cache[rank];

    if (cache.size() >= m_maxCreateCacheSize) {
        clear_cache(rank);
    }

    EntityKeyEntityVector& entities = m_entities[rank];
    EntityKeyEntityVector::iterator iter = std::lower_bound(entities.begin(), entities.end(),
                                                            key, EntityKeyEntityLess());
    Entity entity;
    if (iter == entities.end() || iter->first != key) {
        cache.emplace_back(key, Entity());
        iter = cache.begin()+(cache.size()-1);
        inserted_new_entity = true;
    }
    else {
        inserted_new_entity = false;
    }
 
    if (cache.size() >= m_maxCreateCacheSize) {
        clear_cache(rank);
        iter = std::lower_bound(entities.begin(), entities.end(), key, EntityKeyEntityLess());
    }

    return std::make_pair(iter, inserted_new_entity);
}

stk::mesh::entity_iterator EntityRepository::get_from_cache(const EntityKey& key) const
{
  if (!m_create_cache[key.rank()].empty()) {
    EntityKeyEntityVector& cache = m_create_cache[key.rank()];
    EntityKeyEntityVector::iterator iter =
         std::find_if(cache.begin(), cache.end(), match_EntityKey(key));
    if (iter != cache.end()) {
      return iter;
    }
  }
  return m_create_cache[key.rank()].end();
}

std::pair<stk::mesh::entity_iterator ,bool>
EntityRepository::internal_create_entity( const EntityKey & key)
{
  if (key.rank() > m_entities.size()) {
    m_entities.resize(key.rank());
    m_create_cache.resize(key.rank());
    m_update_cache.resize(key.rank());
    m_destroy_cache.resize(key.rank());
  }

  clear_updated_entity_cache(key.rank());
  clear_destroyed_entity_cache(key.rank());

  entity_iterator ent = get_from_cache(key);
  if (ent != m_create_cache[key.rank()].end()) {
    return std::make_pair(ent, false);
  }

  return add_to_cache(key);
}

Entity EntityRepository::get_entity(const EntityKey &key) const
{
  EntityRank rank = key.rank();
  if (!m_destroy_cache[rank].empty()) {
    const std::vector<EntityKey>& destroyed = m_destroy_cache[rank];
    if (destroyed.size() < 64) {
      std::vector<EntityKey>::const_iterator iter = std::find(destroyed.begin(), destroyed.end(), key);
      if (iter != destroyed.end()) {
        return Entity();
      }
    }
    else {
      clear_destroyed_entity_cache(rank);
    }
  }

  clear_updated_entity_cache(rank);

  ThrowErrorMsgIf( ! key.is_valid(),
      "Invalid key: " << key.rank() << " " << key.id());

  if (rank >= m_entities.size()) {
    return Entity();
  }

  entity_iterator ent = get_from_cache(key);
  if (ent != m_create_cache[rank].end()) {
    return ent->second;
  }

  const EntityKeyEntityVector& entities = m_entities[rank];

  const EntityKeyEntityVector::const_iterator iter = std::lower_bound(entities.begin(), entities.end(), key, EntityKeyEntityLess());

  return (iter != entities.end() && (iter->first==key)) ? iter->second : Entity() ;
}

void EntityRepository::update_entity_key(EntityKey new_key, EntityKey old_key, Entity entity)
{
  EntityRank rank = new_key.rank();
  clear_created_entity_cache(rank);
  clear_destroyed_entity_cache(rank);

  if (m_update_cache[rank].size() >= m_maxUpdateCacheSize) {
    clear_cache(rank);
  }

  m_update_cache[rank].emplace_back(old_key, new_key);
}

void EntityRepository::destroy_entity(EntityKey key, Entity entity)
{ 
  EntityRank rank = key.rank();
  clear_created_entity_cache(rank);
  clear_updated_entity_cache(rank);

  if (m_destroy_cache[rank].size() >= m_maxUpdateCacheSize) {
    clear_cache(rank);
  }

  m_destroy_cache[rank].push_back(key);
} 

} // namespace impl
} // namespace mesh
} // namespace stk

