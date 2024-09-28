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

#include <stk_mesh/baseImpl/EntityKeyMapping.hpp>
#include <stddef.h>                     // for NULL
#include <sstream>                      // for operator<<, basic_ostream, etc
#include <vector>
#include "stk_mesh/base/Entity.hpp"     // for Entity, etc
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssert, etc


namespace stk {
namespace mesh {
namespace impl {

struct EntityKeyEntityLess {
inline bool operator()(const std::pair<EntityKey,Entity>& lhsPair, const std::pair<EntityKey,Entity>& rhsPair) const
{
  return lhsPair.first.m_value < rhsPair.first.m_value;
}
inline bool operator()(const std::pair<EntityKey,Entity>& key_ent_pair, const EntityKey& key) const
{
  return key_ent_pair.first.m_value < key.m_value;
}
inline bool operator()(const EntityKey& key, const std::pair<EntityKey,Entity>& key_ent_pair) const
{
  return key.m_value < key_ent_pair.first.m_value;
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

EntityKeyMapping::EntityKeyMapping()
 : m_entities(stk::topology::NUM_RANKS),
   m_create_cache(stk::topology::NUM_RANKS),
   m_update_cache(stk::topology::NUM_RANKS),
   m_destroy_cache(stk::topology::NUM_RANKS),
   m_maxCreateCacheSize(512),
   m_maxUpdateCacheSize(4096)
{
}

EntityKeyMapping::~EntityKeyMapping()
{
}

void EntityKeyMapping::clear_all_cache()
{
  EntityRank nRanks = static_cast<EntityRank>(m_create_cache.size());
  for(EntityRank rank=stk::topology::BEGIN_RANK; rank<nRanks; ++rank) {
    clear_cache(rank);
  }
}

void EntityKeyMapping::clear_destroyed_entity_cache(EntityRank rank) const
{
  if (!m_destroy_cache[rank].empty()) {
    std::vector<EntityKey>& destroy = m_destroy_cache[rank];
    std::sort(destroy.begin(), destroy.end());
    EntityKeyEntityVector& entities = m_entities[rank];
    size_t destroyIdx = 0;
    EntityKeyEntityVector::iterator start = std::lower_bound(entities.begin(), entities.end(), destroy[0], EntityKeyEntityLess());
    EntityKeyEntityVector::iterator end = std::upper_bound(entities.begin(), entities.end(), destroy.back(), EntityKeyEntityLess());
    size_t startIdx = std::distance(entities.begin(), start);
    size_t endIdx = std::distance(entities.begin(), end);
    size_t keep = startIdx;
    for(size_t i=startIdx; i<endIdx; ++i) {
      if (destroyIdx < destroy.size() && entities[i].first == destroy[destroyIdx]) {
        ++destroyIdx;
        continue;
      }
      if (i > keep) {
        entities[keep] = entities[i];
      }
      ++keep;
    }
    if (endIdx < entities.size()) {
      size_t len = entities.size() - endIdx;
      for(size_t i=0; i<len; ++i) {
        entities[keep+i] = entities[endIdx+i];
      }
      keep += len;
    }
    entities.resize(keep);

    destroy.clear();
    size_t num = std::max(m_entities[stk::topology::NODE_RANK].size(),
                          m_entities[stk::topology::ELEM_RANK].size());
    unsigned possibleCacheSize = num/1000;
    m_maxUpdateCacheSize = std::max(m_maxUpdateCacheSize, possibleCacheSize);
  }
}

void EntityKeyMapping::clear_updated_entity_cache(EntityRank rank) const
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

void EntityKeyMapping::clear_created_entity_cache(EntityRank rank) const
{
  if (!m_create_cache[rank].empty()) {
    std::sort(m_create_cache[rank].begin(), m_create_cache[rank].end());
    stk::util::insert_keep_sorted(m_create_cache[rank], m_entities[rank], EntityKeyEntityLess());
    m_create_cache[rank].clear();
    size_t num = std::max(m_entities[stk::topology::NODE_RANK].size(),
                          m_entities[stk::topology::ELEM_RANK].size());
    unsigned possibleCacheSize = num/1000;
    m_maxCreateCacheSize = std::max(m_maxCreateCacheSize, possibleCacheSize);
  }
}

void EntityKeyMapping::clear_cache(EntityRank rank) const
{
  clear_created_entity_cache(rank);

  clear_updated_entity_cache(rank);

  clear_destroyed_entity_cache(rank);
}

std::pair<stk::mesh::entity_iterator,bool>
EntityKeyMapping::add_to_cache(const EntityKey& key)
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

stk::mesh::entity_iterator EntityKeyMapping::get_from_cache(const EntityKey& key) const
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
EntityKeyMapping::internal_create_entity( const EntityKey & key)
{
  if (key.rank() > entity_rank_count()) {
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

Entity EntityKeyMapping::get_entity(const EntityKey &key) const
{
  if (!key.is_valid()) {
    return Entity();
  }

  EntityRank rank = key.rank();
  if (rank >= entity_rank_count()) {
    return Entity();
  }

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

  entity_iterator ent = get_from_cache(key);
  if (ent != m_create_cache[rank].end()) {
    return ent->second;
  }

  const EntityKeyEntityVector& entities = m_entities[rank];

  const EntityKeyEntityVector::const_iterator iter = std::lower_bound(entities.begin(), entities.end(), key, EntityKeyEntityLess());

  return (iter != entities.end() && (iter->first==key)) ? iter->second : Entity() ;
}

void EntityKeyMapping::update_entity_key(EntityKey new_key, EntityKey old_key, Entity entity)
{
  EntityRank rank = new_key.rank();
  clear_created_entity_cache(rank);
  clear_destroyed_entity_cache(rank);

  if (m_update_cache[rank].size() >= m_maxUpdateCacheSize) {
    clear_cache(rank);
  }

  m_update_cache[rank].emplace_back(old_key, new_key);
}

void EntityKeyMapping::destroy_entity(EntityKey key, Entity entity)
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

