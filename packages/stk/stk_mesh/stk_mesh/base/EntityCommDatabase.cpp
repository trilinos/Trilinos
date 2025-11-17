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

#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/Types.hpp"      // for ConnectivityOrdinal, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg
#include <algorithm>                       // for operator<<
#include <sstream>                      // for operator<<, basic_ostream
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Relation.hpp>   // for Relation
#include <string>                       // for operator<<


namespace stk {
namespace mesh {

EntityCommDatabase::EntityCommDatabase()
 : m_comm_map(),
   m_last_lookup(m_comm_map.end()),
   m_comm_map_change_listener(nullptr),
   m_entityCommInfo(0, EntityCommInfo(InvalidOrdinal, -1)),
   m_removedEntityCommIndices()
{
}

PairIterEntityComm EntityCommDatabase::shared_comm_info( const EntityKey & key ) const
{
  if (!cached_find(key)) return PairIterEntityComm();

  unsigned entityCommIndex = m_last_lookup->second;
  return shared_comm_info_range(m_entityCommInfo.items(entityCommIndex));
}

PairIterEntityComm EntityCommDatabase::comm( const EntityKey & key ) const
{
  if (!cached_find(key)) return PairIterEntityComm();

  return m_entityCommInfo.items(m_last_lookup->second);
}

PairIterEntityComm EntityCommDatabase::comm( const EntityKey & key, const Ghosting & sub ) const
{
  if (!cached_find(key)) return PairIterEntityComm();

  PairIterEntityComm comm_map = m_entityCommInfo.items(m_last_lookup->second);

  const EntityCommInfo s_begin( sub.ordinal() ,     0 );
  const EntityCommInfo s_end(   sub.ordinal() + 1 , 0 );

  const EntityCommInfo* i = comm_map.begin();
  const EntityCommInfo* e = comm_map.end();

  i = std::lower_bound( i , e , s_begin );
  e = std::lower_bound( i , e , s_end );

  return PairIterEntityComm( i , e );
}

PairIterEntityComm EntityCommDatabase::comm(unsigned entityCommIndex) const
{
  return m_entityCommInfo.items(entityCommIndex);
}

int EntityCommDatabase::entity_comm( const EntityKey & key ) const
{
  if (!cached_find(key)) return -1;

  return m_last_lookup->second;
}

// A cached find function that stores the result in m_last_lookup if successful and returns true.
// Otherwise, the find failed and it returns false.
bool EntityCommDatabase::cached_find(const EntityKey& key) const
{
  if (m_last_lookup != m_comm_map.end() && key == m_last_lookup->first) {
    return true;
  }

  map_type::iterator find_it = const_cast<map_type&>(m_comm_map).find(key);
  if (find_it == m_comm_map.end()) {
    return false;
  }
  else {
    m_last_lookup = find_it;
    return true;
  }
}

unsigned EntityCommDatabase::get_new_entity_comm_index()
{
  if (!m_removedEntityCommIndices.empty()) {
    unsigned newEntityCommIndex = m_removedEntityCommIndices.back();
    m_removedEntityCommIndices.pop_back();
    return newEntityCommIndex;
  }

  unsigned newEntityCommIndex = m_entityCommInfo.num_rows();
  m_entityCommInfo.add_row();
  return newEntityCommIndex;
}

int EntityCommDatabase::insert(const EntityKey& key)
{
  if (!cached_find(key)) {
    int newEntityCommIndex = static_cast<int>(get_new_entity_comm_index());
    std::pair<EntityKey,int> keyAndEntityComm(key, newEntityCommIndex);
    m_last_lookup = m_comm_map.insert(keyAndEntityComm).first;
  }
  return m_last_lookup->second;
}

std::pair<int,bool> EntityCommDatabase::insert( const EntityKey & key, const EntityCommInfo & val, int /*owner*/ )
{
  insert(key);

  int entityCommIndex = m_last_lookup->second;
  const bool didInsert = m_entityCommInfo.add_item(entityCommIndex, val);
  std::pair<int,bool> result = std::make_pair(m_last_lookup->second, didInsert);

  return result;
}

bool EntityCommDatabase::erase( const EntityKey & key, const EntityCommInfo & val )
{
  if (!cached_find(key)) return false;

  int entityCommIndex = m_last_lookup->second;
  const bool result = m_entityCommInfo.remove_item(entityCommIndex, val);

  if ( result ) {
    if (m_comm_map_change_listener != nullptr) {
      m_comm_map_change_listener->removedGhost(key, val.ghost_id, val.proc);
    }

    if (comm(entityCommIndex).empty()) {
      m_last_lookup = m_comm_map.erase(m_last_lookup);
      m_removedEntityCommIndices.push_back(entityCommIndex);

      if (m_comm_map_change_listener != nullptr) {
          m_comm_map_change_listener->removedKey(key);
      }
    }
  }

  return result ;
}

bool EntityCommDatabase::erase(unsigned entityCommIndex, const EntityKey& key, unsigned ghostID)
{
  const bool deletingSymmInfo = ghostID == BulkData::SYMM_INFO;

  bool result = m_entityCommInfo.remove_items_if(entityCommIndex, [&](const EntityCommInfo& info) {
    const bool shouldRemove = (info.ghost_id == ghostID) ||
                              (deletingSymmInfo && info.ghost_id >= BulkData::SYMM_INFO) ||
                              (info.ghost_id == BulkData::SYMM_INFO+ghostID);
    if (shouldRemove) {
      if (m_comm_map_change_listener != nullptr) {
        m_comm_map_change_listener->removedGhost(key, info.ghost_id, info.proc);
      }
      return true;
    }
    return false;
  });

  if ( result ) {
    if (comm(entityCommIndex).empty()) {
      cached_find(key);
      m_last_lookup = m_comm_map.erase(m_last_lookup);
      m_removedEntityCommIndices.push_back(entityCommIndex);

      if (m_comm_map_change_listener != nullptr) {
        m_comm_map_change_listener->removedKey(key);
      }
    }
  }

  return result ;
}

bool EntityCommDatabase::erase( const EntityKey & key, const Ghosting & ghost )
{
  return erase(key, ghost.ordinal());
}

bool EntityCommDatabase::erase( const EntityKey & key, unsigned ghostID )
{
  if (!cached_find(key)) return false;

  int entityCommIndex = m_last_lookup->second;

  return erase(entityCommIndex, key, ghostID);
}


bool EntityCommDatabase::comm_clear_ghosting(const EntityKey & key)
{
  if (!cached_find(key)) return false;

  int entityCommIndex = m_last_lookup->second;
  bool did_clear_ghosting = m_entityCommInfo.remove_items_if(entityCommIndex, [&](const EntityCommInfo& info) {
    if (info.ghost_id >= 1) {
      if (m_comm_map_change_listener != nullptr) {
        m_comm_map_change_listener->removedGhost(key, info.ghost_id, info.proc);
      }
      return true;
    }
    return false;
  });

  if (comm(entityCommIndex).empty()) {
    m_last_lookup = m_comm_map.erase(m_last_lookup);
    m_removedEntityCommIndices.push_back(entityCommIndex);
    if (m_comm_map_change_listener != nullptr) {
        m_comm_map_change_listener->removedKey(key);
    }
  }
  return did_clear_ghosting;
}

bool EntityCommDatabase::comm_clear(const EntityKey & key)
{
  if (!cached_find(key)) return false;

  int entityCommIndex = m_last_lookup->second;
  m_entityCommInfo.remove_items(entityCommIndex);
  m_last_lookup = m_comm_map.erase(m_last_lookup);
  m_removedEntityCommIndices.push_back(entityCommIndex);
  bool did_clear = true;
  if (m_comm_map_change_listener != nullptr) {
    m_comm_map_change_listener->removedKey(key);
  }
  return did_clear;
}

PairIterEntityCommListInfo EntityCommDatabase::comm_list_for_rank(EntityRank rank) const
{
  auto start = std::lower_bound(m_entity_comm_list.begin(), m_entity_comm_list.end(), rank, EntityCommListInfoRankLess());
  return PairIterEntityCommListInfo(&(*start),
                                    &(*std::upper_bound(start, m_entity_comm_list.end(), rank, EntityCommListInfoRankLess())));
}

} // namespace mesh
}

