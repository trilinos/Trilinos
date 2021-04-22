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

#ifndef stk_mesh_EntityCommDatabase_hpp
#define stk_mesh_EntityCommDatabase_hpp

//----------------------------------------------------------------------

#include <algorithm>                    // for max
#include <functional>                   // for equal_to
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Types.hpp>      // for PairIterEntityComm, etc
#include <utility>                      // for pair
#include <vector>                       // for vector
#include <unordered_map>
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, hash_value
#include "stk_mesh/base/Ghosting.hpp"
#include "stk_util/util/NamedPair.hpp"
namespace stk { class CommBuffer; }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Relation; } }
namespace stk { namespace mesh { struct Entity; } }



//----------------------------------------------------------------------

namespace stk {
namespace mesh {

class CommMapChangeListener {
public:
    virtual ~CommMapChangeListener(){}
    virtual void removedKey(const EntityKey& key) = 0;
};

// Struct containing things the system must know about an entity to
// handle communication.
struct EntityComm
{
  bool isShared;
  bool isGhost;
  EntityCommInfoVector comm_map;
};

class EntityCommDatabase
{
  typedef std::pair<EntityKey const, EntityComm> map_value;
  typedef std::equal_to<EntityKey> map_predicate;
  typedef std::unordered_map<  EntityKey
                               , EntityComm
                               , stk::mesh::HashValueForEntityKey
                               , map_predicate
                              > map_type;

public:
  EntityCommDatabase() : m_comm_map(), m_last_lookup(m_comm_map.end()), m_comm_map_change_listener(nullptr) {}

  PairIterEntityComm shared_comm_info( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key, const Ghosting & sub ) const;

  std::pair<EntityComm*,bool> insert( const EntityKey & key, const EntityCommInfo & val, int owner );
  bool erase( const EntityKey & key, const EntityCommInfo & val );
  bool erase( const EntityKey & key, const Ghosting & ghost );
  bool comm_clear_ghosting(const EntityKey & key );
  bool comm_clear(const EntityKey & key );

  const EntityComm* entity_comm(const EntityKey& key) const;
        EntityComm* entity_comm(const EntityKey& key);

  void setCommMapChangeListener(CommMapChangeListener* listener)
  { m_comm_map_change_listener = listener; }

private:
  bool cached_find(const EntityKey& key) const;
  const EntityComm* insert(const EntityKey& key);
  void internal_update_shared_ghosted(bool removedSharingProc);

  map_type m_comm_map;
  mutable map_type::iterator m_last_lookup;
  mutable CommMapChangeListener* m_comm_map_change_listener;
};

//----------------------------------------------------------------------
inline
PairIterEntityComm shared_comm_info_range(const EntityCommInfoVector& comm_info_vector) {
    EntityCommInfoVector::const_iterator i = comm_info_vector.begin();
    EntityCommInfoVector::const_iterator end = comm_info_vector.end();
    EntityCommInfoVector::const_iterator e = i;

    while(e != end && e->ghost_id < 1) {
      ++e;
    }

    return PairIterEntityComm( i , e );
}

inline
PairIterEntityComm ghost_info_range(const EntityCommInfoVector& commInfo, const Ghosting & ghosting)
{
  EntityCommInfoVector::const_iterator ghostBegin = commInfo.begin();
  EntityCommInfoVector::const_iterator ghostEnd, end = commInfo.end();
  while(ghostBegin != end && ghostBegin->ghost_id != ghosting.ordinal()) {
    ++ghostBegin;
  } 
  
  if (ghostBegin != end) {
    ghostEnd = ghostBegin+1;
    while(ghostEnd != end && ghostEnd->ghost_id == ghosting.ordinal()) {
      ++ghostEnd;
    } 
    return PairIterEntityComm( ghostBegin , ghostEnd );
  } 

  return PairIterEntityComm( ghostBegin , end );
}

void pack_entity_info(const BulkData& mesh, CommBuffer & buf , const Entity entity );

void unpack_entity_info(
  CommBuffer     & buf,
  const BulkData & mesh ,
  EntityKey      & key ,
  int            & owner ,
  PartVector     & parts ,
  std::vector<Relation> & relations );

void pack_sideset_info(BulkData& mesh, CommBuffer & buf, const Entity entity);

void unpack_sideset_info(CommBuffer & buf, BulkData & mesh, const Entity entity);

/** \brief  Pack an entity's field values into a buffer */
void pack_field_values(const BulkData& mesh, CommBuffer & , Entity );

/** \brief  Unpack an entity's field values from a buffer */
bool unpack_field_values(const BulkData& mesh, CommBuffer & , Entity , std::ostream & error_msg );

} // namespace mesh
} // namespace stk

#endif
