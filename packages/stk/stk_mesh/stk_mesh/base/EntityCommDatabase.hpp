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
#include "stk_mesh/base/HashEntityAndEntityKey.hpp"
#include "stk_util/util/MCSR.hpp"
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
    virtual void removedGhost(const EntityKey& key, unsigned ghostId, int proc) = 0;
    virtual void removedKey(const EntityKey& key) = 0;
};

class EntityCommDatabase
{
  typedef std::unordered_map<EntityKey, int, std::hash<EntityKey>> map_type;

public:
  EntityCommDatabase();

  PairIterEntityComm shared_comm_info( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key, const Ghosting & sub ) const;
  PairIterEntityComm comm( unsigned entityCommIndex) const;

  std::pair<int,bool> insert( const EntityKey & key, const EntityCommInfo & val, int owner );
  bool erase( const EntityKey & key, const EntityCommInfo & val );
  bool erase( const EntityKey & key, const Ghosting & ghost );
  bool comm_clear_ghosting(const EntityKey & key );
  bool comm_clear(const EntityKey & key );

  int entity_comm(const EntityKey& key) const;

  void setCommMapChangeListener(CommMapChangeListener* listener)
  { m_comm_map_change_listener = listener; }

  size_t num_comm_keys() const { return m_comm_map.size(); }
  size_t bucket_count() const { return m_comm_map.bucket_count(); }

private:
  bool cached_find(const EntityKey& key) const;
  unsigned get_new_entity_comm_index();
  int insert(const EntityKey& key);
  void internal_update_shared_ghosted(bool removedSharingProc);

  map_type m_comm_map;
  mutable map_type::iterator m_last_lookup;
  mutable CommMapChangeListener* m_comm_map_change_listener;
  stk::util::MCSR<EntityCommInfo> m_entityCommInfo;
  std::vector<int> m_removedEntityCommIndices;
};

//----------------------------------------------------------------------
inline
PairIterEntityComm shared_comm_info_range(PairIterEntityComm commInfo) {
    const EntityCommInfo* i = commInfo.begin();
    const EntityCommInfo* end = commInfo.end();
    const EntityCommInfo* e = i;

    while(e != end && e->ghost_id < 1) {
      ++e;
    }

    return PairIterEntityComm( i , e );
}

inline
PairIterEntityComm ghost_info_range(PairIterEntityComm commInfo, unsigned ghostingOrdinal)
{
  const EntityCommInfo* ghostBegin = commInfo.begin();
  const EntityCommInfo* ghostEnd, *end = commInfo.end();
  while(ghostBegin != end && ghostBegin->ghost_id != ghostingOrdinal) {
    ++ghostBegin;
  } 
  
  if (ghostBegin != end) {
    ghostEnd = ghostBegin+1;
    while(ghostEnd != end && ghostEnd->ghost_id == ghostingOrdinal) {
      ++ghostEnd;
    } 
    return PairIterEntityComm( ghostBegin , ghostEnd );
  } 

  return PairIterEntityComm( ghostBegin , end );
}

void pack_entity_info(const BulkData& mesh,
                      CommBuffer& buf,
                      const Entity entity,
                      bool onlyPackDownwardRelations = false);

void unpack_entity_info(
  CommBuffer     & buf,
  const BulkData & mesh ,
  EntityKey      & key ,
  int            & owner ,
  PartVector     & parts ,
  RelationVector& relations );

void pack_sideset_info(BulkData& mesh, CommBuffer & buf, const Entity entity);

void unpack_sideset_info(CommBuffer & buf, BulkData & mesh, const Entity entity);

/** \brief  Pack an entity's field values into a buffer */
void pack_field_values(const BulkData& mesh, CommBuffer & , Entity );

/** \brief  Unpack an entity's field values from a buffer */
bool unpack_field_values(const BulkData& mesh, CommBuffer & , Entity , std::ostream & error_msg );

} // namespace mesh
} // namespace stk

#endif
