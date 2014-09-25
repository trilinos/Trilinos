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

#ifndef stk_mesh_EntityCommDatabase_hpp
#define stk_mesh_EntityCommDatabase_hpp

//----------------------------------------------------------------------

#include <algorithm>                    // for max
#include <functional>                   // for equal_to
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Types.hpp>      // for PairIterEntityComm, etc
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "boost/functional/hash/extensions.hpp"  // for hash
#include "boost/unordered/unordered_map.hpp"  // for unordered_map
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, hash_value
#include "stk_util/util/NamedPair.hpp"
#include "stk_util/util/TrackingAllocator.hpp"  // for tracking_allocator
namespace stk { class CommBuffer; }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Ghosting; } }
namespace stk { namespace mesh { class Relation; } }
namespace stk { namespace mesh { struct Entity; } }



//----------------------------------------------------------------------

namespace stk {
namespace mesh {

// Struct containing things the system must know about an entity to
// handle communication. The Entity object itself won't necessarily be
// available, so there's some duplication here (owner_rank).
struct EntityComm
{
  EntityCommInfoVector comm_map;
  int                  owner_rank;
};

class EntityCommDatabase
{
  typedef std::pair<EntityKey const, EntityComm> map_value;
  typedef tracking_allocator< map_value, EntityCommTag> map_allocator;
  typedef boost::hash<EntityKey> map_hash;
  typedef std::equal_to<EntityKey> map_predicate;
  typedef boost::unordered_map<  EntityKey
                               , EntityComm
                               , map_hash
                               , map_predicate
                               , map_allocator
                              > map_type;

public:
  EntityCommDatabase() : m_comm_map(), m_last_lookup(m_comm_map.end()) {}

  PairIterEntityComm shared_comm_info( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key, const Ghosting & sub ) const;
  int owner_rank( const EntityKey & key ) const;

  bool insert( const EntityKey & key, const EntityCommInfo & val, int owner );
  bool erase( const EntityKey & key, const EntityCommInfo & val );
  bool erase( const EntityKey & key, const Ghosting & ghost );
  void comm_clear_ghosting(const EntityKey & key );
  void comm_clear(const EntityKey & key );
  bool change_owner_rank(const EntityKey& key, int owner);

  const EntityComm* entity_comm(const EntityKey& key) const;

private:
  bool cached_find(const EntityKey& key) const;
  void insert(const EntityKey& key);

  map_type m_comm_map;
  mutable map_type::iterator m_last_lookup;
};

//----------------------------------------------------------------------

void pack_entity_info(const BulkData& mesh, CommBuffer & buf , const Entity entity );

void unpack_entity_info(
  CommBuffer     & buf,
  const BulkData & mesh ,
  EntityKey      & key ,
  int            & owner ,
  PartVector     & parts ,
  std::vector<Relation> & relations );

/** \brief  Pack an entity's field values into a buffer */
void pack_field_values(const BulkData& mesh, CommBuffer & , Entity );

/** \brief  Unpack an entity's field values from a buffer */
bool unpack_field_values(const BulkData& mesh, CommBuffer & , Entity , std::ostream & error_msg );

} // namespace mesh
} // namespace stk

#endif
