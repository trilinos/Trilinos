/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_EntityCommDatabase_hpp
#define stk_mesh_EntityCommDatabase_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <vector>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Ghosting.hpp>

#include <boost/unordered_map.hpp>

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
  typedef boost::unordered_map<EntityKey, EntityComm> map_type;

public:
  EntityCommDatabase() : m_comm_map(), m_last_lookup(m_comm_map.end()) {}

  PairIterEntityComm sharing( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key, const Ghosting & sub ) const;
  int owner_rank( const EntityKey & key ) const;

  bool insert( const EntityKey & key, const EntityCommInfo & val, int owner );
  bool erase( const EntityKey & key, const EntityCommInfo & val );
  bool erase( const EntityKey & key, const Ghosting & ghost );
  void comm_clear_ghosting(const EntityKey & key );
  void comm_clear(const EntityKey & key );
  bool change_owner_rank(const EntityKey& key, int owner);

private:
  bool cached_find(const EntityKey& key) const;
  void insert(const EntityKey& key);

  map_type m_comm_map;
  mutable map_type::iterator m_last_lookup;
};

inline
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

inline
void EntityCommDatabase::insert(const EntityKey& key)
{
  if (!cached_find(key)) {
    m_last_lookup = m_comm_map.insert(std::make_pair(key, EntityComm())).first;
  }
}

inline
int EntityCommDatabase::owner_rank( const EntityKey & key ) const
{
  if (!cached_find(key)) return InvalidProcessRank;

  return m_last_lookup->second.owner_rank;
}

inline
PairIterEntityComm EntityCommDatabase::sharing( const EntityKey & key ) const
{
  if (!cached_find(key)) return PairIterEntityComm();

  const EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;

  EntityCommInfoVector::const_iterator i = comm_map.begin();
  EntityCommInfoVector::const_iterator e = comm_map.end();

  e = std::lower_bound( i , e , EntityCommInfo(1,     // ghost id, 1->aura
                                               0 ) ); // proc

  // Contains everything up the first aura comm (IE, only contains shared comms)
  return PairIterEntityComm( i , e );
}

inline
PairIterEntityComm EntityCommDatabase::comm( const EntityKey & key ) const
{
  if (!cached_find(key)) return PairIterEntityComm();

  const EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;
  return PairIterEntityComm(comm_map);
}

inline
PairIterEntityComm EntityCommDatabase::comm( const EntityKey & key, const Ghosting & sub ) const
{
  if (!cached_find(key)) return PairIterEntityComm();

  const EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;

  const EntityCommInfo s_begin( sub.ordinal() ,     0 );
  const EntityCommInfo s_end(   sub.ordinal() + 1 , 0 );

  EntityCommInfoVector::const_iterator i = comm_map.begin();
  EntityCommInfoVector::const_iterator e = comm_map.end();

  i = std::lower_bound( i , e , s_begin );
  e = std::lower_bound( i , e , s_end );

  return PairIterEntityComm( i , e );
}

inline
bool EntityCommDatabase::insert( const EntityKey & key, const EntityCommInfo & val, int owner )
{
  TraceIfWatching("stk::mesh::EntityComm::insert", LOG_ENTITY, key);
  DiagIfWatching(LOG_ENTITY, key, "owner " << owner);

  insert(key);
  EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;
  m_last_lookup->second.owner_rank = owner;

  std::vector< EntityCommInfo >::iterator i =
    std::lower_bound( comm_map.begin() , comm_map.end() , val );

  const bool result = ((i == comm_map.end()) || (val != *i));

  if ( result ) {
    comm_map.insert( i , val );
  }

  return result;
}

inline
bool EntityCommDatabase::erase( const EntityKey & key, const EntityCommInfo & val )
{
  TraceIfWatching("stk::mesh::EntityComm::erase(comm)", LOG_ENTITY, key);

  if (!cached_find(key)) return false;

  EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;

  std::vector< EntityCommInfo >::iterator i =
    std::lower_bound( comm_map.begin() , comm_map.end() , val );

  const bool result = ( (i != comm_map.end()) && (val == *i) ) ;

  if ( result ) {
    comm_map.erase( i );
    if (comm_map.empty()) {
      m_last_lookup = m_comm_map.erase(m_last_lookup);
    }
  }

  return result ;
}

inline
bool EntityCommDatabase::erase( const EntityKey & key, const Ghosting & ghost )
{
  TraceIfWatching("stk::mesh::EntityComm::erase(ghost)", LOG_ENTITY, key);

  if (!cached_find(key)) return false;

  EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;

  const EntityCommInfo s_begin( ghost.ordinal() ,     0 );
  const EntityCommInfo s_end(   ghost.ordinal() + 1 , 0 );

  EntityCommInfoVector::iterator i = comm_map.begin();
  EntityCommInfoVector::iterator e = comm_map.end();

  i = std::lower_bound( i , e , s_begin );
  e = std::lower_bound( i , e , s_end );

  const bool result = i != e ;

  if ( result ) {
    comm_map.erase( i , e );
    if (comm_map.empty()) {
      m_last_lookup = m_comm_map.erase(m_last_lookup);
    }
  }

  // if there is no more comm info, just remove it from the map?

  return result ;
}

inline
void EntityCommDatabase::comm_clear_ghosting(const EntityKey & key)
{
  TraceIfWatching("stk::mesh::EntityComm::comm_clear_ghosting", LOG_ENTITY, key);

  if (!cached_find(key)) return;

  EntityCommInfoVector & comm_map = m_last_lookup->second.comm_map;

  std::vector< EntityCommInfo >::iterator j = comm_map.begin();
  while ( j != comm_map.end() && j->ghost_id == 0 ) { ++j ; }
  comm_map.erase( j , comm_map.end() );

  if (comm_map.empty()) {
    m_last_lookup = m_comm_map.erase(m_last_lookup);
  }
}

inline
void EntityCommDatabase::comm_clear(const EntityKey & key)
{
  TraceIfWatching("stk::mesh::EntityComm::comm_clear", LOG_ENTITY, key);

  if (!cached_find(key)) return;

  m_last_lookup = m_comm_map.erase(m_last_lookup);
}

inline
bool EntityCommDatabase::change_owner_rank(const EntityKey& key, int owner)
{
  // Do not add key to map, only update rank if it was already in the map
  if (cached_find(key)) {
    TraceIfWatching("stk::mesh::EntityComm::change_owner_rank", LOG_ENTITY, key);
    DiagIfWatching(LOG_ENTITY, key, "new owner " << owner);

    const int orig_owner = m_last_lookup->second.owner_rank;
    m_last_lookup->second.owner_rank = owner;
    return orig_owner != owner;
  }
  return false;
}


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
