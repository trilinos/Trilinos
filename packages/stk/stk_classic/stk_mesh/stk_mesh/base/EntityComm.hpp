/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_EntityComm_hpp
#define stk_mesh_EntityComm_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <vector>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Ghosting.hpp>

#include <boost/unordered_map.hpp>

//----------------------------------------------------------------------

namespace stk_classic {
namespace mesh {

class EntityComm
{
public:
  typedef boost::unordered_map<EntityKey, EntityCommInfoVector> map_type;

  PairIterEntityComm sharing( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key ) const;
  PairIterEntityComm comm( const EntityKey & key, const Ghosting & sub ) const;

  bool insert( const EntityKey & key, const EntityCommInfo & val );
  bool erase( const EntityKey & key, const EntityCommInfo & val );
  bool erase( const EntityKey & key, const Ghosting & ghost );
  void comm_clear_ghosting(const EntityKey & key );
  void comm_clear(const EntityKey & key );
  void comm_swap(const EntityKey & key1, const EntityKey & key2);

private:
  map_type m_comm_map;
};

inline PairIterEntityComm EntityComm::sharing( const EntityKey & key ) const
{
  map_type::const_iterator it = m_comm_map.find(key);
  if (it == m_comm_map.cend()) {
    return PairIterEntityComm();
  }
  const EntityCommInfoVector & m_comm = it->second;

  EntityCommInfoVector::const_iterator i = m_comm.begin();
  EntityCommInfoVector::const_iterator e = m_comm.end();

  e = std::lower_bound( i , e , EntityCommInfo(1,     // ghost id, 1->aura
                                               0 ) ); // proc

  // Contains everything up the first aura comm (IE, only contains shared comms)
  return PairIterEntityComm( i , e );
}

inline PairIterEntityComm EntityComm::comm( const EntityKey & key ) const
{
  map_type::const_iterator it = m_comm_map.find(key);
  if (it == m_comm_map.cend()) {
    return PairIterEntityComm();
  }
  const EntityCommInfoVector & m_comm = it->second;
  return PairIterEntityComm(m_comm);
}

inline PairIterEntityComm EntityComm::comm( const EntityKey & key, const Ghosting & sub ) const
{
  map_type::const_iterator it = m_comm_map.find(key);
  if (it == m_comm_map.cend()) {
    return PairIterEntityComm();
  }
  const EntityCommInfoVector & m_comm = it->second;

  const EntityCommInfo s_begin( sub.ordinal() ,     0 );
  const EntityCommInfo s_end(   sub.ordinal() + 1 , 0 );

  EntityCommInfoVector::const_iterator i = m_comm.begin();
  EntityCommInfoVector::const_iterator e = m_comm.end();

  i = std::lower_bound( i , e , s_begin );
  e = std::lower_bound( i , e , s_end );

  return PairIterEntityComm( i , e );
}

inline bool EntityComm::insert( const EntityKey & key, const EntityCommInfo & val )
{
  TraceIfWatching("stk_classic::mesh::EntityComm::insert", LOG_ENTITY, key());
  EntityCommInfoVector & m_comm = m_comm_map[key];

  std::vector< EntityCommInfo >::iterator i =
    std::lower_bound( m_comm.begin() , m_comm.end() , val );

  const bool result = ((i == m_comm.end()) || (val != *i));

  if ( result ) {
    m_comm.insert( i , val );
  }

  return result ;
}

inline bool EntityComm::erase( const EntityKey & key, const EntityCommInfo & val )
{
  TraceIfWatching("stk_classic::mesh::EntityComm::erase(comm)", LOG_ENTITY, key());
  EntityCommInfoVector & m_comm = m_comm_map[key];

  std::vector< EntityCommInfo >::iterator i =
    std::lower_bound( m_comm.begin() , m_comm.end() , val );

  const bool result = ( (i != m_comm.end()) && (val == *i) ) ;

  if ( result ) {
    m_comm.erase( i );
  }

  return result ;
}

inline bool EntityComm::erase( const EntityKey & key, const Ghosting & ghost )
{
  TraceIfWatching("stk_classic::mesh::EntityComm::erase(ghost)", LOG_ENTITY, key());
  EntityCommInfoVector & m_comm = m_comm_map[key];

  const EntityCommInfo s_begin( ghost.ordinal() ,     0 );
  const EntityCommInfo s_end(   ghost.ordinal() + 1 , 0 );

  EntityCommInfoVector::iterator i = m_comm.begin();
  EntityCommInfoVector::iterator e = m_comm.end();

  i = std::lower_bound( i , e , s_begin );
  e = std::lower_bound( i , e , s_end );

  const bool result = i != e ;

  if ( result ) {
    m_comm.erase( i , e );
  }

  return result ;
}

inline void EntityComm::comm_clear_ghosting(const EntityKey & key)
{
  TraceIfWatching("stk_classic::mesh::EntityComm::comm_clear_ghosting", LOG_ENTITY, key());
  EntityCommInfoVector & m_comm = m_comm_map[key];

  std::vector< EntityCommInfo >::iterator j = m_comm.begin();
  while ( j != m_comm.end() && j->ghost_id == 0 ) { ++j ; }
  m_comm.erase( j , m_comm.end() );
}

inline void EntityComm::comm_clear(const EntityKey & key)
{
  TraceIfWatching("stk_classic::mesh::EntityComm::comm_clear", LOG_ENTITY, key());
  EntityCommInfoVector& commvec = m_comm_map[key];
  commvec.clear();
}

inline void EntityComm::comm_swap(const EntityKey & key1, const EntityKey & key2)
{
  map_type::iterator it1 = m_comm_map.find(key1);
  map_type::iterator it2 = m_comm_map.find(key2);

  if (it1 == m_comm_map.cend() && it2 == m_comm_map.cend()) {
    return;
  }

  EntityCommInfoVector & comm1 = it1 == m_comm_map.cend() ? m_comm_map[key1] : it1->second;
  EntityCommInfoVector & comm2 = it2 == m_comm_map.cend() ? m_comm_map[key2] : it2->second;

  comm1.swap(comm2);
}

/** \brief  Is shared with any other process */
bool in_shared( const Entity & entity );

/** \brief  Is shared with a given process */
bool in_shared( const Entity & entity , unsigned proc );

/** \brief  Is a receive ghost copy of an entity */
bool in_receive_ghost( const Entity & entity );

/** \brief  Is a receive ghost copy of an entity */
bool in_receive_ghost( const Ghosting & ghost , const Entity & entity );

/** \brief  Is sent to a ghost copy of an entity on any process */
bool in_send_ghost( const Entity & entity );

/** \brief  Is sent to a ghost copy of an entity on a given process */
bool in_send_ghost( const Entity & entity , unsigned proc );

/** \brief  Is in ghosting either send to 'p' or receive from 'p' */
bool in_ghost( const Ghosting & ghost , const Entity & entity , unsigned p );

/** \brief  Is in owned closure of the given process,
 *          typically the local process.
 */
bool in_owned_closure( const Entity & entity , unsigned proc );

/** \brief  List of all entity communication processes, sorted */
void comm_procs( const Entity & entity , std::vector<unsigned> & procs );

/** \brief  List of entity communication processes for a given ghost, sorted */
void comm_procs( const Ghosting & ghost ,
                 const Entity & entity , std::vector<unsigned> & procs );

//----------------------------------------------------------------------

void pack_entity_info( CommBuffer & buf , const Entity & entity );

void unpack_entity_info(
  CommBuffer     & buf,
  const BulkData & mesh ,
  EntityKey      & key ,
  unsigned       & owner ,
  PartVector     & parts ,
  std::vector<Relation> & relations );

/** \brief  Pack an entity's field values into a buffer */
void pack_field_values( CommBuffer & , Entity & );

/** \brief  Unpack an entity's field values from a buffer */
bool unpack_field_values( CommBuffer & , Entity & , std::ostream & error_msg );

} // namespace mesh
} // namespace stk_classic

#endif
