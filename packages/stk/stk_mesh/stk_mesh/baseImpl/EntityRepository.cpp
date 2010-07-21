/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <stdexcept>

#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {
namespace impl {

namespace {

const MetaData & metadata_from_entity(const Entity *e) {
  return e->bucket().mesh().mesh_meta_data();
}


}

EntityRepository::~EntityRepository() {
  try {
    while ( ! m_entities.empty() ) {
      internal_expunge_entity( m_entities.begin() );
    }
  } catch(...){}
}


void EntityRepository::internal_expunge_entity( EntityMap::iterator i )
{
  const bool ok_ptr = i->second != NULL ;
  const bool ok_key = ok_ptr ? i->first == i->second->key() : true ;

  if ( ! ok_ptr || ! ok_key ) {
    std::ostringstream msg ;
    msg << "stk::mesh::baseImple::EntityRepository::internal_expunge_entity( " ;
    print_entity_key( msg , metadata_from_entity(i->second) , i->first );
    if ( ! ok_ptr ) {
      msg << "NULL" ;
    }
    else {
      msg << " != " ;
      print_entity_key( msg , metadata_from_entity(i->second) , i->second->key() );
    }
    msg << ") FAILED" ;
    throw std::runtime_error( msg.str() );
  }

  delete i->second ;
  i->second = NULL ;
  m_entities.erase( i );
}

std::pair<Entity*,bool>
EntityRepository::internal_create_entity( const EntityKey & key )
{
  EntityMap::value_type tmp(key,NULL);

  const std::pair< EntityMap::iterator , bool >
    insert_result = m_entities.insert( tmp );

  std::pair<Entity*,bool>
    result( insert_result.first->second , insert_result.second );

  if ( insert_result.second )  { // A new entity
    insert_result.first->second = result.first = new Entity( key );
  }
  else if ( result.first->marked_for_destruction() ) {
    // resurrection
    result.first->m_entityImpl.log_resurrect();
    result.second = true;
  }

  return result ;
}

Entity * EntityRepository::get_entity(const EntityKey &key) const
{
  const bool valid_key = entity_key_valid( key );

  const EntityMap::const_iterator i = m_entities.find( key );

  if ( ! valid_key ) {
    static const char method[] = "stk::mesh::BulkData::get_entity" ;
    std::ostringstream msg ;
    msg << method << "( " ;
    // QUESTION: If the key is invalid, then will 'i' ever be a valid iterator?
    if (i != m_entities.end()) {
      print_entity_key( msg , metadata_from_entity(i->second) , key );
    } else {
      const unsigned type = entity_rank(key);
      const EntityId id   = entity_id(key);
      msg << "Invalid key: " << type << " " << id;
    }
    msg << " INVALID KEY" ;
    msg << " ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }

  return i != m_entities.end() ? i->second : NULL ;
}

void EntityRepository::clean_changes() {

  for ( EntityMap::iterator
      i = m_entities.begin() ; i != m_entities.end() ; )
  {
    const EntityMap::iterator j = i ;
    ++i ;

    if ( j->second->m_entityImpl.marked_for_destruction() ) {
      // Clear out the entities destroyed in the previous modification.
      // They were retained for change-logging purposes.
      internal_expunge_entity( j );
    }
    else {
      j->second->m_entityImpl.log_clear();
    }
  }
}

bool EntityRepository::erase_ghosting( Entity & e, const Ghosting & ghosts) const {
  return e.m_entityImpl.erase( ghosts );
}


bool EntityRepository::erase_comm_info( Entity & e, const EntityCommInfo & comm_info) const {
  return e.m_entityImpl.erase( comm_info );
}

bool EntityRepository::insert_comm_info( Entity & e, const EntityCommInfo & comm_info) const {
  return e.m_entityImpl.insert( comm_info );
}

void EntityRepository::destroy_later( Entity & e, Bucket* nil_bucket ) {
  change_entity_bucket( *nil_bucket, e, 0);
  e.m_entityImpl.log_deleted(); //important that this come last
}

void EntityRepository::change_entity_bucket( Bucket & b, Entity & e,
                                             unsigned ordinal) {
  const bool modified_parts = ! e.m_entityImpl.is_bucket_valid() ||
                              ! b.equivalent( e.bucket() );
  if ( modified_parts ) {
    e.m_entityImpl.log_modified();
  }
  e.m_entityImpl.set_bucket_and_ordinal( &b, ordinal);
}

} // namespace impl
} // namespace mesh
} // namespace stk

