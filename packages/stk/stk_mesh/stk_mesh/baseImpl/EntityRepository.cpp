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

EntityRepository::~EntityRepository()
{
  try {
    while ( ! m_entities.empty() ) {
      internal_expunge_entity( m_entities.begin() );
    }
  } catch(...){}
}

void EntityRepository::internal_expunge_entity( EntityMap::iterator i )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::internal_expunge_entity", LOG_ENTITY, i->first);

  ThrowErrorMsgIf( i->second == NULL,
                   "For key " << entity_rank(i->first) << " " <<
                   entity_id(i->first) << ", value was NULL");

  ThrowErrorMsgIf( i->first != i->second->key(),
    "Key " << print_entity_key(MetaData::get( *i->second ), i->first) <<
    " != " << print_entity_key(i->second));

  delete i->second ;
  i->second = NULL ;
  m_entities.erase( i );
}

std::pair<Entity*,bool>
EntityRepository::internal_create_entity( const EntityKey & key )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::internal_create_entity", LOG_ENTITY, key);

  EntityMap::value_type tmp(key,NULL);

  const std::pair< EntityMap::iterator , bool >
    insert_result = m_entities.insert( tmp );

  std::pair<Entity*,bool>
    result( insert_result.first->second , insert_result.second );

  if ( insert_result.second )  { // A new entity
    insert_result.first->second = result.first = new Entity( key );
  }
  else if ( EntityLogDeleted == result.first->log_query() ) {
    // resurrection
    result.first->m_entityImpl.log_resurrect();
    result.second = true;
  }

  return result ;
}

void EntityRepository::log_created_parallel_copy( Entity & entity )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::log_created_parallel_copy", LOG_ENTITY, entity.key());

  entity.m_entityImpl.log_created_parallel_copy();
}

Entity * EntityRepository::get_entity(const EntityKey &key) const
{
  ThrowErrorMsgIf( ! entity_key_valid( key ),
      "Invalid key: " << entity_rank(key) << " " << entity_id(key));

  const EntityMap::const_iterator i = m_entities.find( key );

  return i != m_entities.end() ? i->second : NULL ;
}

void EntityRepository::clean_changes()
{
  TraceIf("stk::mesh::impl::EntityRepository::clean_changes", LOG_ENTITY);

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

bool EntityRepository::erase_ghosting( Entity & e, const Ghosting & ghosts) const
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::erase_ghosting", LOG_ENTITY, e.key());

  return e.m_entityImpl.erase( ghosts );
}

bool EntityRepository::erase_comm_info( Entity & e, const EntityCommInfo & comm_info) const
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::erase_comm_info", LOG_ENTITY, e.key());

  return e.m_entityImpl.erase( comm_info );
}

bool EntityRepository::insert_comm_info( Entity & e, const EntityCommInfo & comm_info) const
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::insert_comm_info", LOG_ENTITY, e.key());

  return e.m_entityImpl.insert( comm_info );
}

void EntityRepository::destroy_later( Entity & e, Bucket* nil_bucket )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::destroy_later", LOG_ENTITY, e.key());

  ThrowErrorMsgIf( e.log_query() == EntityLogDeleted,
                   "double deletion of entity: " << print_entity_key( e ));

  change_entity_bucket( *nil_bucket, e, 0);
  e.m_entityImpl.log_deleted(); //important that this come last
}

void EntityRepository::change_entity_bucket( Bucket & b, Entity & e,
                                             unsigned ordinal)
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::change_entity_bucket", LOG_ENTITY, e.key());
  DiagIfWatching(LOG_ENTITY, e.key(), "New bucket: " << b << ", ordinal: " << ordinal);

  const bool modified_parts = ! e.m_entityImpl.is_bucket_valid() ||
                              ! b.equivalent( e.bucket() );
  if ( modified_parts ) {
    e.m_entityImpl.log_modified_and_propagate();
  }
  e.m_entityImpl.set_bucket_and_ordinal( &b, ordinal);
}

Bucket * EntityRepository::get_entity_bucket( Entity & e ) const
{
  // Note, this allows for returning NULL bucket
  return e.m_entityImpl.bucket_ptr();
}

bool EntityRepository::destroy_relation( Entity & e_from,
                                         Entity & e_to,
                                         const RelationIdentifier local_id )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::destroy_relation", LOG_ENTITY, e_from.key());

  bool caused_change_fwd = e_from.m_entityImpl.destroy_relation(e_to, local_id);

  // Relationships should always be symmetrical
  if ( caused_change_fwd ) {
    bool caused_change_inv = e_to.m_entityImpl.destroy_relation(e_from, local_id);
    ThrowErrorMsgIf( !caused_change_inv,
        " Internal error - could not destroy inverse relation of " <<
        print_entity_key( e_from ) << " to " << print_entity_key( e_to ) <<
        " with local relation id of " << local_id);
  }

  // It is critical that the modification be done AFTER the relations are
  // changed so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    e_to.m_entityImpl.log_modified_and_propagate();
    e_from.m_entityImpl.log_modified_and_propagate();
  }

  return caused_change_fwd;
}

void EntityRepository::declare_relation( Entity & e_from,
                                         Entity & e_to,
                                         const RelationIdentifier local_id,
                                         unsigned sync_count )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::declare_relation", LOG_ENTITY, e_from.key());

  bool caused_change_fwd =
    e_from.m_entityImpl.declare_relation( e_to, local_id, sync_count);

  // Relationships should always be symmetrical
  if ( caused_change_fwd ) {

    // the setup for the converse relationship works slightly differently
    bool is_converse = true;
    bool caused_change_inv =
      e_to.m_entityImpl.declare_relation( e_from, local_id, sync_count,
                                          is_converse );

    ThrowErrorMsgIf( !caused_change_inv,
        " Internal error - could not create inverse relation of " <<
        print_entity_key( e_from ) << " to " << print_entity_key( e_to ));
  }

  // It is critical that the modification be done AFTER the relations are
  // added so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    e_to.m_entityImpl.log_modified_and_propagate();
    e_from.m_entityImpl.log_modified_and_propagate();
  }
}

} // namespace impl
} // namespace mesh
} // namespace stk
