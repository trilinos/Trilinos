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
#include <stk_mesh/baseImpl/EntityImpl.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {
namespace impl {

boost::fast_pool_allocator<EntityImpl>& entity_allocator()
{
  static boost::fast_pool_allocator<EntityImpl> entity_pool_allocator;
  return entity_pool_allocator;
}

Entity EntityRepository::allocate_entity(bool use_pool)
{
  if (use_pool) {
    static EntityImpl tmp_entity;
    EntityImpl* new_entity = entity_allocator().allocate();
    entity_allocator().construct(new_entity, tmp_entity);
    Entity rv = {new_entity};
    return rv;
  }
  else {
    Entity rv = {new EntityImpl};
    return rv;
  }
}

void EntityRepository::internal_destroy_entity(Entity entity)
{
  if (m_use_pool) {
    entity_allocator().destroy(entity.m_entityImpl);
    entity_allocator().deallocate(entity.m_entityImpl, 1);
  }
  else {
    delete entity.m_entityImpl;
  }
}

EntityRepository::~EntityRepository()
{
  if (m_use_pool) {
    boost::singleton_pool<boost::fast_pool_allocator_tag, sizeof(EntityImpl)>::release_memory();
  }
  else {
    while ( !m_entities.empty() ) {
      internal_expunge_entity( m_entities.begin() );
    }
  }
}

void EntityRepository::internal_expunge_entity( EntityMap::iterator i )
{
  ThrowAssertMsg( i->second != Entity(),
                  "For key " << entity_rank(i->first) << " " <<
                  entity_id(i->first) << ", value was NULL");

  ThrowAssertMsg( i->first == i->second.key(),
                  "Key (" << i->first.rank() << ", " << i->first.id() << ")" <<
                  " != (" << i->second.key().rank() << ", " << i->second.key().rank() << ")");


  Entity deleted_entity = i->second;
  internal_destroy_entity(deleted_entity);
  m_entities.erase( i );
}

Entity
EntityRepository::internal_allocate_entity(EntityKey entity_key)
{
  Entity new_entity = allocate_entity(m_use_pool);
  new_entity.set_key(entity_key);
  return new_entity;
}

std::pair<Entity ,bool>
EntityRepository::internal_create_entity( const EntityKey & key )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::internal_create_entity", LOG_ENTITY, key);

  EntityMap::value_type tmp(key, Entity());

  const std::pair< EntityMap::iterator , bool >
    insert_result = m_entities.insert( tmp );

  std::pair<Entity ,bool>
    result( insert_result.first->second , insert_result.second );

  if ( insert_result.second )  { // A new entity
    Entity new_entity = internal_allocate_entity(key);
    insert_result.first->second = result.first = new_entity;
  }

  return result ;
}

void EntityRepository::log_created_parallel_copy( Entity entity )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::log_created_parallel_copy", LOG_ENTITY, entity.key());

  entity.m_entityImpl->created_parallel_copy();
}

Entity EntityRepository::get_entity(const EntityKey &key) const
{
  ThrowErrorMsgIf( ! entity_key_valid( key ),
      "Invalid key: " << entity_rank(key) << " " << entity_id(key));

  const EntityMap::const_iterator i = m_entities.find( key );

  return i != m_entities.end() ? i->second : Entity() ;
}

void EntityRepository::clean_changes()
{
  TraceIf("stk::mesh::impl::EntityRepository::clean_changes", LOG_ENTITY);

  for ( EntityMap::iterator i = m_entities.begin() ; i != m_entities.end() ; ++i ) {
    i->second.m_entityImpl->clear_state();
  }
}

void EntityRepository::change_entity_bucket( Bucket & b, Entity e,
                                             unsigned ordinal)
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::change_entity_bucket", LOG_ENTITY, e.key());
  DiagIfWatching(LOG_ENTITY, e.key(), "New bucket: " << b << ", ordinal: " << ordinal);

  const bool modified_parts = ! e.m_entityImpl->is_bucket_valid() ||
                              ! b.in_same_partition( e.bucket() );
  if ( modified_parts ) {
    e.m_entityImpl->modified();
  }
  e.m_entityImpl->set_bucket_and_ordinal( &b, ordinal);
}

Bucket * EntityRepository::get_entity_bucket( Entity e ) const
{
  // Note, this allows for returning NULL bucket
  return e.m_entityImpl->bucket_ptr();
}

bool EntityRepository::destroy_relation( Entity e_from,
                                         Entity e_to,
                                         const RelationIdentifier local_id )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::destroy_relation", LOG_ENTITY, e_from.key());

  bool caused_change_fwd = e_from.m_entityImpl->destroy_relation(e_to, local_id);

  // Relationships should always be symmetrical
  if ( caused_change_fwd ) {
    bool caused_change_inv = e_to.m_entityImpl->destroy_relation(e_from, local_id);
    ThrowErrorMsgIf( !caused_change_inv,
        " Internal error - could not destroy inverse relation of " <<
        print_entity_key( e_from ) << " to " << print_entity_key( e_to ) <<
        " with local relation id of " << local_id);
  }

  // It is critical that the modification be done AFTER the relations are
  // changed so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    e_to.m_entityImpl->modified();
    e_from.m_entityImpl->modified();
  }

  return caused_change_fwd;
}

void EntityRepository::declare_relation( Entity e_from,
                                         Entity e_to,
                                         const RelationIdentifier local_id,
                                         unsigned sync_count )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::declare_relation", LOG_ENTITY, e_from.key());

  bool caused_change_fwd =
    e_from.m_entityImpl->declare_relation( e_to, local_id, sync_count);

  // Relationships should always be symmetrical
  if ( caused_change_fwd ) {

    // the setup for the converse relationship works slightly differently
    bool is_converse = true;
    bool caused_change_inv =
      e_to.m_entityImpl->declare_relation( e_from, local_id, sync_count,
                                          is_converse );

    ThrowErrorMsgIf( !caused_change_inv,
        " Internal error - could not create inverse relation of " <<
        print_entity_key( e_from ) << " to " << print_entity_key( e_to ));
  }

  // It is critical that the modification be done AFTER the relations are
  // added so that the propagation can happen correctly.
  if ( caused_change_fwd ) {
    e_to.m_entityImpl->modified();
    e_from.m_entityImpl->modified();
  }
}

void EntityRepository::update_entity_key(EntityKey new_key, Entity entity)
{
  EntityKey old_key = entity.key();

  ThrowAssert(m_entities.find(new_key) == m_entities.end());
  m_entities.insert(std::make_pair(new_key, entity));
  entity.m_entityImpl->update_key(new_key);

  EntityMap::iterator old_itr = m_entities.find( old_key );
  ThrowAssert(old_itr != m_entities.end());
  m_entities.erase(old_itr);
}

} // namespace impl
} // namespace mesh
} // namespace stk
