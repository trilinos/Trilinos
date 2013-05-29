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

Entity EntityRepository::allocate_entity(bool use_pool)
{
    Entity retval;
    retval = Entity::InvalidEntity;
    return retval;
}

EntityRepository::~EntityRepository()
{
  while ( !m_entities.empty() ) {
    internal_expunge_entity( m_entities.begin() );
  }
}

void EntityRepository::internal_expunge_entity( EntityMap::iterator i )
{
  ThrowAssertMsg( i->second != Entity(),
                  "For key " << i->first.rank() << " " <<
                  i->first.id() << ", value was NULL");

  m_entities.erase( i );
}

//Entity
//EntityRepository::internal_allocate_entity(EntityKey entity_key)
//{
//  Entity new_entity = allocate_entity(m_use_pool);
//  // This makes no sense (not yet added to m_mesh!)and breaks things!
//  m_mesh.set_entity_key(new_entity, entity_key);
//  return new_entity;
//}

std::pair<Entity ,bool>
EntityRepository::internal_create_entity( const EntityKey & key )
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::internal_create_entity", LOG_ENTITY, key);

  bool inserted_new_entity = false;
  EntityMap::iterator iter = m_entities.lower_bound(key);

  if (iter == m_entities.end() || iter->first != key) {
    Entity next_entity = {Entity::InvalidEntity};
    next_entity.set_bulk_data_id(m_mesh.bulk_data_id());
    next_entity.set_local_offset(m_mesh.generate_next_local_offset());
    m_mesh.set_entity_key(next_entity, key);

    iter = m_entities.insert(iter, std::make_pair(key, next_entity));
    inserted_new_entity = true;
  }

  return std::make_pair(iter->second, inserted_new_entity);
}

Entity EntityRepository::get_entity(const EntityKey &key) const
{
  ThrowErrorMsgIf( ! key.is_valid(),
      "Invalid key: " << key.rank() << " " << key.id());

  const EntityMap::const_iterator i = m_entities.find( key );

  return i != m_entities.end() ? i->second : Entity() ;
}

void EntityRepository::change_entity_bucket( Bucket & b, Entity e,
                                             unsigned ordinal)
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::change_entity_bucket", LOG_ENTITY, b.mesh().entity_key(e));
  DiagIfWatching(LOG_ENTITY, b.mesh().entity_key(e), "New bucket: " << b << ", ordinal: " << ordinal);

  const bool modified_parts = ! (b.mesh().bucket_ptr(e) != NULL) ||
                              ! b.in_same_partition( b.mesh().bucket(e) );
  if ( modified_parts ) {
    b.mesh().modified(e);
  }
  b.mesh().set_mesh_index(e, &b, ordinal);
}

void EntityRepository::update_entity_key(EntityKey new_key, EntityKey old_key, Entity entity)
{
  ThrowAssert(m_entities.find(new_key) == m_entities.end());
  m_entities.insert(std::make_pair(new_key, entity));

  EntityMap::iterator old_itr = m_entities.find( old_key );
  ThrowAssert(old_itr != m_entities.end());
  m_entities.erase(old_itr);
}

} // namespace impl
} // namespace mesh
} // namespace stk
