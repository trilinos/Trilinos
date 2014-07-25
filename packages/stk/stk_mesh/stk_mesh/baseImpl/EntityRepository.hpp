/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_baseImpl_EntityRepository_hpp
#define stk_mesh_baseImpl_EntityRepository_hpp

#include <stddef.h>                     // for size_t
#include <map>                          // for map, map<>::value_compare
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <utility>                      // for pair
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Types.hpp"      // for EntityRank
namespace stk { namespace mesh { class Bucket; } }
namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace mesh {
namespace impl {

class EntityRepository {

public:

  typedef std::map<EntityKey,Entity> EntityMap;

    typedef EntityMap::const_iterator const_iterator;
    typedef EntityMap::iterator iterator;

    EntityRepository(BulkData &mesh)
      : m_mesh(mesh), m_entities() {}

    ~EntityRepository();

    Entity get_entity( const EntityKey &key ) const;

    const_iterator begin() const { return m_entities.begin(); }
    const_iterator end() const { return m_entities.end(); }

    const_iterator begin_rank(EntityRank ent_rank) const { return m_entities.lower_bound(EntityKey(ent_rank, 0)); }
    const_iterator end_rank(EntityRank ent_rank) const { return m_entities.upper_bound(EntityKey(static_cast<EntityRank>(ent_rank+1), 0)); }

    // Return a pair: the relevant entity, and whether it had to be created
    // or not. If there was already an active entity with the specified key, the second item in the pair
    // will be false; otherwise it will be true (even if the Entity was present
    // but marked as destroyed).
    std::pair<Entity ,bool>
    internal_create_entity( const EntityKey & key, size_t preferred_offset = 0 );

    void change_entity_bucket( Bucket & b, Entity e, unsigned ordinal);

    void update_entity_key(EntityKey new_key, EntityKey old_key, Entity entity);

    void destroy_entity(EntityKey key, Entity entity);

  private:
    void internal_expunge_entity( EntityMap::iterator i);

    // Entity internal_allocate_entity(EntityKey entity_key);
    Entity allocate_entity();

    BulkData &m_mesh;
    EntityMap m_entities;

    //disable copy constructor and assignment operator
    EntityRepository(const EntityRepository &);
    EntityRepository & operator =(const EntityRepository &);
};

inline
void EntityRepository::destroy_entity(EntityKey key, Entity entity)
{
  m_entities.erase(key);
}

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_EntityRepository_hpp
