/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_baseImpl_EntityRepository_hpp
#define stk_mesh_baseImpl_EntityRepository_hpp

#include <stk_mesh/base/Trace.hpp>

// We will use tr1 if we can (not on PGI or pathscale); otherwise, fall back to std map.
#if defined(__PGI) || defined(__PATHSCALE__)
  #define STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1 0
#else
  #define STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1 0
#endif

#if STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1
  #include <tr1/unordered_map>
#else
  #include <map>
#endif

#include <stk_mesh/base/Entity.hpp>

namespace stk {
namespace mesh {
namespace impl {

class EntityRepository {

public:

#if STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1
  struct stk_entity_rep_hash : public std::unary_function< EntityKey, std::size_t >
  {
    std::size_t operator()(const EntityKey& x) const
    {
      return (std::size_t)(x);
    }
  };

  typedef std::tr1::unordered_map<EntityKey, Entity , stk_entity_rep_hash, std::equal_to<EntityKey> > EntityMap;
#else
  typedef std::map<EntityKey,Entity> EntityMap;
#endif

    typedef EntityMap::const_iterator const_iterator;
    typedef EntityMap::iterator iterator;

    EntityRepository(BulkData &mesh)
      : m_mesh(mesh), m_entities() {}

    ~EntityRepository();

    Entity get_entity( const EntityKey &key ) const;

    const_iterator begin() const { return m_entities.begin(); }
    const_iterator end() const { return m_entities.end(); }

#if STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1
    //The following begin_rank/end_rank methods use std::map::lower_bound, which won't work
    //if we're using an unordered map.
#else
    const_iterator begin_rank(EntityRank ent_rank) const { return m_entities.lower_bound(EntityKey(ent_rank, 0)); }
    const_iterator end_rank(EntityRank ent_rank) const { return m_entities.upper_bound(EntityKey(ent_rank+1, 0)); }
#endif

    // Return a pair: the relevant entity, and whether it had to be created
    // or not. If there was already an active entity, the second item in the
    // will be false; otherwise it will be true (even if the Entity was present
    // but marked as destroyed).
    std::pair<Entity ,bool>
      internal_create_entity( const EntityKey & key );

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
