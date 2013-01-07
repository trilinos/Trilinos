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

#include <boost/pool/pool_alloc.hpp>

namespace stk {
namespace mesh {
namespace impl {

class EntityRepository {

#if STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1
  struct stk_entity_rep_hash : public std::unary_function< EntityKey, std::size_t >
  {
    std::size_t operator()(const EntityKey& x) const
    {
      return (std::size_t)(x.raw_key());
    }
  };

  typedef std::tr1::unordered_map<EntityKey, Entity , stk_entity_rep_hash, std::equal_to<EntityKey> > EntityMap;
#else
  typedef std::map<EntityKey,Entity> EntityMap;
#endif

  public:

    typedef EntityMap::const_iterator iterator;

    EntityRepository(bool use_pool)
      : m_entities(), m_use_pool(use_pool) {}

    ~EntityRepository();

    Entity get_entity( const EntityKey &key ) const;

    iterator begin() const { return m_entities.begin(); }
    iterator end() const { return m_entities.end(); }

    void clean_changes();

    // Return a pair: the relevant entity, and whether it had to be created
    // or not. If there was already an active entity, the second item in the
    // will be false; otherwise it will be true (even if the Entity was present
    // but marked as destroyed).
    std::pair<Entity ,bool>
      internal_create_entity( const EntityKey & key );

    /** \brief Log that this entity was created as a parallel copy
      *        of an existing entity.
      */
    void log_created_parallel_copy( Entity e );

    /**
     * The client knows that this entity should be marked as modified. In
     * general clients shouldn't need to call this because EntityRepository
     * knows when it performs operations that modify entities. BulkData should
     * be the only caller of this method.
     */
    void log_modified(Entity e) const;

    bool set_entity_owner_rank( Entity e, unsigned owner_rank);
    void set_entity_sync_count( Entity e, size_t count);

    void change_entity_bucket( Bucket & b, Entity e, unsigned ordinal);
    Bucket * get_entity_bucket ( Entity e ) const;

    bool destroy_relation( Entity e_from,
                           Entity e_to,
                           const RelationIdentifier local_id);

    void declare_relation( Entity e_from,
                           Entity e_to,
                           const RelationIdentifier local_id,
                           unsigned sync_count );

    void update_entity_key(EntityKey new_key, EntityKey old_key);

    void destroy_entity(Entity entity);

  private:
    void internal_expunge_entity( EntityMap::iterator i);
    void internal_destroy_entity(Entity entity);

    Entity internal_allocate_entity(EntityKey entity_key);
    Entity allocate_entity(bool use_pool);

    EntityMap m_entities;
    bool m_use_pool;

    //disable copy constructor and assignment operator
    EntityRepository(const EntityRepository &);
    EntityRepository & operator =(const EntityRepository &);
};

/*---------------------------------------------------------------*/

inline
void EntityRepository::set_entity_sync_count( Entity e, size_t count)
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::set_entity_sync_count", LOG_ENTITY, e.key());

  e.m_entityImpl->set_sync_count(count);
}

inline
bool EntityRepository::set_entity_owner_rank( Entity e, unsigned owner_rank)
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::set_entity_owner_rank", LOG_ENTITY, e.key());
  DiagIfWatching(LOG_ENTITY, e.key(), "new owner: " << owner_rank);

  bool changed = e.m_entityImpl->set_owner_rank(owner_rank);
  if ( changed ) {
    e.m_entityImpl->modified();
  }
  return changed;
}

inline
void EntityRepository::log_modified( Entity e ) const
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::log_modified", LOG_ENTITY, e.key());

  e.m_entityImpl->modified();
}

inline
void EntityRepository::destroy_entity(Entity entity)
{
  m_entities.erase(entity.key());
  internal_destroy_entity(entity);
}

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_EntityRepository_hpp
