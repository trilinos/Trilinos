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

#if STK_MESH_ENTITYREPOSITORY_MAP_TYPE_TR1
  struct stk_entity_rep_hash : public std::unary_function< EntityKey, std::size_t >
  {
    inline std::size_t
    operator()(const EntityKey& x) const
    {
      return (std::size_t)(x.raw_key());
    }
  };

  typedef std::tr1::unordered_map<EntityKey, Entity*, stk_entity_rep_hash, std::equal_to<EntityKey> > EntityMap;
#else
  typedef std::map<EntityKey,Entity*> EntityMap;
#endif

  public:

    typedef EntityMap::const_iterator iterator;

    EntityRepository() : m_entities() {}
    ~EntityRepository();

    Entity * get_entity( const EntityKey &key ) const;

    iterator begin() const { return m_entities.begin(); }
    iterator end() const { return m_entities.end(); }

    void clean_changes();

    // Return a pair: the relevant entity, and whether it had to be created
    // or not. If there was already an active entity, the second item in the
    // will be false; otherwise it will be true (even if the Entity was present
    // but marked as destroyed).
    std::pair<Entity*,bool>
      internal_create_entity( const EntityKey & key );

    /** \brief Log that this entity was created as a parallel copy
      *        of an existing entity.
      */
    void log_created_parallel_copy( Entity & e );

    /**
     * The client knows that this entity should be marked as modified. In
     * general clients shouldn't need to call this because EntityRepository
     * knows when it performs operations that modify entities. BulkData should
     * be the only caller of this method.
     */
    inline void log_modified(Entity & e) const;

    inline void set_entity_owner_rank( Entity & e, unsigned owner_rank);
    inline void set_entity_sync_count( Entity & e, size_t count);

    inline void comm_clear( Entity & e) const;
    inline void comm_clear_ghosting( Entity & e) const;

    bool erase_ghosting( Entity & e, const Ghosting & ghosts) const;
    bool erase_comm_info( Entity & e, const EntityCommInfo & comm_info) const;

    bool insert_comm_info( Entity & e, const EntityCommInfo & comm_info) const;

    void change_entity_bucket( Bucket & b, Entity & e, unsigned ordinal);
    Bucket * get_entity_bucket ( Entity & e ) const;
    void destroy_later( Entity & e, Bucket* nil_bucket );

    bool destroy_relation( Entity & e_from,
                           Entity & e_to,
                           const RelationIdentifier local_id);

    void declare_relation( Entity & e_from,
                           Entity & e_to,
                           const RelationIdentifier local_id,
                           unsigned sync_count );

  private:
    void internal_expunge_entity( EntityMap::iterator i);

    EntityMap m_entities;

    //disabel copy constructor and assignment operator
    EntityRepository(const EntityRepository &);
    EntityRepository & operator =(const EntityRepository &);
};

/*---------------------------------------------------------------*/

void EntityRepository::set_entity_sync_count( Entity & e, size_t count)
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::set_entity_sync_count", LOG_ENTITY, e.key());

  e.m_entityImpl.set_sync_count(count);
}

void EntityRepository::set_entity_owner_rank( Entity & e, unsigned owner_rank)
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::set_entity_owner_rank", LOG_ENTITY, e.key());
  DiagIfWatching(LOG_ENTITY, e.key(), "new owner: " << owner_rank);

  bool changed = e.m_entityImpl.set_owner_rank(owner_rank);
  if ( changed ) {
    e.m_entityImpl.log_modified_and_propagate();
  }
}

void EntityRepository::comm_clear( Entity & e) const
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::comm_clear", LOG_ENTITY, e.key());

  e.m_entityImpl.comm_clear();
}

void EntityRepository::comm_clear_ghosting( Entity & e) const
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::comm_clear_ghosting", LOG_ENTITY, e.key());

  e.m_entityImpl.comm_clear_ghosting();
}

void EntityRepository::log_modified( Entity & e ) const
{
  TraceIfWatching("stk::mesh::impl::EntityRepository::log_modified", LOG_ENTITY, e.key());

  e.m_entityImpl.log_modified_and_propagate();
}

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_EntityRepository_hpp
