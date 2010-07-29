/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_baseImpl_EntityRepository_hpp
#define stk_mesh_baseImpl_EntityRepository_hpp

#include <map>
#include <stk_mesh/base/Entity.hpp>

namespace stk {
namespace mesh {
namespace impl {

class EntityRepository {
    typedef std::map<EntityKey,Entity*> EntityMap;
  public:
    typedef std::map<EntityKey,Entity*>::const_iterator iterator;

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

    void destroy_relation( Entity & e_from, Entity & e_to);

    void declare_relation( Entity & e_from,
                           Entity & e_to,
                           const unsigned local_id,
                           unsigned sync_count );

  private:
    void internal_expunge_entity( EntityMap::iterator i);

    EntityMap                           m_entities ;

    //disabel copy constructor and assignment operator
    EntityRepository(const EntityRepository &);
    EntityRepository & operator =(const EntityRepository &);
};

/*---------------------------------------------------------------*/

void EntityRepository::set_entity_sync_count( Entity & e, size_t count) {
  e.m_entityImpl.set_sync_count(count);
}

void EntityRepository::set_entity_owner_rank( Entity & e, unsigned owner_rank) {
  e.m_entityImpl.set_owner_rank(owner_rank);
}

void EntityRepository::comm_clear( Entity & e) const {
  e.m_entityImpl.comm_clear();
}

void EntityRepository::comm_clear_ghosting( Entity & e) const {
  e.m_entityImpl.comm_clear_ghosting();
}

} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_EntityRepository_hpp
