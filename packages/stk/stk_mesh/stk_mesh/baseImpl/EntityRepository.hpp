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

    Entity * get_entity( const EntityKey key ) const ;

    iterator begin() const { return m_entities.begin(); }
    iterator end() const { return m_entities.end(); }

    void clean_changes();

    std::pair<Entity*,bool>
      internal_create_entity( const EntityKey & key );

  private:
    void internal_expunge_entity( EntityMap::iterator i);

    EntityMap                           m_entities ;

    //disabel copy constructor and assignment operator
    EntityRepository(const EntityRepository &);
    EntityRepository & operator =(const EntityRepository &);
};


} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_EntityRepository_hpp
