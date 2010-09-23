/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_base_Entity_hpp
#define stk_mesh_base_Entity_hpp

#include <utility>
#include <vector>

#include <stk_mesh/base/Types.hpp>

#include <stk_mesh/baseImpl/EntityImpl.hpp>

namespace stk {
namespace mesh {

namespace impl {

class EntityRepository;
class BucketRepository;

}

/** \addtogroup stk_mesh_module
 * \{
 */

//----------------------------------------------------------------------



//----------------------------------------------------------------------
/** \brief  A fundamental unit within the discretization of a problem domain,
 *          including but not limited to nodes, edges, sides, and elements.
 *
 *  Entities are distributed among parallel processors.
 *  A given entity may reside on more than one processor;
 *  however, it is owned by exactly one of the processors
 *  on which it resides.
 *
 *  Note that an Entity's state comprises:
 *   - existence - Whether this entity has been created and is not destroyed
 *   - owner - The rank of the owning process
 *   - part-membership - The set of parts this entity belongs to
 *   - relations - Relationships between other Entities
 *  When any of the above changes, the Entity's log state may change
 */
class Entity {
public:

  /** \brief  Query the current state of the entity log */
  EntityModificationLog log_query() const { return m_entityImpl.log_query(); }

  /** \brief  The rank of this entity. */
  EntityRank entity_rank() const { return m_entityImpl.entity_rank(); }

  /** \brief  Identifier for this entity which is globally unique
   *          for a given entity type.
   */
  EntityId identifier() const { return m_entityImpl.identifier(); }

  /** \brief  The globally unique key ( entity type + identifier )
   *          of this entity.
   */
  const EntityKey & key() const { return m_entityImpl.key(); }

  /** \brief  The bucket which holds this mesh entity's field data */
  Bucket & bucket() const { return m_entityImpl.bucket(); }

  /** \brief  The ordinal for this entity within its bucket. */
  unsigned bucket_ordinal() const { return m_entityImpl.bucket_ordinal(); }

  /** \brief  The mesh bulk data synchronized_count when this entity's
   *          part membership was most recently modified.
   *
   *  If ( mesh.synchronized_state() == false &&
   *       mesh.synchronized_count() == entity.synchronized_count() )
   *  then entity was modified during this modification phase.
   */
  size_t synchronized_count() const { return m_entityImpl.synchronized_count(); }

  //------------------------------------
  /** \brief  All \ref stk::mesh::Relation "Entity relations"
   *          for which this entity is a member. The relations are ordered
   *          from lowest entity-rank to highest entity-rank.
   */
  PairIterRelation relations() const { return m_entityImpl.relations(); }

  /** \brief  \ref stk::mesh::Relation "Entity relations" for which this
   *          entity is a member, the other entity is of a given type.
   */
  PairIterRelation relations( EntityRank type ) const { return m_entityImpl.relations(type); }

  //------------------------------------
  /** \brief  Parallel processor rank of the processor which owns this entity */
  unsigned owner_rank() const { return m_entityImpl.owner_rank(); }

  /** \brief  Parallel processes which share this entity. */
  PairIterEntityComm sharing() const { return m_entityImpl.sharing(); }

  /** \brief  Complete communicaiton list for this entity */
  PairIterEntityComm comm() const { return m_entityImpl.comm(); }

  /** \brief  Subset communicaiton list for this entity */
  PairIterEntityComm comm( const Ghosting & sub ) const { return m_entityImpl.comm( sub ); }

  //------------------------------------
private:

  impl::EntityImpl m_entityImpl;

  ~Entity();
  explicit Entity( const EntityKey & arg_key );

  Entity(); ///< Default constructor not allowed
  Entity( const Entity & ); ///< Copy constructor not allowed
  Entity & operator = ( const Entity & ); ///< Assignment operator not allowed

#ifndef DOXYGEN_COMPILE
  friend class impl::EntityRepository ;
  friend class impl::EntityImpl ;

#endif /* DOXYGEN_COMPILE */
};

/** \brief  Comparison operator for entities compares the entities' keys */
class EntityLess {
public:
  ~EntityLess() {}
  EntityLess() {}
  EntityLess( const EntityLess & ) {}
  EntityLess & operator = ( const EntityLess & ) { return *this ; }

  /** \brief  Comparison operator */
  bool operator()(const Entity& lhs, const Entity& rhs) const
  { return lhs.key() < rhs.key(); }

  bool operator()(const Entity& lhs, const EntityKey & rhs) const
  { return lhs.key() < rhs ; }

  /** \brief  Comparison operator */
  bool operator()(const Entity* lhs, const Entity* rhs) const
  {
    const EntityKey lhs_key = lhs ? lhs->key() : EntityKey() ;
    const EntityKey rhs_key = rhs ? rhs->key() : EntityKey() ;
    return lhs_key < rhs_key ;
  }

  bool operator()(const Entity* lhs, const Entity& rhs) const
  {
    const EntityKey lhs_key = lhs ? lhs->key() : EntityKey();
    return lhs_key < rhs.key() ;
  }

  bool operator()(const Entity& lhs, const Entity* rhs) const
  {
    const EntityKey rhs_key = rhs ? rhs->key() : EntityKey();
    return lhs.key() < rhs_key ;
  }

  bool operator()(const Entity* lhs, const EntityKey & rhs) const
  {
    const EntityKey lhs_key = lhs ? lhs->key() : EntityKey() ;
    return lhs_key < rhs ;
  }

  bool operator()( const EntityProc & lhs, const EntityProc & rhs) const
  {
    const EntityKey lhs_key = lhs.first ? lhs.first->key() : EntityKey() ;
    const EntityKey rhs_key = rhs.first ? rhs.first->key() : EntityKey() ;
    return lhs_key != rhs_key ? lhs_key < rhs_key : lhs.second < rhs.second ;
  }

  bool operator()( const EntityProc & lhs, const Entity & rhs) const
  {
    const EntityKey lhs_key = lhs.first ? lhs.first->key() : EntityKey() ;
    return lhs_key < rhs.key();
  }

  bool operator()( const EntityProc & lhs, const Entity * rhs) const
  {
    const EntityKey lhs_key = lhs.first ? lhs.first->key() : EntityKey() ;
    const EntityKey rhs_key = rhs       ? rhs->key() : EntityKey() ;
    return lhs_key < rhs_key ;
  }

  bool operator()( const EntityProc & lhs, const EntityKey & rhs) const
  {
    const EntityKey lhs_key = lhs.first ? lhs.first->key() : EntityKey() ;
    return lhs_key < rhs ;
  }

}; //class EntityLess

class EntityEqual
{
public:
  bool operator()(const stk::mesh::Entity* lhs, const stk::mesh::Entity* rhs) const
  {
    const stk::mesh::EntityKey lhs_key = lhs ? lhs->key() : stk::mesh::EntityKey();
    const stk::mesh::EntityKey rhs_key = rhs ? rhs->key() : stk::mesh::EntityKey();
    return lhs_key == rhs_key;
  }

  bool operator()(const stk::mesh::Entity& lhs, const stk::mesh::Entity& rhs) const
  {
    const stk::mesh::EntityKey lhs_key = lhs.key();
    const stk::mesh::EntityKey rhs_key = rhs.key();
    return lhs_key == rhs_key;
  }
};

//----------------------------------------------------------------------

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_base_Entity_hpp */

