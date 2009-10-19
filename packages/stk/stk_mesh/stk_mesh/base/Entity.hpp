#ifndef stk_mesh_Entity_hpp
#define stk_mesh_Entity_hpp

#include <limits>
#include <iosfwd>

#include <stk_util/util/PairIter.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Relation.hpp>

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 * \{
 */

//----------------------------------------------------------------------

/** \brief  Span of a sorted relations for a given domain entity.
 *
 *  The span is sorted by
 *  -# range entity type,
 *  -# relation type,
 *  -# relation identifier, and
 *  -# range entity global identifier.
 */
typedef PairIter< std::vector<Relation>::const_iterator > PairIterRelation ;

//----------------------------------------------------------------------
/** \brief  A fundamental unit within the discretization of a problem domain,
 *          including but not limited to nodes, edges, sides, and elements.
 *
 *  Entities are distributed among parallel processors.
 *  A given entity may reside on more than one processor;
 *  however, it is owned by exactly one of the processors
 *  on which it resides.
 */
class Entity {
public:

  /** \brief  The type (a.k.a. rank) of this entity. */
  EntityType entity_type() const { return stk::mesh::entity_type( m_key ); }

  /** \brief  Identifier for this entity which is globally unique
   *          for a given entity type.
   */
  EntityId identifier() const { return stk::mesh::entity_id( m_key ); }

  /** \brief  The globally unique key ( entity type + identifier )
   *          of this entity.
   */
  const EntityKey & key() const { return m_key ; }

  /** \brief  The bucket which holds this mesh entity's field data */
  const Bucket & bucket() const { return *m_bucket ; }

  /** \brief  The ordinal for this entity within its bucket. */
  unsigned bucket_ordinal() const { return m_bucket_ord ; }

  /** \brief  The mesh bulk data synchronized_count when this entity's
   *          part membership was most recently modified.
   *
   *  If ( mesh.synchronized_state() == false &&
   *       mesh.synchronized_count() == entity.synchronized_count() )
   *  then entity was modified during this modification phase.
   */
  size_t synchronized_count() const { return m_sync_count ; }

  //------------------------------------
  /** \brief  All \ref stk::mesh::Relation "Entity relations"
   *          for which this entity is a member.
   */
  PairIterRelation relations() const { return PairIterRelation( m_relation ); }

  /** \brief  \ref stk::mesh::Relation "Entity relations" for which this
   *          entity is a member, the other entity is of a given type,
   *          and the relation is of the given kind (topological = 0).
   */
  PairIterRelation relations( unsigned type , unsigned kind = 0 ) const ;

  //------------------------------------
  /** \brief  Parallel processor rank of the processor which owns this entity */
  unsigned owner_rank() const { return m_owner_rank ; }

  /** \brief  Parallel processors which share this entity. */
  const PairIterEntityProc & sharing() const { return m_sharing ; }

  //------------------------------------


private:

  const EntityKey       m_key ;        ///< Globally unique key
  std::vector<Relation> m_relation ;   ///< This entity's relationships
  Bucket *              m_bucket ;     ///< Bucket for the entity's field data
  unsigned              m_bucket_ord ; ///< Ordinal within the bucket
  unsigned              m_owner_rank ; ///< Owner processors' rank
  size_t                m_sync_count ; ///< Last membership change
  PairIterEntityProc    m_sharing ;    ///< Span of entity sharing vector

  ~Entity();
  explicit Entity( const EntityKey & arg_key );

  Entity(); ///< Default constructor not allowed
  Entity( const Entity & ); ///< Copy constructor not allowed
  Entity & operator = ( const Entity & ); ///< Assignment operator not allowed

#ifndef DOXYGEN_COMPILE
  friend class BulkData ;
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

//----------------------------------------------------------------------

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_Entity_hpp */

