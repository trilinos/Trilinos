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
#include <iosfwd>
#include <string>

#include <stk_mesh/base/Types.hpp>

#include <stk_mesh/baseImpl/EntityImpl.hpp>

#include <boost/range.hpp>

#ifdef SIERRA_MIGRATION
#include <stk_mesh/base/Relation.hpp>

namespace stk {
namespace mesh {
typedef RelationVector::const_iterator   RelationIterator;
typedef boost::iterator_range<RelationIterator> RelationRange;
class Entity;

//On a *strictly temporary* basis, we need to stick the following
//fmwk stuff on an entity, just to help us through the sierra migration.
//Move along folks, there's nothing to see here.
struct fmwk_attributes {
  // Each member has an explanation of why it can't be handled by Entity.

  // Relations of this mesh object that can't be managed by STK such as PARENT/CHILD
  RelationVector aux_relations;

  // Not a supported STK_Mesh concept
  const void*           shared_attr;

  // Cannot just use the id embedded in the entity-key because of negative global-ids
  // for temporaries.
  int                   global_id;

  // Can't use STK_Mesh's notion of this entity's owner because STK_Mesh is being used
  // in serial mode and therefore every process will thinks it owns everything.
  int                   owner;

  // Not a supported STK_Mesh concept
  unsigned short        connect_count;
};

}
}

namespace sierra {
namespace Fmwk {

class MeshObjRoster;
class MeshObjSharedAttr;
class MeshBulkData;

extern const unsigned int INVALID_LOCAL_ID;
extern const stk::mesh::RelationIterator INVALID_RELATION_ITR;

namespace detail {
bool set_attributes( stk::mesh::Entity & , const int , const MeshObjSharedAttr*, const int);
bool set_attributes( stk::mesh::Entity & , const MeshObjSharedAttr*, const int);
void unset_shared_attr(stk::mesh::Entity& );
}

namespace roster_only {
void destroy_meshobj(stk::mesh::Entity*);
void set_shared_attr(stk::mesh::Entity&, const MeshObjSharedAttr*);
}

const MeshObjSharedAttr * get_shared_attr(const stk::mesh::Entity&);
bool insert_relation( stk::mesh::Entity * const, const stk::mesh::Relation::RelationType, stk::mesh::Entity * const, const unsigned, const unsigned, const bool, MeshBulkData &);
bool remove_relation(stk::mesh::Entity &, const stk::mesh::RelationIterator, MeshBulkData &);
bool verify_relations(const stk::mesh::Entity &);
}
}
#endif

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
  ~Entity() {}

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
  Bucket * bucket_ptr() const { return m_entityImpl.bucket_ptr(); }

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
  PairIterRelation node_relations() const { return m_entityImpl.node_relations(); }

#ifdef SIERRA_MIGRATION
  RelationIterator node_relation(unsigned ordinal) const { return m_entityImpl.node_relation(ordinal); }
#endif

  //------------------------------------
  /** \brief  Parallel processor rank of the processor which owns this entity */
  unsigned owner_rank() const {
#ifdef SIERRA_MIGRATION
    return 0;
#else
    return m_entityImpl.owner_rank();
#endif
  }

  /** \brief  Parallel processes which share this entity. */
  PairIterEntityComm sharing() const { return m_entityImpl.sharing(); }

  /** \brief  Complete communicaiton list for this entity */
  PairIterEntityComm comm() const { return m_entityImpl.comm(); }

  /** \brief  Subset communicaiton list for this entity */
  PairIterEntityComm comm( const Ghosting & sub ) const { return m_entityImpl.comm( sub ); }

//  RelationVector& rel_vec() { return m_entityImpl.rel_vec(); }
  void compress_relation_capacity();

private:

  void set_key(const EntityKey& arg_key) { m_entityImpl.set_key(arg_key); }

  impl::EntityImpl m_entityImpl;

  explicit Entity( const EntityKey & arg_key );

  Entity();
  Entity & operator = ( const Entity & ); ///< Assignment operator not allowed

#ifndef DOXYGEN_COMPILE
  friend class impl::EntityRepository ;
  friend class impl::EntityImpl ;
#endif /* DOXYGEN_COMPILE */



// Issue: We began the migration of Fmwk MeshObj with an implementation that
// had a stk entity pointer underneath the MeshObj that could be used to manage
// the state that stk entities can manage while leaving the rest to MeshObj.
// This proved to be a performance killer since many operations required an
// extra dereference (the entity pointer) compared to before, causing lots of
// additional cache misses even for simple operations.
//
// The solution is similar to what we did with Relation; we will use the preprocessor
// to add in the Fmwk pieces of the class API when the class is being compiled
// with Sierra.
#ifdef SIERRA_MIGRATION
 public:
  friend class sierra::Fmwk::MeshObjRoster;
  // These are free functions to facilitate the stk migration:
  friend const sierra::Fmwk::MeshObjSharedAttr * sierra::Fmwk::get_shared_attr(const Entity &);
  friend bool sierra::Fmwk::detail::set_attributes( Entity &, const int, const sierra::Fmwk::MeshObjSharedAttr *, const int);
  friend bool sierra::Fmwk::detail::set_attributes( Entity &, const sierra::Fmwk::MeshObjSharedAttr *, const int);
  friend void sierra::Fmwk::detail::unset_shared_attr(Entity &);
  friend bool sierra::Fmwk::insert_relation( Entity * const, const stk::mesh::Relation::RelationType, Entity * const, const unsigned, const unsigned, const bool, sierra::Fmwk::MeshBulkData &);
  friend bool sierra::Fmwk::remove_relation(Entity &, const stk::mesh::RelationIterator, sierra::Fmwk::MeshBulkData &);
  friend bool sierra::Fmwk::verify_relations(const Entity &);
  friend void sierra::Fmwk::roster_only::destroy_meshobj(stk::mesh::Entity*);
  friend void sierra::Fmwk::roster_only::set_shared_attr(stk::mesh::Entity&, const sierra::Fmwk::MeshObjSharedAttr*);

  typedef unsigned DerivedType; ///< Derived type identifier, the admissible values may be extended

  /**
   * Predefined derived type identifiers.
   */
  enum ObjectTypeEnum {
    NODE = 0, EDGE = 1, FACE = 2, ELEMENT = 3, CONSTRAINT = 4, NUM_TYPES = 5, BASE_CLASS = 0x00ff
  };
  static std::string TypeToString (ObjectTypeEnum type);



  template <class SharedAttr>
  void init_fmwk(
    const int         id,
    const SharedAttr* attr,
    const int         owner,
    const int         parallel_rank,
    const int         parallel_size)
  {
    ThrowAssertMsg(aux_relations().capacity() == 0, "Leftover memory found in relation vector");

    m_fmwk_attrs->global_id     = id;
    m_fmwk_attrs->shared_attr   = attr;
    m_fmwk_attrs->owner         = owner;
    m_fmwk_attrs->connect_count = 0;
    m_local_id                  = sierra::Fmwk::INVALID_LOCAL_ID;

    if (attr->locally_owned() && owner_processor_rank() == -1) {
      m_fmwk_attrs->owner = parallel_rank;
      ThrowAssert(owner_processor_rank() < parallel_size);
    }

    internal_verify_initialization_invariant();
  }

  /**
   * Get global identifier
   */
  int global_id() const {
    return m_fmwk_attrs->global_id;
  }

  /**
     At MeshObj construction, local_id is set to INVALID_LOCAL_ID.
     local_id is set to a valid value when MeshObjRoster::compact_roster is called,
     which occurs when optimize_roster is called, or when any of the commit_global_*
     methods are called.

     When set to valid values, local_ids are contiguous and start at 0 on
     the local MPI process for MeshObjs of each derived-type.
     They are in the same order as the global_ids.
     i.e., if meshobjA.local_id() < meshobjB.local_id(),
     then meshobjA.global_id() < meshobjB.global_id() if meshobjA and meshobjB have
     the same derived-type.
  */
  inline unsigned local_id() const { return m_local_id; }
  void set_local_id(unsigned int l_id) { m_local_id = l_id; }

  int owner_processor_rank() const { return m_fmwk_attrs->owner; }
  void set_owner_processor_rank(int owner) { m_fmwk_attrs->owner = owner; }

  /**
   * Number of connections to this mesh object.
   * The connection count should be used to track inter-mesh object
   * relationships to this mesh object that are not expressed by a
   * corresponding relationship owned by this mesh object.
   */
  unsigned size_connection() const {
    return m_fmwk_attrs->connect_count;
  }

  /**
   * Increment the connection count.
   */
  unsigned inc_connection() {
    ++m_fmwk_attrs->connect_count;
    ThrowAssert(m_fmwk_attrs->connect_count /* Did not roll over */);
    return m_fmwk_attrs->connect_count;
  }

  /**
   * Decrement the connection count.
   */
  unsigned dec_connection() {
    ThrowAssert(m_fmwk_attrs->connect_count /* Will not roll-under */);
    --m_fmwk_attrs->connect_count;
    return m_fmwk_attrs->connect_count;
  }

  RelationIterator aux_relation_begin() const { return m_fmwk_attrs->aux_relations.begin(); }
  RelationIterator aux_relation_end() const { return m_fmwk_attrs->aux_relations.end(); }

  RelationVector& aux_relations() { return m_fmwk_attrs->aux_relations; }

  /**
   * iterator to first relationship within the collection that mananges
   * relations of the given type.
   */
  RelationIterator internal_begin_relation(const Relation::RelationType relation_type) const {
    if (internal_is_handled_generically(relation_type)) {
      return relations().first;
    }
    else {
      return aux_relation_begin();
    }
  }

  /**
   * iterator to 'end' relationship within the collection that mananges
   * relations of the given type.
   */
  RelationIterator internal_end_relation(const Relation::RelationType relation_type) const {
    if (internal_is_handled_generically(relation_type)) {
      return relations().second;
    }
    else {
      return aux_relation_end();
    }
  }

  void set_shared_attr(const void* attr) { m_fmwk_attrs->shared_attr = attr; }
  const void* get_shared_attr() const { return m_fmwk_attrs->shared_attr; }

  // TODO: Refactor clients so that this method is no longer needed
  void set_relation_orientation(RelationIterator rel, unsigned orientation);

 private:

  /**
   * Reserve storage for a total of 'num' relationships
   * If 'num' is more that the current number of relationships
   * then additional storage is allocated but not used.
   * If 'num' is less that the current number of relationships
   * then no action is taken.
   */
  void reserve_relation(const unsigned num);

  /**
   * Update the back-relation
   * If back_rel_flag then insert the back-pointer if it does not exist
   * else remove the back-pointer and increment the related object's counter.
   *
   * In general, clients should not be calling this. MeshObj should be managing
   * its own back relations.
   */
  bool update_relation(const RelationIterator ir, const bool back_rel_flag) const;

  RelationIterator find_relation(const Relation& relation) const;

  void erase_and_clear_if_empty(RelationIterator rel_itr);

  void internal_verify_meshobj_invariant() const;

  void internal_swap_in_real_entity(const int globalId);

  void internal_verify_initialization_invariant() {
    // If this MeshObj has a proper ID (fully initialized), then the id should match
    // the id in the entity-key; otherwise they should not match.
    ThrowAssert( !(m_fmwk_attrs->global_id < 0 && key().id() == static_cast<uint64_t>(m_fmwk_attrs->global_id)) &&
                 !(m_fmwk_attrs->global_id > 0 && key().id() != static_cast<uint64_t>(m_fmwk_attrs->global_id)) );

  }

  unsigned stk_entity_rank() const { return key().rank(); }

  bool internal_is_handled_generically(const Relation::RelationType relation_type) const
  {
    return relation_type == Relation::USES || relation_type == Relation::USED_BY;
  }

  //
  // Members needed to support Fmwk_MeshObj API
  //

  fmwk_attributes* m_fmwk_attrs;

  // Not a supported STK_Mesh concept, but fmwk requires it. We have it
  // here instead of in m_fmwk_attrs for performance reasons (saves a
  // memory hop when accessing local_id).
  unsigned m_local_id;
#endif
};

#ifdef SIERRA_MIGRATION

inline
Relation::RelationType
back_relation_type(const Relation::RelationType relType)
{
  /* %TRACE[NONE]% */  /* %TRACE% */
  switch(relType) {
  case Relation::USES:
    return Relation::USED_BY;
  case Relation::USED_BY:
    return Relation::USES;
  case Relation::CHILD:
    return Relation::PARENT;
  case Relation::PARENT:
    return Relation::CHILD;
  default:
    return relType;
  }
}

// Made publicly available so that MeshObj.C can use it
template <class Iterator>
bool
verify_relation_ordering(Iterator begin, Iterator end)
{
  for (Iterator itr = begin; itr != end; ) {
    Iterator prev = itr;
    ++itr;
    if (itr != end) {

      if (itr->entity_rank() < prev->entity_rank()) {
        return false ;
      }

      if (itr->entity_rank() == prev->entity_rank()) {

        if (itr->getRelationType() < prev->getRelationType()) {
          return false ;
        }

        if (itr->getRelationType() == prev->getRelationType()) {

          if (itr->getOrdinal() < prev->getOrdinal()) {
            return false ;
          }
        }
      }
    }
  }
  return true ;
}

template <class Range>
bool
verify_relation_ordering(Range range)
{
  return verify_relation_ordering(boost::const_begin(range), boost::const_end(range));
}

#endif

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

std::string print_entity_key(const Entity& entity);

std::string print_entity_key(const Entity* entity);

inline
Entity::Entity()
  : m_entityImpl()
#ifdef SIERRA_MIGRATION
   , m_fmwk_attrs(NULL)
#endif
{}

inline
Entity::Entity( const EntityKey & arg_key )
  : m_entityImpl( arg_key )
#ifdef SIERRA_MIGRATION
   , m_fmwk_attrs(NULL)
#endif
{}

/** \} */

} // namespace mesh
} // namespace stk

#endif /* stk_mesh_base_Entity_hpp */
