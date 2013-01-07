/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_base_RelationAndEntity_hpp
#define stk_mesh_base_RelationAndEntity_hpp

#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Trace.hpp>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/util/PairIter.hpp>

#include <boost/static_assert.hpp>
#include <boost/range.hpp>
#include <boost/type_traits/is_pod.hpp>

#include <utility>
#include <vector>
#include <iosfwd>
#include <string>
#include <algorithm>

namespace stk {
namespace mesh {

class Entity;

namespace impl {
class EntityImpl;
}

/** \addtogroup stk_mesh_module
 *  \{
 */
//----------------------------------------------------------------------
/** \brief  A relation between two mesh entities with a
 *          relation <b> identifier </b> and <b> kind </b>.
 *
 *  Each entity owns a collection of relations to other entities.
 *  Each of these relations has a referenced entity, direction,
 *  relation identifier, and kind.
 *  - The direction specifies whether the relation is from the
 *    entity owning the relation to the referenced entity or
 *    conversely from the referenced entity to the owning entity.
 *  - The identifier provides a local numbering convention for relations,
 *    for example the local numbering convention for the nodes of an element.
 *  - Relations can be given a <em> kind </em> raw_relation_id to differentiate
 *    between topological relations (kind = 0) and other kinds of
 *    relationships such as constraints.
 *
 *  Relations are ordered by their
 *  - kind of relation,
 *  - type of the referenced entity,
 *  - if the relation is forward or converse,
 *  - local relation identifier, and then
 *  - referenced entity's identifier.
 */
class Relation {
public:
  typedef uint32_t raw_relation_id_type ;
  typedef uint32_t attribute_type;

  /** \brief  Constructor */
  Relation();

  /** \brief  Construct a relation from a referenced entity and local identifier
   */
  Relation( Entity entity , RelationIdentifier identifier );

  attribute_type   attribute() const { return m_attribute; }
  void set_attribute(attribute_type attr) const  { m_attribute = attr; }

  /** \brief  The encoded relation raw_relation_id */
  static raw_relation_id_type raw_relation_id( unsigned rank , unsigned id );

  /** \brief  The encoded relation raw_relation_id */
  raw_relation_id_type raw_relation_id() const { return m_raw_relation.value ; }

  /** \brief  The rank of the referenced entity */
  unsigned entity_rank() const ;

  /** \brief  The local relation identifier */
  RelationIdentifier relation_ordinal() const ;

  /** \brief  The referenced entity */
  Entity entity() const;

  /** \brief  Equality operator */
  bool operator == ( const Relation & r ) const;

  /** \brief  Inequality operator */
  bool operator != ( const Relation & r ) const
  { return !(*this == r); }

  /** \brief  Ordering operator */
  bool operator < ( const Relation & r ) const ;

private:

  enum {
    rank_digits = 8  ,
    id_digits   = 24 ,
    id_mask     = ~(0u) >> rank_digits
#ifdef SIERRA_MIGRATION
    ,
    fwmk_relation_type_digits = 8,
    fmwk_orientation_digits   = 24,
    fmwk_orientation_mask     = ~(0u) >> fwmk_relation_type_digits
#endif
  };

  BOOST_STATIC_ASSERT(( static_cast<unsigned>(EntityKey::rank_digits) == static_cast<unsigned>(rank_digits) ));

  union RawRelationType {
  public:
    raw_relation_id_type value ;

    struct {
      raw_relation_id_type identifier  : id_digits ;
      raw_relation_id_type entity_rank : rank_digits ;
    } normal_view ;

    struct {
      raw_relation_id_type entity_rank : rank_digits ;
      raw_relation_id_type identifier  : id_digits ;
    } reverse_view ;

    RawRelationType( raw_relation_id_type v ) : value(v) {}
    RawRelationType() : value(0) {}
    RawRelationType( const RawRelationType & rhs ) : value( rhs.value ) {}
    RawRelationType & operator = ( const RawRelationType & rhs )
      { value = rhs.value ; return *this ; }
  };

  RawRelationType        m_raw_relation ;
  mutable attribute_type m_attribute ;
  impl::EntityImpl*      m_target_entity ;

// Issue: Framework supports relation types (parent, child, etc) that STK_Mesh
// does not support, so these relations will have to be managed by framework
// until:
// A) STK_Mesh is rewritten to support these extra relation types
// B) An extra data-structure is added to manage these relation types and
// this data structure works together with STK_mesh under the generic API
// to manage all relation-types.
// C) We transition to using a mesh that supports everything framework
// supports and this problem goes away.
//
// The problem is that, with framework managing some relations and STK_Mesh
// managing others, the type of the relation descriptor is different depending
// on what type of relation you're dealing with. This can be addressed with
// templates, but this makes the code very ugly. Instead...
//
// Solution: Have framework and STK_Mesh use the same type as its relation_descriptor.
// The code below is designed to make this class compatible with the fmwk
// Relation class.
#ifdef SIERRA_MIGRATION
 public:
    /**
   * Predefined identifiers for mesh object relationship types.
   */
  enum RelationType {
    USES	= 0 ,
    USED_BY	= 1 ,
    CHILD	= 2 ,
    PARENT	= 3 ,
    EMBEDDED	= 0x00ff , // 4
    CONTACT	= 0x00ff , // 5
    AUXILIARY   = 0x00ff ,
    INVALID     = 10
  };

  enum {
    POLARITY_MASK       = 0x80,
    POLARITY_POSITIVE   = 0x80,
    POLARITY_NEGATIVE   = 0x00,
    POLARITY_IDENTITY   = 0x80
  };

  static bool polarity(unsigned orient) {
    return (orient & POLARITY_MASK) == POLARITY_POSITIVE;
  }

  static unsigned permutation(unsigned orient) {
    return orient & ~POLARITY_MASK;
  }

  /**
   * Construct filled-out relation, fmwk-style
   */
  Relation(Entity obj, const unsigned relation_type, const unsigned ordinal, const unsigned orient = 0);

  void setMeshObj(Entity object);

  RelationType getRelationType() const {
    return static_cast<RelationType>(attribute() >> fmwk_orientation_digits);
  }

  void setRelationType(RelationType relation_type) {
    set_attribute( (relation_type << fmwk_orientation_digits) | getOrientation() );
  }

  RelationIdentifier getOrdinal() const {
    return relation_ordinal();
  }

  void setOrdinal(RelationIdentifier ordinal) {
    m_raw_relation = Relation::raw_relation_id( entity_rank(), ordinal );
  }

  attribute_type getOrientation() const {
    return attribute() & fmwk_orientation_mask;
  }

  void setOrientation(attribute_type orientation) {
    set_attribute( (getRelationType() << fmwk_orientation_digits) | orientation );
  }

  /**
   * Query polarity of the related mesh object.
   * A 'true' polarity indicates that the related mesh object is aligned.
   * For element-edge or face-edge a 'true' polarity indicates that the
   * nodes of the edge are compatibly ordered with the nodes of the
   * element or faces.
   * For element-face a 'true' polarity indicates that the ordering
   * of face-nodes is compatible with the ordering of element nodes,
   * i.e. the face's normal defined by a clockwise ordering is outward.
   */
  bool polarity() const {
    return (getOrientation() & POLARITY_MASK) == POLARITY_POSITIVE;
  }

  unsigned permutation() const {
    return getOrientation() & ~POLARITY_MASK;
  }

private:
  bool has_fmwk_state() const { return getRelationType() != INVALID; }
#endif // SIERRA_MIGRATION
};

//----------------------------------------------------------------------

inline
Relation::raw_relation_id_type
Relation::raw_relation_id( unsigned rank , unsigned id )
{
  ThrowAssertMsg( id <= id_mask,
                  "For args rank " << rank << ", id " << id << ": " <<
                  "id " << " > id_mask=" << id_mask );

  return ( raw_relation_id_type(rank) << id_digits ) | id ;
}

inline
unsigned Relation::entity_rank() const
{ return m_raw_relation.value >> id_digits; }

inline
RelationIdentifier Relation::relation_ordinal() const
{ return unsigned( m_raw_relation.value & id_mask ); }

struct LessRelation {
  bool operator() ( const Relation & lhs , const Relation & rhs ) const
    { return lhs < rhs ; }

  bool operator() ( const Relation & lhs , Relation::raw_relation_id_type rhs ) const
    { return lhs.raw_relation_id() < rhs ; }
};

//----------------------------------------------------------------------
/** \brief  Query which mesh entities have a relation
 *          to all of the input mesh entities.
 */
void get_entities_through_relations(
  const std::vector<Entity> & entities ,
        std::vector<Entity> & entities_related );

/** \brief  Query which mesh entities have a relation
 *          to all of the input mesh entities of the given
 *          mesh rank.
 */
void get_entities_through_relations(
  const std::vector<Entity> & entities ,
        EntityRank             entities_related_rank ,
        std::vector<Entity> & entities_related );

//----------------------------------------------------------------------
/** \brief  Query if a member entity of the given entity type
 *          has an induced membership.
 */
bool membership_is_induced( const Part & part , unsigned entity_rank );

/** \brief  Induce entities' part membership based upon relationships
 *          between entities. Insert the result into 'induced_parts'.
 */
void induced_part_membership( const Part & part ,
                              unsigned entity_rank_from ,
                              unsigned entity_rank_to ,
                              RelationIdentifier relation_identifier ,
                              OrdinalVector & induced_parts,
                              bool include_supersets=true);

/** \brief  Induce entities' part membership based upon relationships
 *          between entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership( const Entity entity_from ,
                              const OrdinalVector       & omit ,
                                    unsigned           entity_rank_to ,
                                    RelationIdentifier relation_identifier ,
                                    OrdinalVector       & induced_parts,
                                    bool include_supersets=true);

/** \brief  Induce an entity's part membership based upon relationships
 *          from other entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership( const Entity entity ,
                              const OrdinalVector & omit ,
                                    OrdinalVector & induced_parts,
                                    bool include_supersets=true);


//----------------------------------------------------------------------

#if 0
/** \brief  Decode and print the relation raw_relation_id */
std::ostream &
print_relation( std::ostream & , Relation::raw__attr_type );

/** \brief  Decode and print the relation raw_relation_id and referenced entity key */
std::ostream &
print_relation( std::ostream & , const MetaData & ,
                Relation::raw__attr_type , EntityKey );
#endif

/** \brief  Print the relation raw_relation_ids and referenced entity's key */
std::ostream & operator << ( std::ostream & , const Relation & );

/** \} */

inline
Relation::Relation() :
  m_raw_relation(),
  m_attribute(),
  m_target_entity(NULL)
{
#ifdef SIERRA_MIGRATION
  setRelationType(INVALID);
#endif
}

inline
bool Relation::operator == ( const Relation & rhs ) const
{
  return m_raw_relation.value == rhs.m_raw_relation.value && m_target_entity == rhs.m_target_entity
#ifdef SIERRA_MIGRATION
    // compared fmwk state too
    && m_attribute == rhs.m_attribute
#endif
    ;
}

inline
bool same_specification(const Relation& lhs, const Relation& rhs)
{
#ifdef SIERRA_MIGRATION
  return  lhs.entity_rank()     == rhs.entity_rank() &&
          lhs.getRelationType() == rhs.getRelationType() &&
          lhs.getOrdinal()      == rhs.getOrdinal();
#else
  return  lhs.entity_rank()     == rhs.entity_rank();
#endif
}

//----------------------------------------------------------------------
// EntityImpl
//----------------------------------------------------------------------

namespace impl {

#ifdef SIERRA_MIGRATION

//On a *strictly temporary* basis, we need to stick the following
//fmwk stuff on an entity, just to help us through the sierra migration.
//Move along folks, there's nothing to see here.
struct fmwk_attributes {
  // Each member has an explanation of why it can't be handled by Entity.

  // Relations of this mesh object that can't be managed by STK such as PARENT/CHILD
  RelationVector        aux_relations;

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

  // Not a supported STK_Mesh concept, but fmwk requires it.
  unsigned              local_id;
};

#endif

/** \addtogroup stk_mesh_module
 * \{
 */

class EntityImpl {
public:

  EntityImpl( const EntityKey & arg_key );
  EntityImpl();

  // Exposed in external interface:
  EntityRank entity_rank() const { return stk::mesh::entity_rank( m_key ); }
  EntityId identifier() const { return stk::mesh::entity_id( m_key ); }
  const EntityKey & key() const { return m_key ; }
  PairIterRelation relations() const { return PairIterRelation(m_relation); }
  PairIterRelation relations( unsigned rank ) const ;
  PairIterRelation node_relations( ) const
  {
    RelationVector::const_iterator i = m_relation.begin();
    RelationVector::const_iterator e = m_relation.end();

    const Relation::raw_relation_id_type hi = Relation::raw_relation_id(1, 0);
    e = std::lower_bound( i , e , hi , LessRelation() );

    return PairIterRelation( i , e );
  }

  RelationVector::const_iterator node_relation(unsigned ordinal) const
  { return m_relation.begin() + ordinal; }

  Bucket & bucket() const
  {
    ThrowAssert(m_bucket); //don't want to return a reference to a null bucket
    return *m_bucket ;
  }

  Bucket* bucket_ptr() const
  {
    return m_bucket; // allow for NULL return value
  }

  bool is_bucket_valid() const { return m_bucket != NULL; }
  unsigned bucket_ordinal() const { return m_bucket_ord ; }
  unsigned owner_rank() const { return m_owner_rank ; }
  size_t synchronized_count() const { return m_sync_count ; }

  // The two relation methods below need to be called symmetically, ideally
  // through EntityRepository which will enforce the symmetry.

  bool destroy_relation( Entity e_to, const RelationIdentifier local_id);
  bool declare_relation( Entity e_to,
                         const RelationIdentifier local_id,
                         unsigned sync_count,
                         bool is_back_relation = false);

  // return true if entity was actually modified
  bool set_owner_rank( unsigned in_owner_rank );

  void set_sync_count( size_t sync_count );

  EntityState state() const { return m_state ; }

  void clear_state();

  /**
   * Mark this entity as modified (only changes from Unchanged
   * to Modified). Propagates the modification to higher-ranking
   * entities related to this entity. In other words, based on our
   * modification model, all entities that have modified_entity in their
   * closure must also be marked as modified.
   */
  void modified();

  /** \brief  Log that this entity was created as a parallel copy. */
  void created_parallel_copy();

  //set_key is only to be used for setting a key on a newly-constructed entity.
  void set_key(EntityKey key);

  void set_bucket_and_ordinal( Bucket * in_bucket, unsigned ordinal );

  //update_key is used to change the key for an entity that has been in use with
  //a different key.
  void update_key(EntityKey key);

//  RelationVector& rel_vec() {return m_relation;}
  void compress_relation_capacity();

 private:

  EntityKey               m_key ;        ///< Globally unique key
  RelationVector          m_relation ;   ///< This entity's relationships
  Bucket                * m_bucket ;     ///< Bucket for the entity's field data
  unsigned                m_bucket_ord ; ///< Ordinal within the bucket
  unsigned                m_owner_rank ; ///< Owner processors' rank
  size_t                  m_sync_count ; ///< Last membership change
  EntityState             m_state ;

#ifdef SIERRA_MIGRATION
 public:
  //
  // Members needed to support Fmwk_MeshObj API
  //
  fmwk_attributes m_fmwk_attrs;
 private:
#endif

//  EntityImpl( const EntityImpl & ); ///< Copy constructor not allowed
  EntityImpl & operator = ( const EntityImpl & ); ///< Assignment operator not allowed
};

//----------------------------------------------------------------------

/** \} */

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
// Entity
//----------------------------------------------------------------------


/** \addtogroup stk_mesh_module
 * \{
 */

namespace stk {
namespace adapt {
struct my_tuple_hash;
}
}

#ifdef SIERRA_MIGRATION

namespace stk {
namespace mesh {
typedef RelationVector::const_iterator   RelationIterator;
typedef boost::iterator_range<RelationIterator> RelationRange;
class Entity;
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
bool set_attributes( MeshBulkData& meshbulk, stk::mesh::Entity , const int , const MeshObjSharedAttr*, const int);
bool set_attributes( MeshBulkData& meshbulk, stk::mesh::Entity , const MeshObjSharedAttr*, const int);
void unset_shared_attr(stk::mesh::Entity );
}

namespace roster_only {
void destroy_meshobj(stk::mesh::Entity, MeshBulkData& meshbulk );
void set_shared_attr(stk::mesh::Entity , const MeshObjSharedAttr*);
}

const MeshObjSharedAttr * get_shared_attr(const stk::mesh::Entity );
bool insert_relation( stk::mesh::Entity , const stk::mesh::Relation::RelationType, stk::mesh::Entity , const unsigned, const unsigned, const bool, MeshBulkData &);
bool remove_relation(stk::mesh::Entity , const stk::mesh::RelationIterator, MeshBulkData &);
bool verify_relations(const stk::mesh::Entity );
}
}
#endif

namespace stk {
namespace mesh {

class Relation;
namespace impl {
class EntityRepository;
}

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
 *  When any of the above changes, the Entity's state may change
 */
class Entity {
public:

  bool is_valid() const { return m_entityImpl != NULL; }

  /** \brief  Query the current state of the entity */
  EntityState state() const { return m_entityImpl->state(); }

  /** \brief  The rank of this entity. */
  EntityRank entity_rank() const { return m_entityImpl->entity_rank(); }

  /** \brief  Identifier for this entity which is globally unique
   *          for a given entity type.
   */
  EntityId identifier() const { return m_entityImpl->identifier(); }

  /** \brief  The globally unique key ( entity type + identifier )
   *          of this entity.
   */
  const EntityKey & key() const { return m_entityImpl->key(); }

  /** \brief  The bucket which holds this mesh entity's field data */
  Bucket & bucket() const { return m_entityImpl->bucket(); }
  Bucket * bucket_ptr() const { return m_entityImpl->bucket_ptr(); }

  /** \brief  The ordinal for this entity within its bucket. */
  unsigned bucket_ordinal() const { return m_entityImpl->bucket_ordinal(); }

  /** \brief  The mesh bulk data synchronized_count when this entity's
   *          part membership was most recently modified.
   *
   *  If ( mesh.synchronized_state() == false &&
   *       mesh.synchronized_count() == entity.synchronized_count() )
   *  then entity was modified during this modification phase.
   */
  size_t synchronized_count() const { return m_entityImpl->synchronized_count(); }

  //------------------------------------
  /** \brief  All \ref stk::mesh::Relation "Entity relations"
   *          for which this entity is a member. The relations are ordered
   *          from lowest entity-rank to highest entity-rank.
   */
  PairIterRelation relations() const { return m_entityImpl->relations(); }

  /** \brief  \ref stk::mesh::Relation "Entity relations" for which this
   *          entity is a member, the other entity is of a given type.
   */
  PairIterRelation relations( EntityRank type ) const { return m_entityImpl->relations(type); }
  PairIterRelation node_relations() const { return m_entityImpl->node_relations(); }

#ifdef SIERRA_MIGRATION
  RelationIterator node_relation(unsigned ordinal) const { return m_entityImpl->node_relation(ordinal); }
#endif

  //------------------------------------
  /** \brief  Parallel processor rank of the processor which owns this entity.
   * IMPORTANT NOTE: when sierra framework is handling communication for stk-mesh,
   * this method will always return zero. If using stk-mesh beneath sierra framework,
   * use the method owner_processor_rank() instead (defined below). */
  unsigned owner_rank() const { return m_entityImpl->owner_rank(); }

//  RelationVector& rel_vec() { return m_entityImpl->rel_vec(); }
  void compress_relation_capacity();

  bool operator==(Entity entity) const;

  bool operator!=(Entity entity) const;

  bool operator<(Entity entity) const;

private:

  void set_key(const EntityKey& arg_key) { m_entityImpl->set_key(arg_key); }

#ifndef DOXYGEN_COMPILE
  friend class impl::EntityRepository ;
  friend class impl::EntityImpl ;
  friend class Relation;
  friend struct adapt::my_tuple_hash;
#endif /* DOXYGEN_COMPILE */

 public:
  impl::EntityImpl* m_entityImpl;

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
  friend const sierra::Fmwk::MeshObjSharedAttr * sierra::Fmwk::get_shared_attr(const Entity );
  friend bool sierra::Fmwk::detail::set_attributes( sierra::Fmwk::MeshBulkData&, Entity , const int, const sierra::Fmwk::MeshObjSharedAttr *, const int);
  friend bool sierra::Fmwk::detail::set_attributes( sierra::Fmwk::MeshBulkData&, Entity , const sierra::Fmwk::MeshObjSharedAttr *, const int);
  friend void sierra::Fmwk::detail::unset_shared_attr(Entity );
  friend bool sierra::Fmwk::insert_relation( Entity , const stk::mesh::Relation::RelationType, Entity , const unsigned, const unsigned, const bool, sierra::Fmwk::MeshBulkData &);
  friend bool sierra::Fmwk::remove_relation(Entity , const stk::mesh::RelationIterator, sierra::Fmwk::MeshBulkData &);
  friend bool sierra::Fmwk::verify_relations(const Entity );
  friend void sierra::Fmwk::roster_only::destroy_meshobj(stk::mesh::Entity, sierra::Fmwk::MeshBulkData& meshbulk );
  friend void sierra::Fmwk::roster_only::set_shared_attr(stk::mesh::Entity , const sierra::Fmwk::MeshObjSharedAttr*);

  typedef unsigned DerivedType; ///< Derived type identifier, the admissible values may be extended

  /**
   * Predefined derived type identifiers. Should match rank enum in MetaData, but we don't
   * want to add a header dependency on MetaData.
   */
  enum ObjectTypeEnum {
    NODE       = 0,
    EDGE       = 1,
    FACE       = 2,
    ELEMENT    = 3,
    CONSTRAINT = 4,
    NUM_TYPES  = 5,
    BASE_CLASS = 0x00ff
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

    m_entityImpl->m_fmwk_attrs.global_id     = id;
    m_entityImpl->m_fmwk_attrs.shared_attr   = attr;
    m_entityImpl->m_fmwk_attrs.owner         = owner;
    m_entityImpl->m_fmwk_attrs.connect_count = 0;
    m_entityImpl->m_fmwk_attrs.local_id      = sierra::Fmwk::INVALID_LOCAL_ID;

    if (attr->locally_owned() && owner_processor_rank() == -1) {
      m_entityImpl->m_fmwk_attrs.owner = parallel_rank;
      ThrowAssert(owner_processor_rank() < parallel_size);
    }

    internal_verify_initialization_invariant();
  }

  /**
   * Get global identifier
   */
  int global_id() const {
    return m_entityImpl->m_fmwk_attrs.global_id;
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
  inline unsigned local_id() const { return m_entityImpl->m_fmwk_attrs.local_id; }
  void set_local_id(unsigned int l_id) { m_entityImpl->m_fmwk_attrs.local_id = l_id; }

  int owner_processor_rank() const { return m_entityImpl->m_fmwk_attrs.owner; }
  void set_owner_processor_rank(int owner) { m_entityImpl->m_fmwk_attrs.owner = owner; }

  /**
   * Number of connections to this mesh object.
   * The connection count should be used to track inter-mesh object
   * relationships to this mesh object that are not expressed by a
   * corresponding relationship owned by this mesh object.
   */
  unsigned size_connection() const {
    return m_entityImpl->m_fmwk_attrs.connect_count;
  }

  /**
   * Increment the connection count.
   */
  unsigned inc_connection() {
    ++m_entityImpl->m_fmwk_attrs.connect_count;
    ThrowAssert(m_entityImpl->m_fmwk_attrs.connect_count /* Did not roll over */);
    return m_entityImpl->m_fmwk_attrs.connect_count;
  }

  /**
   * Decrement the connection count.
   */
  unsigned dec_connection() {
    ThrowAssert(m_entityImpl->m_fmwk_attrs.connect_count /* Will not roll-under */);
    --m_entityImpl->m_fmwk_attrs.connect_count;
    return m_entityImpl->m_fmwk_attrs.connect_count;
  }

  RelationIterator aux_relation_begin() const { return m_entityImpl->m_fmwk_attrs.aux_relations.begin(); }
  RelationIterator aux_relation_end() const { return m_entityImpl->m_fmwk_attrs.aux_relations.end(); }

  RelationVector& aux_relations() { return m_entityImpl->m_fmwk_attrs.aux_relations; }

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

  void set_shared_attr(const void* attr) { m_entityImpl->m_fmwk_attrs.shared_attr = attr; }
  const void* get_shared_attr() const { return m_entityImpl->m_fmwk_attrs.shared_attr; }

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
    ThrowAssert( !(m_entityImpl->m_fmwk_attrs.global_id < 0 && key().id() == static_cast<uint64_t>(m_entityImpl->m_fmwk_attrs.global_id)) &&
                 !(m_entityImpl->m_fmwk_attrs.global_id > 0 && key().id() != static_cast<uint64_t>(m_entityImpl->m_fmwk_attrs.global_id)) );

  }

  unsigned stk_entity_rank() const { return key().rank(); }

  bool internal_is_handled_generically(const Relation::RelationType relation_type) const
  {
    return relation_type == Relation::USES || relation_type == Relation::USED_BY;
  }
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
  bool operator()(const Entity lhs, const Entity rhs) const
  {
    const EntityKey lhs_key = lhs.is_valid() ? lhs.key() : EntityKey();
    const EntityKey rhs_key = rhs.is_valid() ? rhs.key() : EntityKey();
    return lhs_key < rhs_key;
  }

  bool operator()(const Entity lhs, const EntityKey & rhs) const
  {
    const EntityKey lhs_key = lhs.is_valid() ? lhs.key() : EntityKey();
    return lhs_key < rhs;
  }

  bool operator()( const EntityProc & lhs, const EntityProc & rhs) const;

  bool operator()( const EntityProc & lhs, const Entity rhs) const;

  bool operator()( const EntityProc & lhs, const EntityKey & rhs) const;

}; //class EntityLess

class EntityEqual
{
public:
  bool operator()(const Entity lhs, const Entity rhs) const
  {
    const EntityKey lhs_key = lhs.is_valid() ? lhs.key() : EntityKey();
    const EntityKey rhs_key = rhs.is_valid() ? rhs.key() : EntityKey();
    return lhs_key == rhs_key;
  }
};

std::string print_entity_key(const Entity entity);

std::string print_entity_key(const Entity entity);

/** \} */

inline
bool Entity::operator==(Entity entity) const
{ return EntityEqual()(*this, entity); }

inline
bool Entity::operator!=(Entity entity) const
{ return !(EntityEqual()(*this, entity)); }

inline
bool Entity::operator<(Entity entity) const
{ return EntityLess()(*this, entity); }

inline
Relation::Relation( Entity ent , RelationIdentifier id )
  : m_raw_relation( Relation::raw_relation_id( ent.entity_rank() , id ) ),
    m_target_entity(ent.m_entityImpl)
{
#ifdef SIERRA_MIGRATION
  setRelationType(INVALID);
#endif
}

inline
Entity Relation::entity() const
{ Entity rv = {m_target_entity};
  return rv;
}

inline
bool Relation::operator < ( const Relation & rhs ) const
{
  bool result = false;

#ifdef SIERRA_MIGRATION
  if (entity_rank() != rhs.entity_rank()) {
    result = entity_rank() < rhs.entity_rank();
  }
  else if (getRelationType() != rhs.getRelationType()) {
    result = getRelationType() < rhs.getRelationType();
  }
  else if (relation_ordinal() != rhs.relation_ordinal()) {
    result = relation_ordinal() < rhs.relation_ordinal();
  }
#else
  if ( m_raw_relation.value != rhs.m_raw_relation.value ) {
    result = m_raw_relation.value < rhs.m_raw_relation.value ;
  }
#endif
  else {
    const EntityKey lhs_key = m_target_entity     ? m_target_entity->key()     : EntityKey();
    const EntityKey rhs_key = rhs.m_target_entity ? rhs.m_target_entity->key() : EntityKey();
    result = lhs_key < rhs_key ;
  }
  return result ;
}

#ifdef SIERRA_MIGRATION

inline
Relation::Relation(Entity obj, const unsigned relation_type, const unsigned ordinal, const unsigned orient)
  :
  m_raw_relation( Relation::raw_relation_id( obj.entity_rank(), ordinal )),
  m_attribute( (relation_type << fmwk_orientation_digits) | orient ),
  m_target_entity(obj.m_entityImpl)
{
  ThrowAssertMsg( orient <= fmwk_orientation_mask,
                  "orientation " << orient << " exceeds maximum allowed value");
}

inline
void Relation::setMeshObj(Entity object)
{
  if (object.is_valid()) {
    m_raw_relation = Relation::raw_relation_id( object.entity_rank(), relation_ordinal() );
  }
  m_target_entity = object.m_entityImpl;
}

#endif

inline
void Entity::reserve_relation(const unsigned num)
{
  if (num == 0 && aux_relations().empty()) {
    RelationVector tmp;
    aux_relations().swap(tmp); // clear memory of m_relations.
  }
  else {
    aux_relations().reserve(num);
  }
}

inline
void Entity::compress_relation_capacity()
{
  m_entityImpl->compress_relation_capacity();
  if (!m_entityImpl->m_fmwk_attrs.aux_relations.empty()) {
    RelationVector tmp(m_entityImpl->m_fmwk_attrs.aux_relations);
    tmp.swap(m_entityImpl->m_fmwk_attrs.aux_relations);
  }
}

inline
void Entity::erase_and_clear_if_empty(RelationIterator rel_itr)
{
  ThrowAssert(!internal_is_handled_generically(rel_itr->getRelationType()));

  RelationVector& aux_rels = m_entityImpl->m_fmwk_attrs.aux_relations;
  aux_rels.erase(aux_rels.begin() + (rel_itr - aux_rels.begin())); // Need to convert to non-const iterator

  if (aux_rels.empty()) {
    reserve_relation(0);
  }
}

inline
bool EntityLess::operator()( const EntityProc & lhs, const EntityProc & rhs) const
{
  const EntityKey lhs_key = lhs.first.is_valid() ? lhs.first.key() : EntityKey() ;
  const EntityKey rhs_key = rhs.first.is_valid() ? rhs.first.key() : EntityKey() ;
  return lhs_key != rhs_key ? lhs_key < rhs_key : lhs.second < rhs.second ;
}

inline
bool EntityLess::operator()( const EntityProc & lhs, const Entity rhs) const
{
  const EntityKey lhs_key = lhs.first.is_valid() ? lhs.first.key() : EntityKey();
  const EntityKey rhs_key = rhs.is_valid()       ? rhs.key()       : EntityKey();
  return lhs_key < rhs.key();
}

inline
bool EntityLess::operator()( const EntityProc & lhs, const EntityKey & rhs) const
{
  const EntityKey lhs_key = lhs.first.is_valid() ? lhs.first.key() : EntityKey();
  return lhs_key < rhs ;
}

inline
std::ostream& operator<<(std::ostream& out, Entity entity)
{
  return out << entity.identifier();
}

inline
impl::EntityImpl::EntityImpl( const EntityKey & arg_key )
  : m_key(arg_key),
    m_relation(),
    m_bucket( NULL ),
    m_bucket_ord(0),
    m_owner_rank(0),
    m_sync_count(0),
    m_state( Created )
{
  TraceIfWatching("stk::mesh::impl::EntityImpl::EntityImpl", LOG_ENTITY, arg_key);
}

inline
impl::EntityImpl::EntityImpl()
  : m_key(),
    m_relation(),
    m_bucket( NULL ),
    m_bucket_ord(0),
    m_owner_rank(0),
    m_sync_count(0),
    m_state( Created )
{
}

inline
void impl::EntityImpl::set_bucket_and_ordinal( Bucket * in_bucket, unsigned ordinal )
{
  TraceIfWatching("stk::mesh::impl::EntityImpl::set_bucket_and_ordinal", LOG_ENTITY, key());

  m_bucket = in_bucket;
  m_bucket_ord = ordinal;
}

// return true if entity was actually modified
inline
bool impl::EntityImpl::set_owner_rank( unsigned in_owner_rank )
{
  TraceIfWatching("stk::mesh::impl::EntityImpl::set_owner_rank", LOG_ENTITY, key());

  if ( in_owner_rank != m_owner_rank ) {
    m_owner_rank = in_owner_rank;
    return true;
  }
  return false;
}

inline
void impl::EntityImpl::set_sync_count( size_t sync_count )
{
  TraceIfWatching("stk::mesh::impl::EntityImpl::set_sync_count", LOG_ENTITY, key());

  m_sync_count = sync_count;
}

inline
void impl::EntityImpl::clear_state()
{
  TraceIfWatching("stk::mesh::impl::EntityImpl::clear_state", LOG_ENTITY, key());

  m_state = Unchanged;
}

inline
size_t hash_value( Entity entity) {
  return hash_value(entity.key());
}

} // namespace mesh
} // namespace stk

#endif /* stk_mesh_RelationAndEntity_hpp */
