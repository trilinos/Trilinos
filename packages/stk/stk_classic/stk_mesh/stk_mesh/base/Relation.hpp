/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_Relation_hpp
#define stk_mesh_Relation_hpp

#include <iosfwd>
#include <limits>

#include <stk_mesh/base/EntityKey.hpp>

#include <boost/static_assert.hpp>

#ifdef SIERRA_MIGRATION

namespace stk_classic {
namespace mesh {
class Entity;
}
}
#endif

namespace stk_classic {
namespace mesh {

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
  Relation( Entity & entity , RelationIdentifier identifier );

  attribute_type   attribute() const { return m_attribute; }
  void set_attribute(attribute_type attr) const  { m_attribute = attr; }

  /** \brief  The encoded relation raw_relation_id */
  static raw_relation_id_type raw_relation_id( unsigned rank , unsigned id );

  /** \brief  The encoded relation raw_relation_id */
  raw_relation_id_type raw_relation_id() const { return m_raw_relation.value ; }

  /** \brief  The rank of the referenced entity */
  unsigned entity_rank() const ;

  /** \brief  The local relation identifier */
  RelationIdentifier identifier() const ;

  /** \brief  The referenced entity */
  Entity * entity() const { return m_target_entity ; }

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
  Entity               * m_target_entity ;

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
  Relation(Entity *obj, const unsigned relation_type, const unsigned ordinal, const unsigned orient = 0);

  Entity *getMeshObj() const {
    return entity();
  }

  void setMeshObj(Entity *object);

  RelationType getRelationType() const {
    return static_cast<RelationType>(attribute() >> fmwk_orientation_digits);
  }

  void setRelationType(RelationType relation_type) {
    set_attribute( (relation_type << fmwk_orientation_digits) | getOrientation() );
  }

  RelationIdentifier getOrdinal() const {
    return identifier();
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
RelationIdentifier Relation::identifier() const
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
  const std::vector<Entity*> & entities ,
        std::vector<Entity*> & entities_related );

/** \brief  Query which mesh entities have a relation
 *          to all of the input mesh entities of the given
 *          mesh rank.
 */
void get_entities_through_relations(
  const std::vector<Entity*> & entities ,
        EntityRank             entities_related_rank ,
        std::vector<Entity*> & entities_related );

//----------------------------------------------------------------------
/** \brief  Query if a member entity of the given entity type
 *          has an induced membership.
 */
bool membership_is_induced( const Part & part , unsigned entity_rank );

/** \brief  Induce entities' part membership based upon relationships
 *          between entities. Insert the result into 'induced_parts'.
 */
void induced_part_membership( Part & part ,
                              unsigned entity_rank_from ,
                              unsigned entity_rank_to ,
                              RelationIdentifier relation_identifier ,
                              OrdinalVector & induced_parts,
                              bool include_supersets=true);

/** \brief  Induce entities' part membership based upon relationships
 *          between entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership( const Entity           & entity_from ,
                              const OrdinalVector       & omit ,
                                    unsigned           entity_rank_to ,
                                    RelationIdentifier relation_identifier ,
                                    OrdinalVector       & induced_parts,
                                    bool include_supersets=true);

/** \brief  Induce an entity's part membership based upon relationships
 *          from other entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership( const Entity     & entity ,
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

} // namespace mesh
} // namespace stk_classic

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_Relation_hpp */

