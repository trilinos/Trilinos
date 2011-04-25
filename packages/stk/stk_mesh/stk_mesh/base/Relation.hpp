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

namespace stk {
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

  /** \brief  Destructor */
  ~Relation() {}

  /** \brief  Constructor */
  Relation() : m_raw_relation(), m_attribute(), m_entity(NULL) {}

  /** \brief  Copy Constructor */
  Relation( const Relation & r )
    : m_raw_relation( r.m_raw_relation ), m_attribute(r.m_attribute), m_entity(r.m_entity) {}

  /** \brief  Assignment operator */
  Relation & operator = ( const Relation & r )
  {
    if( this != &r ) {
      m_raw_relation = r.m_raw_relation ;
      m_attribute    = r.m_attribute ;
      m_entity       = r.m_entity ;
    }
    return *this ;
  }

  /** \brief  Construct a relation from a referenced entity,
   *          local identifier, kind, and converse flag.
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
  Entity * entity() const { return m_entity ; }

  /** \brief  Equality operator */
  bool operator == ( const Relation & r ) const
    { return m_raw_relation.value == r.m_raw_relation.value && m_entity == r.m_entity ; }

  /** \brief  Inequality operator */
  bool operator != ( const Relation & r ) const
    { return m_raw_relation.value != r.m_raw_relation.value || m_entity != r.m_entity ; }

  /** \brief  Ordering operator */
  bool operator < ( const Relation & r ) const ;

private:

  enum { entity_rank_ok = 1 / (!!(EntityKey::rank_digits == 8)) };
  enum {
    rank_digits = 8  ,
    id_digits   = 24 ,
    id_mask     = ~(0u) >> rank_digits
  };

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

  RawRelationType   m_raw_relation ;
  mutable attribute_type    m_attribute ;
  Entity          * m_entity ;
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
                              PartVector & induced_parts );

/** \brief  Induce entities' part membership based upon relationships
 *          between entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership( const Entity           & entity_from ,
                              const PartVector       & omit ,
                                    unsigned           entity_rank_to ,
                                    RelationIdentifier relation_identifier ,
                                    PartVector       & induced_parts );

/** \brief  Induce an entity's part membership based upon relationships
 *          from other entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership( const Entity     & entity ,
                              const PartVector & omit ,
                                    PartVector & induced_parts );


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

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_Relation_hpp */

