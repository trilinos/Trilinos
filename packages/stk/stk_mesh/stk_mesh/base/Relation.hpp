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
 *  - Relations can be given a <em> kind </em> attribute to differentiate
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
  typedef uintptr_t raw_attr_type ;

  /** \brief  Destructor */
  ~Relation() {}

  /** \brief  Constructor */
  Relation() : m_attr(), m_entity(NULL) {}

  /** \brief  Copy Constructor */
  Relation( const Relation & r )
    : m_attr( r.m_attr ), m_entity(r.m_entity) {}

  /** \brief  Assignment operator */
  Relation & operator = ( const Relation & r )
  {
    if( this != &r ) {
      m_attr   = r.m_attr ;
      m_entity = r.m_entity ;
    }
    return *this ;
  }

  /** \brief  Construct a relation from an encoded relation attribute
   *          and a referenced entity.
   */
  Relation( raw_attr_type attr , Entity & entity );

  /** \brief  Construct a relation from a referenced entity,
   *          local identifier, kind, and converse flag.
   */
  Relation( Entity & entity , RelationIdentifier identifier );

  /** \brief  The encoded relation attribute */
  static raw_attr_type attribute( unsigned rank , unsigned id );

  /** \brief  The encoded relation attribute */
  raw_attr_type attribute() const { return m_attr.value ; }

  /** \brief  The rank of the referenced entity */
  unsigned entity_rank() const ;

  /** \brief  The local relation identifier */
  RelationIdentifier identifier() const ;

  /** \brief  The referenced entity */
  Entity * entity() const { return m_entity ; }

  /** \brief  Equality operator */
  bool operator == ( const Relation & r ) const
    { return m_attr.value == r.m_attr.value && m_entity == r.m_entity ; }

  /** \brief  Inequality operator */
  bool operator != ( const Relation & r ) const
    { return m_attr.value != r.m_attr.value || m_entity != r.m_entity ; }

  /** \brief  Ordering operator */
  bool operator < ( const Relation & r ) const ;

private:

  enum {
    attr_digits   = std::numeric_limits<raw_attr_type>::digits ,
    uint_digits   = std::numeric_limits<unsigned>::digits ,
    rank_digits   = EntityKey::rank_digits ,
    id_max_digits = attr_digits - rank_digits ,
    id_digits     = uint_digits - rank_digits ,
    id_mask       = ~(0u) >> ( uint_digits - id_digits ),
    rank_shift    = id_max_digits
  };

  union AttrType {
  public:
    raw_attr_type value ;

    struct {
      raw_attr_type identifier  : id_digits ;
      raw_attr_type entity_rank : rank_digits ;
    } normal_view ;

    struct {
      raw_attr_type entity_rank : rank_digits ;
      raw_attr_type identifier  : id_digits ;
    } reverse_view ;

    AttrType( raw_attr_type v ) : value(v) {}
    AttrType() : value(0) {}
    AttrType( const AttrType & rhs ) : value( rhs.value ) {}
    AttrType & operator = ( const AttrType & rhs )
      { value = rhs.value ; return *this ; }
  };

  AttrType m_attr ;
  Entity * m_entity ;
};

//----------------------------------------------------------------------

inline
unsigned Relation::entity_rank() const
{ return unsigned( m_attr.value >> rank_shift ); }

inline
RelationIdentifier Relation::identifier() const
{ return unsigned( m_attr.value & id_mask ); }

struct LessRelation {
  bool operator() ( const Relation & lhs , const Relation & rhs ) const
    { return lhs < rhs ; }

  bool operator() ( const Relation & lhs , Relation::raw_attr_type rhs ) const
    { return lhs.attribute() < rhs ; }
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
                                    PartVector       & entity_to_parts );

/** \brief  Induce an entity's part membership based upon relationships
 *          from other entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership( const Entity     & entity ,
                              const PartVector & omit ,
                                    PartVector & induced );


//----------------------------------------------------------------------

#if 0
/** \brief  Decode and print the relation attribute */
std::ostream &
print_relation( std::ostream & , Relation::raw__attr_type );

/** \brief  Decode and print the relation attribute and referenced entity key */
std::ostream &
print_relation( std::ostream & , const MetaData & ,
                Relation::raw__attr_type , EntityKey );
#endif

/** \brief  Print the relation attributes and referenced entity's key */
std::ostream & operator << ( std::ostream & , const Relation & );

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_Relation_hpp */

