
#ifndef stk_mesh_Relation_hpp
#define stk_mesh_Relation_hpp

#include <iosfwd>
#include <limits>

#include <stk_mesh/base/EntityKey.hpp>

namespace stk {
namespace mesh {

typedef uint8_t RelationKind;
typedef uintptr_t relation_attr_type ;


enum { type_digits = 4 ,
       type_mask = (~unsigned(0)) >> ( std::numeric_limits<unsigned>::digits - type_digits ),
       attr_kind_digits = 4 ,
       attr_digits      = std::numeric_limits<relation_attr_type>::digits ,
       attr_kind_shift  = attr_digits - attr_kind_digits ,
       attr_type_shift  = attr_kind_shift - type_digits ,
       attr_id_digits   = attr_type_shift < std::numeric_limits<int>::digits ?
       attr_type_shift : std::numeric_limits<int>::digits ,
       attr_kind_mask   = (~unsigned(0)) >>
       ( std::numeric_limits<unsigned>::digits - attr_kind_digits ),
       attr_id_mask     = (~unsigned(0)) >>
       ( std::numeric_limits<unsigned>::digits - attr_id_digits )
};

inline
unsigned relation_identifier( relation_attr_type attr )
{ return attr_id_mask & attr ; }

inline
unsigned relation_kind( relation_attr_type attr )
{ return attr >> attr_kind_shift ; }

inline
unsigned relation_entity_type( relation_attr_type attr )
{ return type_mask & ( attr >> attr_type_shift ); }

inline
relation_attr_type
relation_attr( unsigned entity_type ,
               unsigned identifier ,
               unsigned relation_kind )
{
  enum { kind_digits = attr_kind_digits ,
         kind_mask   = attr_kind_mask ,
         type_mask   = type_mask };

  relation_attr_type
    attr( ( ( relation_kind & kind_mask ) << kind_digits ) |
            ( entity_type   & type_mask ) );
  attr <<= attr_type_shift ;
  attr |=  attr_id_mask & identifier ;
  return attr ;
}


//----------------------------------------------------------------------
// Not quite so trivial, perhaps should not be inlined.

// inline
// entity_key_type entity_key( unsigned type , entity_id_type id )
// {
//   entity_key_type key( type );
//   key <<= EntityKey::type_shift ;
//   key |=  EntityKey::check_bit ;
//   key |=  EntityKey::id_mask & static_cast<entity_key_type>(id);
//   return key ;
// }

// inline
// relation_attr_type
// relation_attr( unsigned entity_type ,
//                unsigned identifier ,
//                unsigned relation_kind )
// {
//   enum { kind_digits = EntityKey::attr_kind_digits ,
//          kind_mask   = EntityKey::attr_kind_mask ,
//          type_mask   = EntityKey::type_mask };

//   relation_attr_type
//     attr( ( ( relation_kind & kind_mask ) << kind_digits ) |
//             ( entity_type   & type_mask ) );
//   attr <<= EntityKey::attr_type_shift ;
//   attr |=  EntityKey::attr_id_mask & identifier ;
//   return attr ;
// }

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
    if(this != &r) {
      m_attr = r.m_attr ;
      m_entity = r.m_entity ;
    }
    return *this ;
  }

  /** \brief  Construct a relation from an encoded relation attribute
   *          and a referenced entity.
   */
  Relation( relation_attr_type attr , Entity & entity );

  /** \brief  Construct a relation from a referenced entity,
   *          local identifier, kind, and converse flag.
   */
  Relation( Entity & entity ,
            unsigned identifier ,
            unsigned kind = 0 );

  /** \brief  The encoded relation attribute */
  relation_attr_type attribute() const { return m_attr ; }

  /** \brief  The kind of relation */
  unsigned   kind() const { return relation_kind( m_attr ); }

  /** \brief  The type of the referenced entity */
  unsigned entity_type() const { return relation_entity_type( m_attr ); }

  /** \brief  The local relation identifier */
  unsigned   identifier() const { return relation_identifier( m_attr ); }

#if 0
  /** \brief  If the relation is from the owning to the referenced entity */
  bool       forward() const ;

  /** \brief  If the relation is from the referenced to the owning entity */
  bool       converse() const ;
#endif

  /** \brief  The referenced entity */
  Entity * entity() const { return m_entity ; }

  /** \brief  Equality operator */
  bool operator == ( const Relation & r ) const
    { return m_attr == r.m_attr && m_entity == r.m_entity ; }

  /** \brief  Inequality operator */
  bool operator != ( const Relation & r ) const
    { return m_attr != r.m_attr || m_entity != r.m_entity ; }

  /** \brief  Ordering operator */
  bool operator < ( const Relation & r ) const ;

private:
  relation_attr_type m_attr ;
  Entity           * m_entity ;
};

//----------------------------------------------------------------------
/** \brief  Query if a member entity of the given entity type
 *          has an induced membership.
 */
bool membership_is_induced( const Part & part , unsigned entity_type );

/** \brief  Induce entities' part membership based upon relationships
 *          between entities. Insert the result into 'induced_parts'.
 */
void induced_part_membership( Part & part ,
                              unsigned entity_type_from ,
                              unsigned entity_type_to ,
                              unsigned relation_identifier ,
                              unsigned relation_kind ,
                              PartVector & induced_parts );

/** \brief  Induce entities' part membership based upon relationships
 *          between entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership( const Entity     & entity_from ,
                              const PartVector & omit ,
                                    unsigned     entity_type_to ,
                                    unsigned     relation_identifier ,
                                    unsigned     relation_kind ,
                                    PartVector & entity_to_parts );

/** \brief  Induce an entity's part membership based upon relationships
 *          from other entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership( const Entity     & entity ,
                              const PartVector & omit ,
                                    PartVector & induced );

//----------------------------------------------------------------------

/** \brief  Decode and print the relation attribute */
std::ostream &
print_relation( std::ostream & , relation_attr_type );

/** \brief  Decode and print the relation attribute and referenced entity key */
std::ostream &
print_relation( std::ostream & , const MetaData & ,
                relation_attr_type , EntityKey );

/** \brief  Print the relation attributes and referenced entity's key */
std::ostream & operator << ( std::ostream & , const Relation & );

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_Relation_hpp */

