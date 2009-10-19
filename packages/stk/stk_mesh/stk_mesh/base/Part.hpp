#ifndef stk_mesh_Part_hpp
#define stk_mesh_Part_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <string>
#include <vector>

#include <stk_util/util/CSet.hpp>
#include <stk_mesh/base/Types.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------
/** \ingroup stk_mesh_relations
 *  \brief  A defined entity-relationship between
 *          \ref stk::mesh::Part "parts".
 *          An internal class that should never need to be
 *          <em> directly </em> used within application code.
 *
 *  <b> If </b> an
 *     \ref stk::mesh::Entity "entity" <b> e1 </b> is a member of
 *     \ref stk::mesh::Part "part" <em> m_root </em> and
 *     there exists a
 *     \ref stk::mesh::Relation "relation"
 *     from entity <b> e1 </b> to entity <b> e2 </b> that
 *     is in the domain of the
 *     \ref stk_mesh_relations "relation stencil"
 *     <em> m_function </em>
 *  <b> then </b> entity <b> e2 </b> is a member of part <em> m_target </em>.
 */
struct PartRelation {
  /** \brief  relation domain part */
  Part * m_root ;

  /** \brief  relation range part */
  Part * m_target ;

  /** \brief  \ref stk_mesh_relations "relation stencil" */
  relation_stencil_ptr m_function ;

#ifndef DOXYGEN_COMPILE

  ~PartRelation() {}

  PartRelation() : m_root( NULL ), m_target( NULL ), m_function( NULL ) {}

  PartRelation( const PartRelation & rhs )
    : m_root( rhs.m_root ),
      m_target( rhs.m_target ),
      m_function( rhs.m_function ) {}

  PartRelation & operator = ( const PartRelation & rhs )
  {
    m_root = rhs.m_root ;
    m_target = rhs.m_target ;
    m_function = rhs.m_function ;
    return *this ;
  }

#endif /* DOXYGEN_COMPILE */

};

//----------------------------------------------------------------------
/** \brief  An application-defined subset of a problem domain.
 *
 *  An application may define parts corresponding to geometric subdomains,
 *  modeling subdomains, material subdomains, parallel distributed subdomains,
 *  entities of the same type of discretization
 *  (e.g., hexahedrons, tetrahedrons, ...),
 *  or any other application need for subsetting.
 *
 *  A Part is created, owned, and modified by a
 *  \ref stk::mesh::MetaData "meta data manager".
 */
class Part {
public:

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns this part
   */
  MetaData & mesh_meta_data() const { return *m_mesh_meta_data ; }

  /** \brief  The primary entity type for this part.
   *
   *   For example, the primary purpose of an Element part 
   *   is to define a collection of elements.  However, the
   *   nodes of those elements are also members of an element part.
   *   Return std::numeric_limits<unsigned>::max() if no primary entity type.
   */
  unsigned primary_entity_type() const { return m_entity_rank ; }

  /** \brief  Application-defined text name of this part */
  const std::string & name() const { return m_name ; }

  /** \brief  Internally generated ordinal of this part that is unique
   *          within the owning \ref stk::mesh::MetaData "meta data manager".
   */
  unsigned mesh_meta_data_ordinal() const { return m_universe_ordinal ; }

  /** \brief  Parts that are supersets of this part. */
  const PartVector & supersets() const { return m_supersets ; }

  /** \brief  Parts that are subsets of this part. */
  const PartVector & subsets() const { return m_subsets ; }

  /** \brief  Parts for which this part is defined as the intersection.  */
  const PartVector & intersection_of() const { return m_intersect ; }

  /** \brief  PartRelations for which this part is a member, root or target */
  const std::vector<PartRelation> & relations() const { return m_relations ; }

  /** \brief  Equality comparison */
  bool operator == ( const Part & rhs ) const { return this == & rhs ; }

  /** \brief  Inequality comparison */
  bool operator != ( const Part & rhs ) const { return this != & rhs ; }

  /** \brief  Query attribute that has been attached to this part */
  template<class A>
  const A * attribute() const { return m_attribute.template get<A>(); }

private:

  /* \brief  A part is owned by a MetaData, as such only the owning
   *         MetaData can create, delete, or modify a part.
   *         The owner-modifies rule is enforced by all non-const
   *         methods being private and the MetaData be a friend.
   */
  friend class ::stk::mesh::MetaData ;

  /** \brief  Allow the unit test driver access */
  friend class ::stk::mesh::UnitTestMetaData ;

  /** \brief  Construct a universal part */
  explicit Part( MetaData * );

  /** \brief  Declare a subset part on a universal part */
  Part & declare_part( const std::string & arg_name , EntityType arg_rank );

  /** \brief  Declare an intersecton part on a universal part */
  Part & declare_part( const PartVector & part_intersect );

  /** \brief  Declare a superset <-> subset relationship */
  void declare_subset( Part & subset );

#ifndef DOXYGEN_COMPILE

  /** Construct a subset part within a given mesh.
   *  Is used internally by the two 'declare_part' methods.
   */
  Part( MetaData * , const std::string & , EntityType , size_t );

  ~Part();
  Part();
  Part( const Part & );
  Part & operator = ( const Part & );

  const std::string         m_name ;
  CSet                      m_attribute ;
  PartVector                m_subsets ;
  PartVector                m_supersets ;
  PartVector                m_intersect ;
  std::vector<PartRelation> m_relations ;
  MetaData          * const m_mesh_meta_data ;
  const unsigned            m_universe_ordinal ;
  const unsigned            m_entity_rank ;

#endif /* DOXYGEN_COMPILE */

};

//----------------------------------------------------------------------
/** \brief  Ordering operator for parts. */
struct PartLess {
 
  inline bool operator()( const Part & lhs , const Part & rhs ) const
    { return lhs.mesh_meta_data_ordinal() < rhs.mesh_meta_data_ordinal(); }

  inline bool operator()( const Part & lhs , const Part * rhs ) const
    { return lhs.mesh_meta_data_ordinal() < rhs->mesh_meta_data_ordinal(); }

  inline bool operator()( const Part * lhs , const Part & rhs ) const
    { return lhs->mesh_meta_data_ordinal() < rhs.mesh_meta_data_ordinal(); }
 
  inline bool operator()( const Part * lhs , const Part * rhs ) const
    { return lhs->mesh_meta_data_ordinal() < rhs->mesh_meta_data_ordinal(); }
};

/** \brief  Order a collection of parts: invoke sort and then unique */
void order( PartVector & );

/** \brief  Insert a part into a properly ordered collection of parts.
 *          Returns true if this is a new insertion.
 */
bool insert( PartVector & , Part & );

/** \brief  Remove a part from a properly ordered collection of parts. */
void remove( PartVector & , Part & );

/** \brief  Find a part by name in a collection of parts. */
Part * find( const PartVector & , const std::string & );

/** \brief  Query containment within properly ordered PartVector */
bool contain( const PartVector & , const Part & );

/** \brief  Query containment for two properly ordered PartVector */
bool contain( const PartVector & , const PartVector & );

/** \brief  Query cardinality of intersection of two PartVectors */
size_t intersect( const PartVector & , const PartVector & );

/** \brief  Generate the intersection of two PartVectors */
size_t intersect( const PartVector & , const PartVector & , PartVector & );

/** \brief  Query if two parts intersect; i.e.,
 *          if one is a subset of the other or they share a common subset
 */
bool intersect( const Part & , const Part & );

//----------------------------------------------------------------------
/** \brief  Print a part's information including supersets, subsets, and
 *          intersection.  Each line starts with the given leader string.
 */
std::ostream & print( std::ostream & , const char * const , const Part & );

/** \} */


} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

