/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_Part_hpp
#define stk_mesh_Part_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <string>
#include <vector>

#include <stk_util/util/CSet.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/PartRelation.hpp>
#include <stk_mesh/baseImpl/PartImpl.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

namespace impl {
  class PartRepository;
} // namespace impl

/** \addtogroup stk_mesh_module
 *  \{
 */


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
 *  \ref stk::mesh::PartRepository "Part manager".
 */
class Part {
public:

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns the PartRepository which created this part.
   */
  MetaData & mesh_meta_data() const { return m_partImpl.mesh_meta_data(); }

  /** \brief  The primary entity type for this part.
   *
   *   For example, the primary purpose of an Element part
   *   is to define a collection of elements.  However, the
   *   nodes of those elements are also members of an element part.
   *   Return InvalidEntityRank if no primary entity type.
   */
  unsigned primary_entity_rank() const { return m_partImpl.primary_entity_rank(); }

  /** \brief  Application-defined text name of this part */
  const std::string & name() const { return m_partImpl.name(); }

  /** \brief  Internally generated ordinal of this part that is unique
   *          within the owning \ref stk::mesh::MetaData "meta data manager".
   */
  unsigned mesh_meta_data_ordinal() const { return m_partImpl.mesh_meta_data_ordinal(); }

  /** \brief  Parts that are supersets of this part. */
  const PartVector & supersets() const { return m_partImpl.supersets(); }

  /** \brief  Parts that are subsets of this part. */
  const PartVector & subsets() const { return m_partImpl.subsets(); }

  /** \brief  Parts for which this part is defined as the intersection.  */
  const PartVector & intersection_of() const { return m_partImpl.intersection_of(); }

  /** \brief  PartRelations for which this part is a member, root or target */
  const std::vector<PartRelation> & relations() const { return m_partImpl.relations(); }

  /** \brief  Equality comparison */
  bool operator == ( const Part & rhs ) const { return this == & rhs ; }

  /** \brief  Inequality comparison */
  bool operator != ( const Part & rhs ) const { return this != & rhs ; }

  /** \brief  Query attribute that has been attached to this part */
  template<class A>
  const A * attribute() const { return m_partImpl.attribute<A>(); }

private:

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns the PartRepository which created this part.
   */
  MetaData & meta_data() const { return m_partImpl.mesh_meta_data(); }


  impl::PartImpl m_partImpl;

  /* \brief  A part is owned by a PartRepository, as such only the owning
   *         PartRepository can create, delete, or modify a part.
   *         The owner-modifies rule is enforced by the implementation being
   *         a private data object on the Part and the PartRepository is a
   *         friend.
   */
  friend class ::stk::mesh::impl::PartRepository ;
  friend class ::stk::mesh::MetaData ;

#ifndef DOXYGEN_COMPILE

  /** Construct a subset part within a given mesh.
   *  Is used internally by the two 'declare_part' methods on PartRepository.
   */
  Part( MetaData * arg_meta_data , const std::string & arg_name, EntityRank arg_rank, size_t arg_ordinal)
    : m_partImpl(arg_meta_data,arg_name,arg_rank,arg_ordinal)
  { }

  ~Part() {}
  Part();
  Part( const Part & );
  Part & operator = ( const Part & );

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

