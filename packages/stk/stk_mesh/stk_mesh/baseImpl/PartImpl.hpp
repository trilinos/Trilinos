/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_PartImpl_hpp
#define stk_mesh_PartImpl_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <string>
#include <vector>

#include <stk_util/util/CSet.hpp>
#include <stk_mesh/base/Types.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

class PartImpl {
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
   *   Return InvalidEntityRank if no primary entity type.
   */
  unsigned primary_entity_rank() const { return m_entity_rank ; }

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
  bool operator == ( const PartImpl & rhs ) const { return this == & rhs ; }

  /** \brief  Inequality comparison */
  bool operator != ( const PartImpl & rhs ) const { return this != & rhs ; }

  /** \brief  Query attribute that has been attached to this part */
  template<class A>
  const A * attribute() const { return m_attribute.template get<A>(); }

  explicit PartImpl( MetaData * );

  void add_part_to_subset( Part & part);
  void add_part_to_superset( Part & part );
  void add_relation( PartRelation relation );
  void set_intersection_of( const PartVector & );

  template<class T>
  const T * declare_attribute_with_delete( const T *);
  template<class T>
  const T * declare_attribute_no_delete( const T *);
  template<class T>
  bool remove_attribute( const T *);

  /** Construct a subset part within a given mesh.
   *  Is used internally by the two 'declare_part' methods.
   */
  PartImpl( MetaData * meta, const std::string & name,
            EntityRank rank, size_t ordinal);

  void set_primary_entity_rank( EntityRank entity_rank );

private:

#ifndef DOXYGEN_COMPILE

  PartImpl();
  PartImpl( const PartImpl & );
  PartImpl & operator = ( const PartImpl & );

  const std::string         m_name ;
  CSet                      m_attribute ;
  PartVector                m_subsets ;
  PartVector                m_supersets ;
  PartVector                m_intersect ;
  std::vector<PartRelation> m_relations ;
  MetaData          * const m_mesh_meta_data ;
  const unsigned            m_universe_ordinal ;
  EntityRank                m_entity_rank ;

#endif /* DOXYGEN_COMPILE */

};

template<class T>
inline
const T *
PartImpl::declare_attribute_with_delete( const T * a )
{
  return m_attribute.template insert_with_delete<T>( a );
}

template<class T>
inline
const T *
PartImpl::declare_attribute_no_delete( const T * a )
{
  return m_attribute.template insert_no_delete<T>( a );
}

template<class T>
inline
bool 
PartImpl::remove_attribute( const T * a )
{
  return m_attribute.template remove<T>( a );
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

