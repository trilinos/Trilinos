/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_base_FieldBase_hpp
#define stk_mesh_base_FieldBase_hpp

#include <iosfwd>
#include <string>
#include <vector>

#include <Shards_Array.hpp>

#include <stk_util/util/SimpleArrayOps.hpp>
#include <stk_util/util/CSet.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/base/FieldState.hpp>

#include <stk_mesh/base/FieldRestriction.hpp>
#include <stk_mesh/baseImpl/FieldBaseImpl.hpp>

namespace stk {
namespace mesh {

class UnitTestFieldImpl;

namespace impl {

class FieldRepository;

}

//----------------------------------------------------------------------
/** \ingroup stk_stk_mesh_module
 *  \brief  Field base class with an anonymous data type and
 *          anonymous multi-dimension.
 *
 * This class is the base class for all Fields. It defines the member-function
 * API for Field objects. Note that the entire Field API is much broader than
 * what's defined in this class. See Field.hpp for more.
 */
// Implementation details:
//   Simply wraps a FieldBaseImpl object which is kept as a member, all calls are
//   inlined and passed through to the impl object. This design is analogous
//   to the "P-impl" pattern.
class FieldBase
{
  public:
  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns this field
   */
  MetaData & mesh_meta_data() const { return m_impl.meta_data(); }
  /** \brief  Internally generated ordinal of this field that is unique
   *          within the owning \ref stk::mesh::MetaData "meta data manager".
   */
  unsigned mesh_meta_data_ordinal() const { return m_impl.ordinal(); }

  /** \brief  Application-defined text name of this field */
  const std::string & name() const { return m_impl.name() ; }

  /** \brief  Query if the type is Type */
  template<class Type> bool type_is() const
  { return m_impl.type_is<Type>(); }

  /** \brief  Return the \ref stk::mesh::DataTraits "data traits"
   *          for this field's type
   */
  const DataTraits & data_traits() const { return m_impl.data_traits() ; }

  /** \brief  Number of states of this field */
  unsigned number_of_states() const { return m_impl.number_of_states() ; }

  /** \brief  FieldState of this field */
  FieldState state() const { return m_impl.state() ; }

  /** \brief  Multi-dimensional array rank of this field,
   *          which is zero for a scalar field.
   */
  unsigned rank() const { return m_impl.rank(); }

  /** \brief  Multi-dimensional
   *          \ref shards::ArrayDimTag "array dimension tags"
   *          of this field.
   */
  const shards::ArrayDimTag * const * dimension_tags() const
  { return m_impl.dimension_tags() ; }

  /** \brief  Maximum field data allocation size declared for this
   *          field for the given entity rank.
   */
  unsigned max_size( EntityRank entity_rank) const {
    return m_impl.max_size( entity_rank );
  }

  //----------------------------------------

  /** \brief  Query attribute that has been attached to this field */
  template<class A>
  const A * attribute() const { return m_impl.attribute<A>(); }

  typedef FieldRestriction Restriction;

  /** \brief  A fields' restrictions are maintained in a std::vector */
  typedef FieldRestrictionVector RestrictionVector;

  /** \brief  Vector of field restriction which is volatile until the owning
   *          \ref stk::mesh::MetaData "meta data manager" is committed.
   */
  const RestrictionVector &restrictions() const {
    return m_impl.restrictions();
  }

  /** \brief  Query a field restriction, result is volatile until the owning
   *          \ref stk::mesh::MetaData "meta data manager" is committed.
   */
  const Restriction & restriction( unsigned entity_rank , const Part & part ) const {
    return m_impl.restriction( entity_rank, part);
  }

  //----------------------------------------

  FieldBase * field_state( FieldState fstate) const {
    return m_impl.field_state(fstate);
  }

private:

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns this field
   */
  MetaData & meta_data() const { return m_impl.meta_data(); }

  friend class ::stk::mesh::MetaData ;
  friend class ::stk::mesh::impl::FieldRepository ;
  friend class ::stk::mesh::impl::FieldBaseImpl ;

  /** \brief  Allow the unit test driver access */
  friend class ::stk::mesh::UnitTestFieldImpl ;

  FieldBase(
      MetaData                   * arg_mesh_meta_data ,
      unsigned                     arg_ordinal ,
      const std::string          & arg_name ,
      const DataTraits           & arg_traits ,
      unsigned                     arg_rank,
      const shards::ArrayDimTag  * const * arg_dim_tags,
      unsigned                     arg_number_of_states ,
      FieldState                   arg_this_state
      )
    : m_impl(
        arg_mesh_meta_data,
        arg_ordinal,
        arg_name,
        arg_traits,
        arg_rank,
        arg_dim_tags,
        arg_number_of_states,
        arg_this_state
        )
  {}

  // WORKAROUND 5/19/2010 [DGB]: intel 10.? and pgi do not link if this is made virtual
  //  virtual ~FieldBase();
  ~FieldBase();

  impl::FieldBaseImpl  m_impl;

  //the following functions are declared but not defined
  FieldBase();
  FieldBase( const FieldBase & );
  FieldBase & operator = ( const FieldBase & );
};

/** \brief  Print the field type, text name, and number of states. */
std::ostream & operator << ( std::ostream & , const FieldBase & );

/** \brief  Print field and field restrictions on new lines. */
std::ostream & print( std::ostream & ,
                      const char * const , const FieldBase & );

} //namespace mesh
} //namespace stk

#endif //stk_mesh_base_FieldBase_hpp
