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

  struct FieldMetaData
  {
    int m_size;                                           // num bytes per entity, 0 means bucket does not have this field
    const FieldRestriction::size_type * m_stride;
    unsigned char* m_data;
  };
#ifdef __IBMCPP__
  // The IBM compiler is easily confused by complex template types...
  typedef std::vector<FieldMetaData>                                     FieldMetaDataVector;
  typedef std::vector<FieldMetaDataVector>                               FieldMetaDataVectorVector;
#else
  typedef std::vector<FieldMetaData, tracking_allocator<FieldMetaData, FieldDataTag> >             FieldMetaDataVector;
  typedef std::vector<FieldMetaDataVector, tracking_allocator<FieldMetaDataVector, FieldDataTag> > FieldMetaDataVectorVector;
#endif

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
  unsigned field_array_rank() const { return m_impl.field_array_rank(); }

  stk::topology::rank_t entity_rank() const { return m_impl.entity_rank(); }

  /** \brief  Multi-dimensional
   *          \ref shards::ArrayDimTag "array dimension tags"
   *          of this field.
   */
  const shards::ArrayDimTag * const * dimension_tags() const
  { return m_impl.dimension_tags() ; }

  /** \brief  Maximum field data allocation size declared for this
   *          field for the given entity rank.
   */
  unsigned max_size( EntityRank ent_rank) const {
    return m_impl.max_size( ent_rank );
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

  //----------------------------------------

  FieldBase * field_state( FieldState fstate) const {
    return m_impl.field_state(fstate);
  }

  const void* get_initial_value() const { return m_impl.get_initial_value(); }

  void* get_initial_value() { return m_impl.get_initial_value(); }

  unsigned get_initial_value_num_bytes() const { return m_impl.get_initial_value_num_bytes(); }

  virtual ~FieldBase() {}

  virtual std::ostream& print_data(std::ostream& out, void* data, unsigned size_per_entity) const { return out; }

  stk::mesh::BulkData& get_mesh() const
  { return *m_mesh; }


  //
  //  NKC, don't like that this is non-const.  Bulk data needs to modify this, but no other classes should, maybe 
  //  this should be private and have BulkData as friend?  Or maybe restructure in a different way so that bucket
  //  changes are cleanly propogated to field.
  //
  inline FieldMetaDataVectorVector& get_meta_data_for_field() {
    return m_field_meta_data;
  }
  inline const FieldMetaDataVectorVector& get_meta_data_for_field() const{
    return m_field_meta_data;
  }


private:

  //  Associate this field with a bulk data.
  //    Note, a field can be assocaited with one and only one bulk data object
  void set_mesh(stk::mesh::BulkData* bulk);

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns this field
   */
  MetaData & meta_data() const { return m_impl.meta_data(); }

  friend class ::stk::mesh::MetaData ;
  friend class ::stk::mesh::impl::FieldRepository;
  friend class ::stk::mesh::impl::FieldBaseImpl ;

  /** \brief  Allow the unit test driver access */
  friend class ::stk::mesh::UnitTestFieldImpl ;

  FieldMetaDataVectorVector m_field_meta_data;

protected:
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
    : m_mesh(NULL),
      m_impl(
        arg_mesh_meta_data,
        stk::topology::INVALID_RANK,
        arg_ordinal,
        arg_name,
        arg_traits,
        arg_rank,
        arg_dim_tags,
        arg_number_of_states,
        arg_this_state
        )
  {}

  FieldBase(
      MetaData                   * arg_mesh_meta_data ,
      stk::topology::rank_t        arg_entity_rank ,
      unsigned                     arg_ordinal ,
      const std::string          & arg_name ,
      const DataTraits           & arg_traits ,
      unsigned                     arg_rank,
      const shards::ArrayDimTag  * const * arg_dim_tags,
      unsigned                     arg_number_of_states ,
      FieldState                   arg_this_state
      )
    : m_mesh(NULL),
      m_impl(
        arg_mesh_meta_data,
        arg_entity_rank,
        arg_ordinal,
        arg_name,
        arg_traits,
        arg_rank,
        arg_dim_tags,
        arg_number_of_states,
        arg_this_state
        )
  {}

private:

  stk::mesh::BulkData* m_mesh;
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
                      const char * const ,
                      const FieldBase & );


std::ostream & print_restrictions( std::ostream & ,
                                   const char * const ,
                                   const FieldBase & );


} //namespace mesh
} //namespace stk

#endif //stk_mesh_base_FieldBase_hpp
