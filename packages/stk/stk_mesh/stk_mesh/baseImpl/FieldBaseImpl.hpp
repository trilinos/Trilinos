/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#ifndef stk_mesh_baseImpl_FieldBaseImpl_hpp
#define stk_mesh_baseImpl_FieldBaseImpl_hpp

#include <vector>
#include <Shards_Array.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/base/FieldState.hpp>
#include <stk_mesh/base/FieldRestriction.hpp>

#include <stk_util/util/CSet.hpp>

namespace stk {
namespace mesh {

class FieldBase;
class MetaData;

namespace impl {

class FieldBaseImpl {
public:

  MetaData & meta_data() const {
    return *m_meta_data ;
  }

  unsigned ordinal() const {
    return m_ordinal ;
  }

  const std::string & name() const {
    return m_name ;
  }

  template<class Type>
  bool type_is() const {
    return m_data_traits.type_info == typeid(Type) ;
  }

  const DataTraits & data_traits() const {
    return m_data_traits ;
  }

  unsigned number_of_states() const {
    return m_num_states ;
  }

  FieldState state() const {
    return m_this_state ;
  }

  unsigned rank() const {
    return m_rank ;
  }

  const shards::ArrayDimTag * const * dimension_tags() const {
    return m_dim_tags ;
  }

  //not inline
  unsigned max_size( unsigned entity_rank) const ;

  template<class A>
  const A * attribute() const {
    return m_attribute.template get<A>();
  }

  template<class A>
    const A * declare_attribute_no_delete(const A * a) {
      return m_attribute.template insert_no_delete<A>(a);
    }

  template<class A>
    const A * declare_attribute_with_delete(const A * a) {
      return m_attribute.template insert_with_delete<A>(a);
    }

  const FieldRestriction & restriction( unsigned entity_rank , const Part & part) const ;

  const FieldRestrictionVector &restrictions() const ;

  FieldBase * field_state(FieldState fstate) const {
    return m_field_states[fstate];
  }

  void insert_restriction( const char       * arg_method ,
                           EntityRank         arg_entity_rank ,
                           const Part       & arg_part ,
                           const unsigned   * arg_stride );

  void verify_and_clean_restrictions( const char       * arg_method ,
                                      const PartVector & arg_all_parts );



  void set_field_states( FieldBase ** field_states);

  FieldBaseImpl(
      MetaData                   * arg_mesh_meta_data ,
      unsigned                     arg_ordinal ,
      const std::string          & arg_name ,
      const DataTraits           & arg_traits ,
      unsigned                     arg_rank,
      const shards::ArrayDimTag  * const * arg_dim_tags,
      unsigned                     arg_number_of_states ,
      FieldState                   arg_this_state
      );

  ~FieldBaseImpl() {}

private:

  FieldRestrictionVector & restrictions();

  const std::string            m_name ;                    ///< Name of the field
  CSet                         m_attribute ;               ///< User's attributes
  const DataTraits           & m_data_traits ;             ///< Data type traits
  MetaData             * const m_meta_data ;               ///< Owner of this field
  const unsigned               m_ordinal ;                 ///< Ordinal in the field set
  const unsigned               m_num_states ;              ///< Number of states
  const FieldState             m_this_state ;              ///< Field state of this field
  unsigned                     m_rank ;                    ///< Number of dimensions
  FieldRestrictionVector       m_dim_map ;                 ///< Only valid on StateNone
  FieldBase                  * m_field_states[ MaximumFieldStates ];
  const shards::ArrayDimTag  * m_dim_tags[ MaximumFieldDimension ];

  //disallow copy and default constructors
  FieldBaseImpl();
  FieldBaseImpl( const FieldBase & );
  FieldBaseImpl & operator = ( const FieldBaseImpl & );
};


/** \brief  Print the field type, text name, and number of states. */
std::ostream & operator << ( std::ostream & , const FieldBaseImpl & );

/** \brief  Print field and field restrictions on new lines. */
std::ostream & print( std::ostream & ,
                      const char * const , const FieldBase & );


} // namespace impl
} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_FieldBaseIml_hpp
