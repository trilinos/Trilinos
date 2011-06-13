/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/base/Trace.hpp>

#include <stk_mesh/baseImpl/FieldRepository.hpp>

#include <stk_util/util/string_case_compare.hpp>

#include <cstring>
#include <sstream>
#include <stdexcept>

namespace stk {
namespace mesh {
namespace impl {

namespace {

std::string print_field_type(const DataTraits                  & arg_traits ,
                             unsigned                            arg_rank ,
                             const shards::ArrayDimTag * const * arg_tags )
{
  std::ostringstream oss;
  oss << "FieldBase<" ;
  oss << arg_traits.name ;
  for ( unsigned i = 0 ; i < arg_rank ; ++i ) {
    oss << "," << arg_tags[i]->name();
  }
  oss << ">" ;
  return oss.str();
}

// Check for compatibility:
// 1) Scalar type must match
// 2) Number of states must match
// 3) Dimension must be different by at most one rank,
//    where the tags match for the smaller rank.
void verify_field_type( const char                        * arg_method ,
                        const FieldBase                   & arg_field ,
                        const DataTraits                  & arg_traits ,
                        unsigned                            arg_rank ,
                        const shards::ArrayDimTag * const * arg_dim_tags ,
                        unsigned                    arg_num_states )
{

  const bool ok_traits = arg_traits.is_void
                      || & arg_traits == & arg_field.data_traits();

  const bool ok_number_states =
    ! arg_num_states || arg_num_states == arg_field.number_of_states();

  bool ok_dimension = ! arg_rank || arg_rank     == arg_field.rank() ||
                                    arg_rank + 1 == arg_field.rank() ||
                                    arg_rank - 1 == arg_field.rank() ;

  const unsigned check_rank = arg_rank < arg_field.rank() ?
                              arg_rank : arg_field.rank() ;

  for ( unsigned i = 0 ; i < check_rank && ok_dimension ; ++i ) {
    ok_dimension = arg_dim_tags[i] == arg_field.dimension_tags()[i] ;
  }

  ThrowErrorMsgIf( ! ok_traits || ! ok_number_states || ! ok_dimension,
                   arg_method << " FAILED: Existing field = " <<
                   print_field_type( arg_field.data_traits() ,
                                     arg_field.rank() ,
                                     arg_field.dimension_tags() ) <<
                   "[ name = \"" << arg_field.name() <<
                   "\" , #states = " << arg_field.number_of_states() << " ]" <<
                   " Expected field info = " <<
                   print_field_type( arg_traits , arg_rank , arg_dim_tags ) <<
                   "[ #states = " << arg_num_states << " ]");
}

} //unamed namespace

//----------------------------------------------------------------------

FieldBase * FieldRepository::get_field(
  const char                        * arg_method ,
  const std::string                 & arg_name ,
  const DataTraits                  & arg_traits ,
  unsigned                            arg_rank ,
  const shards::ArrayDimTag * const * arg_dim_tags ,
  unsigned                            arg_num_states ) const
{
  FieldBase * f = NULL ;

  for ( std::vector<FieldBase*>::const_iterator
        j =  m_fields.begin() ;
        j != m_fields.end() && NULL == f ; ++j ) {
    if ( equal_case( (*j)->name() , arg_name ) ) {

      f = *j ;

      verify_field_type( arg_method , *f , arg_traits ,
                         arg_rank , arg_dim_tags , arg_num_states );
    }
  }
  return f ;
}

FieldBase * FieldRepository::declare_field(
  const std::string                 & arg_name ,
  const DataTraits                  & arg_traits ,
  unsigned                            arg_rank ,
  const shards::ArrayDimTag * const * arg_dim_tags ,
  unsigned                            arg_num_states ,
  MetaData                          * arg_meta_data )
{
  TraceIf("stk::mesh::impl::FieldRepository::declare_field", LOG_FIELD);

  static const char* reserved_state_suffix[6] = {
    "_STKFS_OLD",
    "_STKFS_N",
    "_STKFS_NM1",
    "_STKFS_NM2",
    "_STKFS_NM3",
    "_STKFS_NM4"
  };

  // Check that the name does not have a reserved suffix

  for ( unsigned i = 0 ; i < 6 ; ++i ) {
    const int len_name   = arg_name.size();
    const int len_suffix = std::strlen( reserved_state_suffix[i] );
    const int offset     = len_name - len_suffix ;
    if ( 0 <= offset ) {
      const char * const name_suffix = arg_name.c_str() + offset ;
      ThrowErrorMsgIf( equal_case( name_suffix , reserved_state_suffix[i] ),
                       "For name = \"" << name_suffix <<
                       "\" CANNOT HAVE THE RESERVED STATE SUFFIX \"" <<
                       reserved_state_suffix[i] << "\"" );
    }
  }

  // Check that the field of this name has not already been declared

  FieldBase * f[ MaximumFieldStates ] ;

  f[0] = get_field(
                "FieldRepository::declare_field" ,
                arg_name ,
                arg_traits ,
                arg_rank ,
                arg_dim_tags ,
                arg_num_states
                );

  if ( NULL != f[0] ) {
    for ( unsigned i = 1 ; i < arg_num_states ; ++i ) {
      f[i] = f[0]->m_impl.field_state(static_cast<FieldState>(i));
    }
  }
  else {
    // Field does not exist then create it

    std::string field_names[ MaximumFieldStates ];

    field_names[0] = arg_name ;

    if ( 2 == arg_num_states ) {
      field_names[1] = arg_name ;
      field_names[1].append( reserved_state_suffix[0] );
    }
    else {
      for ( unsigned i = 1 ; i < arg_num_states ; ++i ) {
        field_names[i] = arg_name ;
        field_names[i].append( reserved_state_suffix[i] );
      }
    }

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {

      f[i] = new FieldBase(
          arg_meta_data ,
          m_fields.size() ,
          field_names[i] ,
          arg_traits ,
          arg_rank,
          arg_dim_tags,
          arg_num_states ,
          static_cast<FieldState>(i)
          );

      m_fields.push_back( f[i] );
    }

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {
      f[i]->m_impl.set_field_states( f );
    }
  }

  return f[0] ;
}

void FieldRepository::verify_and_clean_restrictions(
    const char       * arg_method ,
    const PartVector & arg_all_parts )
{
  TraceIf("stk::mesh::impl::FieldRepository::verify_and_clean_restrictions", LOG_FIELD);

  for ( FieldVector::iterator f = m_fields.begin() ; f != m_fields.end() ; ++f ) {
    (*f)->m_impl.verify_and_clean_restrictions( arg_method , arg_all_parts );
  }
}

FieldRepository::~FieldRepository() {
  try {
    FieldVector::iterator j = m_fields.begin();
    for ( ; j != m_fields.end() ; ++j ) { delete *j ; }
    m_fields.clear();
  } catch(...) {}
}

} // namespace impl
} // namespace mesh
} // namespace stk
