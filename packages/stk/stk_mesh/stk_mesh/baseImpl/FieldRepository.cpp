// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_mesh/baseImpl/FieldRepository.hpp>
#include <cstring>                      // for NULL, strlen
#include <iosfwd>                       // for ostringstream
#include <sstream>                      // for operator<<, basic_ostream, etc
#include <stk_util/util/string_case_compare.hpp>  // for equal_case
#include "Shards_Array.hpp"             // for ArrayDimTag
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/FieldState.hpp"  // for ::MaximumFieldStates, etc
#include "stk_mesh/baseImpl/FieldBaseImpl.hpp"  // for FieldBaseImpl
#include "stk_util/util/ReportHandler.hpp"  // for ThrowErrorMsgIf

namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace impl {

namespace {

std::string print_field_type(const DataTraits                  & arg_traits ,
                             unsigned                            arg_rank ,
                             const shards::ArrayDimTag * const * arg_tags )
{
  ThrowRequireMsg(arg_rank < 8, "Invalid field rank: " << arg_rank);

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
void verify_field_type( const FieldBase                   & arg_field ,
                        const DataTraits                  & arg_traits ,
                        unsigned                            arg_rank ,
                        const shards::ArrayDimTag * const * arg_dim_tags ,
                        unsigned                    arg_num_states )
{

  const bool ok_traits = arg_traits.is_void
                      || & arg_traits == & arg_field.data_traits();

  const bool ok_number_states =
    ! arg_num_states || arg_num_states == arg_field.number_of_states();

  bool ok_dimension = ! arg_rank || arg_rank     == arg_field.field_array_rank() ||
                                    arg_rank + 1 == arg_field.field_array_rank() ||
                                    arg_rank - 1 == arg_field.field_array_rank() ;

  const unsigned check_rank = arg_rank < arg_field.field_array_rank() ?
                              arg_rank : arg_field.field_array_rank() ;

  for ( unsigned i = 0 ; i < check_rank && ok_dimension ; ++i ) {
    ok_dimension = arg_dim_tags[i] == arg_field.dimension_tags()[i] ;
  }

  ThrowErrorMsgIf( ! ok_traits || ! ok_number_states || ! ok_dimension,
                   " verify_field_type FAILED: Existing field = " <<
                   print_field_type( arg_field.data_traits() ,
                                     arg_field.field_array_rank() ,
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
  stk::topology::rank_t               arg_entity_rank ,
  const std::string                 & arg_name ,
  const DataTraits                  & arg_traits ,
  unsigned                            arg_array_rank ,
  const shards::ArrayDimTag * const * arg_dim_tags ,
  unsigned                            arg_num_states ) const
{
  for ( std::vector<FieldBase*>::const_iterator
        j =  m_rankedFields[arg_entity_rank].begin();
        j != m_rankedFields[arg_entity_rank].end(); ++j ) {
    if ( equal_case( (*j)->name() , arg_name ) ) {

      FieldBase* f = *j ;

      verify_field_type( *f , arg_traits , arg_array_rank , arg_dim_tags , arg_num_states );

      return f;
    }
  }
  return NULL;
}

FieldBase * FieldRepository::declare_field(
  const std::string                 & arg_name ,
  stk::topology::rank_t               arg_entity_rank ,
  const DataTraits                  & arg_traits ,
  unsigned                            arg_rank ,
  const shards::ArrayDimTag * const * arg_dim_tags ,
  unsigned                            arg_num_states ,
  MetaData                          * arg_meta_data )
{
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

  FieldBase * f[ MaximumFieldStates ] = {nullptr};

  f[0] = get_field(
                arg_entity_rank ,
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
          arg_entity_rank,
          m_fields.size() ,
          field_names[i] ,
          arg_traits ,
          arg_rank,
          arg_dim_tags,
          arg_num_states ,
          static_cast<FieldState>(i)
          );

      add_field( f[i] );
    }

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {
      f[i]->m_impl.set_field_states( f );
    }
  }

  return f[0] ;
}

void FieldRepository::verify_and_clean_restrictions(const Part& superset, const Part& subset)
{
  for ( FieldVector::iterator f = m_fields.begin() ; f != m_fields.end() ; ++f ) {
    (*f)->m_impl.verify_and_clean_restrictions( superset, subset );
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
