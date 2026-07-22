// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/FieldState.hpp"  // for ::MaximumFieldStates, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowErrorMsgIf

namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace impl {

namespace {

std::string print_field_type(const DataTraits & arg_traits)
{
  std::ostringstream oss;
  oss << "Field<" << arg_traits.name << ">";
  return oss.str();
}

} //unamed namespace

// Check for compatibility:
// 1) Scalar type must match
// 2) Number of states must match
// 3) Dimension must be different by at most one rank,
//    where the tags match for the smaller rank.
void
FieldRepository::verify_field_type(const FieldBase  & arg_field,
                                   const DataTraits & arg_traits,
                                   unsigned           arg_num_states) const
{

  const bool ok_traits = arg_traits.is_void || &arg_traits == &arg_field.data_traits();

  const bool ok_number_states = not arg_num_states || arg_num_states == arg_field.number_of_states();

  STK_ThrowErrorMsgIf(not ok_traits || not ok_number_states,
                      " verify_field_type FAILED: Existing field = " <<
                      print_field_type(arg_field.data_traits()) << "[ name = \"" << arg_field.name() <<
                      "\" , #states = " << arg_field.number_of_states() << " ]" <<
                      " Expected field info = " << print_field_type(arg_traits) <<
                      "[ #states = " << arg_num_states << " ]");
}

//----------------------------------------------------------------------
FieldBase * FieldRepository::get_field(stk::topology::rank_t   arg_entity_rank,
                                       const std::string     & arg_name,
                                       const DataTraits      & arg_traits,
                                       unsigned                arg_num_states) const
{
  for (FieldBase * field : m_rankedFields[arg_entity_rank]) {
    if (equal_case(field->name(), arg_name)) {
      verify_field_type(*field, arg_traits, arg_num_states);
      return field;
    }
  }

  return nullptr;
}

void FieldRepository::verify_and_clean_restrictions(const Part& superset, const Part& subset)
{
  stk::mesh::EntityRank partRank = subset.primary_entity_rank();
  for ( FieldVector::iterator f = m_fields.begin() ; f != m_fields.end() ; ++f ) {
    if (partRank == stk::topology::INVALID_RANK || partRank == (*f)->entity_rank())
    {
      (*f)->verify_and_clean_restrictions( superset, subset );
    }
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
