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

#include <stk_mesh/base/FieldBase.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <vector>                       // for vector, etc
#include "Shards_Array.hpp"             // for ArrayDimTag
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/FieldRestriction.hpp"  // for FieldRestriction
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequireMsg
namespace stk { namespace mesh { class BulkData; } }


namespace stk {
namespace mesh {

std::ostream & operator << ( std::ostream & s , const FieldBase & field )
{
  s << "Field<" ;
  s << field.data_traits().name ;
  for ( unsigned i = 0 ; i < field.field_array_rank() ; ++i ) {
    s << "," << field.dimension_tags()[i]->name();
  }
  s << ">" ;

  s << "[ name: \"" ;
  s << field.name() ;
  s << "\" , #states: " ;
  s << field.number_of_states();
  s << " ]" ;
  return s ;
}

std::ostream & print( std::ostream & s ,
                      const char * const b ,
                      const FieldBase & field )
{
  s << b << field << std::endl;
  std::string indent = b;
  indent += "  ";
  print_restrictions(s, indent.c_str(), field);
  return s ;
}

std::ostream & print_restrictions(std::ostream & s ,
                                  const char * const b ,
                                  const FieldBase & field )
{
  const std::vector<FieldBase::Restriction> & rMap = field.restrictions();

  for ( std::vector<FieldBase::Restriction>::const_iterator
        i = rMap.begin() ; i != rMap.end() ; ++i ) {
    s << b;
    i->print( s, i->selector(), field.field_array_rank() );
    s << std::endl;
  }
  return s;
}

void FieldBase::set_mesh(stk::mesh::BulkData* bulk)
{
  if (m_mesh == NULL || bulk == NULL) {
    m_mesh = bulk;
  }
  else {
    ThrowRequireMsg(bulk == m_mesh, "Internal Error: Trying to use field " << name() << " on more than one bulk data");
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

