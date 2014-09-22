// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <string.h>                     // for NULL
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_util/util/string_case_compare.hpp>  // for not_equal_case
#include <string>                       // for operator<<, string
#include <typeinfo>                     // for type_info
#include <vector>                       // for vector, etc
#include "stk_mesh/base/PropertyBase.hpp"  // for Property
#include "stk_mesh/base/Types.hpp"      // for PropertyBase
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowErrorMsgIf


namespace stk {
namespace mesh {

//----------------------------------------------------------------------

Property<void>::~Property() {}

PropertyBase *
MetaData::get_property_base( const std::string    & name ,
                             const std::type_info & type ,
                                   unsigned         size ) const
{
  PropertyBase * p = NULL ;
  {
    std::vector< PropertyBase * >::const_iterator i ;
    for ( i = m_properties.begin() ;
          i != m_properties.end() && not_equal_case( (*i)->name() , name ) ;
          ++i );

    if ( i != m_properties.end() ) {
      const bool error_type =         ( (*i)->m_type != type );
      const bool error_size = size && ( (*i)->m_size != size );

      ThrowErrorMsgIf( error_type,
                       "For property name " << name << ": " <<
                       " actual_type(" << (*i)->m_type.name() <<
                       ") != request_type(" << type.name() << ")");

      ThrowErrorMsgIf( error_size,
                       "For property name " << name << ": " <<
                       " actual_size(" << (*i)->m_size <<
                       ") != request_size(" << size << ")") ;
      p = *i;
    }
  }
  return p ;
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

