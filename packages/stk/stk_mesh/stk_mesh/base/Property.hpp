/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_Property_hpp
#define stk_mesh_Property_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <string>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/PropertyBase.hpp>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */

/** \brief  Print the text name for a field. */
std::ostream & operator << ( std::ostream & , const PropertyBase & );

/** \brief  Print field and field restrictions on new lines. */
std::ostream & print( std::ostream & ,
                      const char * const , const PropertyBase & );

/** \brief  Query pointer to property data for a given part */
template< typename property_type >
inline
const typename property_type::data_type *
property_data( const property_type & prop , const Part & part )
{
  const PropertyBase & prop_base = dynamic_cast< const PropertyBase & > (prop);
  const MetaData * prop_meta_data = & MetaData::get(prop_base);
  const MetaData * part_meta_data = & MetaData::get(part);
  ThrowRequireMsg( prop_meta_data == part_meta_data,
                   "MetaData mismatch between property and part" );
  return prop.data( part.mesh_meta_data_ordinal() );
}

/** \brief  Query pointer to property data for a given part */
template< typename property_type >
inline
typename property_type::data_type *
property_data( property_type & prop , const Part & part )
{
  const PropertyBase & prop_base = dynamic_cast< const PropertyBase & > (prop);
  const MetaData * prop_meta_data = & MetaData::get(prop_base);
  const MetaData * part_meta_data = & MetaData::get(part);
  ThrowRequireMsg( prop_meta_data == part_meta_data,
                   "MetaData mismatch between property and part" );
  return prop.data( part.mesh_meta_data_ordinal() );
}

/** \} */

} // namespace mesh
} // namespace stk

#endif /* stk_mesh_Property_hpp */

