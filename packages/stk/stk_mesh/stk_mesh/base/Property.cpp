/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <string.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_util/util/string_case_compare.hpp>
#include <stk_mesh/base/MetaData.hpp>

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

      if ( error_type || error_size ) {
        std::ostringstream msg ;
        msg << "stk::mesh::MetaData::get_property( " << name << " ) FAILED:" ;
        if ( error_type ) {
          msg << " actual_type(" << (*i)->m_type.name();
          msg << ") != request_type(" << type.name() << ")" ;
        }
        if ( error_type ) {
          msg << " actual_size(" << (*i)->m_size ;
          msg << ") != request_size(" << size << ")" ;
        }
        throw std::runtime_error( msg.str() );
      }
      p = *i;
    }
  }
  return p ;
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

