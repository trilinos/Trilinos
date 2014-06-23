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

