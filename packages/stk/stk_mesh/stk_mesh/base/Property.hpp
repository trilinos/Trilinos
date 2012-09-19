/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
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

