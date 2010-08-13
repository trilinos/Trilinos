/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <stdexcept>

#include <stk_mesh/femImpl/FiniteElementMeshImpl.hpp>

namespace stk {
namespace mesh {
namespace impl {

void verify_spatial_dimension( const char * method ,
                               unsigned spatial_dimension )
{
  if ( spatial_dimension < 1 || 3 < spatial_dimension ) {
    std::ostringstream msg ;
    msg << method << " ERRONEOUS spatial dimension = "
        << spatial_dimension ;
    throw std::logic_error( msg.str() );
  }
}

std::vector<std::string>
finite_element_mesh_entity_rank_names( unsigned spatial_dimension )
{
  static const char method[] =
    "stk::mesh::femImpl::finite_element_mesh_entity_rank_names" ;

  verify_spatial_dimension( method , spatial_dimension );

  std::vector< std::string > names ;

  names.reserve( spatial_dimension + 1 );

  names.push_back( std::string( "NODE" ) );

  if ( 1 < spatial_dimension ) { names.push_back( std::string("EDGE") ); }
  if ( 2 < spatial_dimension ) { names.push_back( std::string("FACE") ); }

  names.push_back( std::string("ELEMENT") );
  names.push_back( std::string("PATCH") );

  return names ;
}

}
}
}


