/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/fem/EntityRanks.hpp>

namespace stk {
namespace mesh {

// DEPRECATED: 09/15/10 FEM TopologicalMetaData refactor
const std::vector<std::string> & fem_entity_rank_names()
{
  static std::vector<std::string> names ;
  if ( names.empty() ) {
    names.resize( 6 );
    names[0].assign( "NODE" );
    names[1].assign( "EDGE" );
    names[2].assign( "FACE" );
    names[3].assign( "ELEMENT" );
    names[4].assign( "PARTICLE" );
    names[5].assign( "CONSTRAINT" );
  }
  return names ;
}

}//namespace mesh
}//namespace stk

