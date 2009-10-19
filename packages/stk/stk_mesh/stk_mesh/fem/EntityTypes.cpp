
#include <stk_mesh/fem/EntityTypes.hpp>

namespace stk {
namespace mesh {

const std::vector<std::string> & fem_entity_type_names()
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

