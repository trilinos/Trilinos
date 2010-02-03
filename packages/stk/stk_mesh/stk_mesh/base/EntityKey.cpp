
#include <sstream>
#include <stdexcept>

#include <stk_util/util/StaticAssert.hpp>
#include <stk_mesh/base/EntityKey.hpp>

namespace stk {
namespace mesh {

EntityKey::EntityKey( unsigned entity_rank ,
                      EntityKey::raw_key_type entity_id )
  : key( ( raw_key_type(entity_rank) << id_digits ) | entity_id )
{
  enum { OK = StaticAssert< sizeof(EntityKey) ==
                            sizeof(EntityKey::raw_key_type) >::OK };

  const bool error_rank = ( rank() != entity_rank ) ;
  const bool error_id   = ( id()   != entity_id );

  if ( error_rank || error_id ) {
    std::ostringstream msg ;
    msg << "stk::mesh::EntityKey::EntityKey( "
        << entity_rank << " , " << entity_id << " ) FAILED" ;
    if ( error_rank ) {
      msg << " : RANK OUT OF RANGE" ;
    }
    if ( error_id ) {
      msg << " : ID OUT OF RANGE" ;
    }
    throw std::runtime_error( msg.str() );
  }
}


}
}

