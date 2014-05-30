/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>
#include <stdexcept>

#include <stk_util/util/StaticAssert.hpp>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/base/EntityKey.hpp>

namespace stk_classic {
namespace mesh {

EntityKey::EntityKey( EntityRank entity_rank ,
                      EntityKey::raw_key_type entity_id )
  : key( ( raw_key_type(entity_rank) << id_digits ) | entity_id )
{
  enum { OK = StaticAssert< sizeof(EntityKey) ==
                            sizeof(EntityKey::raw_key_type) >::OK };

  ThrowAssertMsg( rank() == entity_rank,
                  "entity_rank out of range, entity_rank= " << entity_rank << " rank() = " << rank() << " entity_id= " << entity_id << " id() = " << id() );

  ThrowAssertMsg( id() == entity_id,
                  "entity_id out of range, entity_rank= " << entity_rank << " rank() = " << rank() << " entity_id= " << entity_id << " id() = " << id() );
}


}
}

