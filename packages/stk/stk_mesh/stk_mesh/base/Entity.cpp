/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stddef.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

Entity::Entity( const EntityKey & arg_key )
  : m_entityImpl( arg_key )
{}


Entity::~Entity()
{}

std::string print_entity_key(const Entity& entity)
{
  return print_entity_key(entity.bucket().mesh().mesh_meta_data(),
                          entity.key());
}

std::string print_entity_key(const Entity* entity)
{
  if (entity == NULL) {
    return "NULL ENTITY";
  }
  else {
    return print_entity_key(*entity);
  }
}

//
//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

