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
#include <stk_mesh/base/Bucket.hpp>

#include <boost/mpl/assert.hpp>

#ifdef SIERRA_MIGRATION
namespace {
static const stk::mesh::RelationVector dummy_vector;
}

namespace sierra {
namespace Fmwk {

  //const unsigned int INVALID_LOCAL_ID = stk::mesh::GetInvalidLocalId();
const stk::mesh::RelationIterator INVALID_RELATION_ITR = dummy_vector.end(); // Some STL implementation use POD for iterators

} // namespace Fmwk
} // namespace Sierra
#endif

namespace stk {
namespace mesh {

std::ostream & operator << ( std::ostream &os , const Entity &entity )
{
  os << entity.m_value;
  return os;
}

// TODO - Activate once we move to intel-12.1
//BOOST_MPL_ASSERT(( boost::is_pod<Entity> ));

//----------------------------------------------------------------------

#ifdef SIERRA_MIGRATION

std::string Entity::TypeToString (Entity::ObjectTypeEnum type)
{
  if(type == NODE      ) return "NODE";
  if(type == EDGE      ) return "EDGE";
  if(type == FACE      ) return "FACE";
  if(type == ELEMENT   ) return "ELEMENT";
  if(type == CONSTRAINT) return "CONSTRANT";
  if(type == BASE_CLASS) return "BASE_CLASS";
  return "UNKNOWN";
}

// ---------------------------------------------------------------------

#endif // SIERRA_MIGRATION

BOOST_STATIC_ASSERT(( static_cast<int>(stk::topology::NODE_RANK) == static_cast<int>(Entity::NODE) ));
BOOST_STATIC_ASSERT(( static_cast<int>(stk::topology::EDGE_RANK) == static_cast<int>(Entity::EDGE) ));
BOOST_STATIC_ASSERT(( static_cast<int>(stk::topology::FACE_RANK) == static_cast<int>(Entity::FACE) ));
BOOST_STATIC_ASSERT(( static_cast<int>(stk::topology::ELEMENT_RANK) == static_cast<int>(Entity::ELEMENT) ));

} // namespace mesh
} // namespace stk
