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
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  os << entity.bulk_data_id() << "[" << entity.local_offset() << "]";
#else
  os << entity.m_value;
#endif
  return os;
}

#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS

bool Entity::is_valid() const
{
  return (m_value != 0) && BulkData::get(*this).is_valid(*this);
}

EntityState Entity::state() const
{
  return BulkData::get(*this).state(*this);
}

EntityRank Entity::entity_rank() const
{
  return BulkData::get(*this).entity_rank(*this);
}

EntityId Entity::identifier() const
{
  return BulkData::get(*this).identifier(*this);
}

const EntityKey Entity::key() const
{
  return BulkData::get(*this).entity_key(*this);
}

Bucket & Entity::bucket() const
{
  return BulkData::get(*this).bucket(*this);
}

Bucket * Entity::bucket_ptr() const
{
  return BulkData::get(*this).bucket_ptr(*this);
}

Bucket::size_type Entity::bucket_ordinal() const
{
  return BulkData::get(*this).bucket_ordinal(*this);
}

size_t Entity::synchronized_count() const
{
  return BulkData::get(*this).synchronized_count(*this);
}

int Entity::owner_rank() const
{
  return BulkData::get(*this).parallel_owner_rank(*this);
}

#endif

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

BOOST_STATIC_ASSERT(( (int)MetaData::NODE_RANK == (int)Entity::NODE ));
BOOST_STATIC_ASSERT(( (int)MetaData::EDGE_RANK == (int)Entity::EDGE ));
BOOST_STATIC_ASSERT(( (int)MetaData::FACE_RANK == (int)Entity::FACE ));
BOOST_STATIC_ASSERT(( (int)MetaData::ELEMENT_RANK == (int)Entity::ELEMENT ));

} // namespace mesh
} // namespace stk
