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

const unsigned int INVALID_LOCAL_ID = stk::mesh::InvalidLocalId;
const stk::mesh::RelationIterator INVALID_RELATION_ITR = dummy_vector.end(); // Some STL implementation use POD for iterators

#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
unsigned get_derived_type(const stk::mesh::Entity entity)
{
  return stk::mesh::BulkData::get(entity).entity_rank(entity);
}
#endif

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

#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS

void Entity::compress_relation_capacity()
{
  return BulkData::get(*this).compress_relation_capacity(*this);
}

int Entity::global_id() const
{
  return BulkData::get(*this).global_id(*this);
}

unsigned Entity::local_id() const
{
  return BulkData::get(*this).local_id(*this);
}

void Entity::set_local_id(unsigned int l_id)
{
  BulkData::get(*this).set_local_id(*this, l_id);
}

int Entity::owner_processor_rank() const
{
  return BulkData::get(*this).parallel_owner_rank(*this);
}

void Entity::set_owner_processor_rank(int owner)
{
  BulkData::get(*this).set_parallel_owner_rank(*this, owner);
}

void Entity::set_owner_rank(int owner)
{
  BulkData::get(*this).set_parallel_owner_rank(*this, owner);
}

void Entity::erase_and_clear_if_empty(RelationIterator rel_itr)
{
  BulkData::get(*this).erase_and_clear_if_empty(*this, rel_itr);
}

void Entity::internal_verify_initialization_invariant()
{
  BulkData::get(*this).internal_verify_initialization_invariant(*this);
}

unsigned Entity::size_connection() const
{
  return BulkData::get(*this).get_connect_count(*this);
}

unsigned Entity::inc_connection()
{
  BulkData& bulk_data = BulkData::get(*this);
  int count = bulk_data.get_connect_count(*this) + 1;
  bulk_data.set_connect_count(*this, count);
  return count;
}

unsigned Entity::dec_connection()
{
  BulkData& bulk_data = BulkData::get(*this);
  int count = bulk_data.get_connect_count(*this) - 1;
  bulk_data.set_connect_count(*this, count);
  return count;
}

RelationIterator Entity::aux_relation_begin() const
{
  return BulkData::get(*this).aux_relations(*this).begin();
}

RelationIterator Entity::aux_relation_end() const
{
  return BulkData::get(*this).aux_relations(*this).end();
}

RelationVector& Entity::aux_relations()
{
  return BulkData::get(*this).aux_relations(*this);
}

const RelationVector& Entity::aux_relations() const
{
  return BulkData::get(*this).aux_relations(*this);
}

void Entity::set_shared_attr(const void* attr)
{
  BulkData::get(*this).set_shared_attr(*this, attr);
}

const void* Entity::get_shared_attr() const
{
  return BulkData::get(*this).get_shared_attr(*this);
}

void Entity::set_relation_orientation(RelationIterator rel, unsigned orientation)
{
  BulkData::get(*this).set_relation_orientation(*this, rel, orientation);
}

void Entity::set_relation_orientation(Entity meshObj, ConnectivityOrdinal ord, unsigned orientation)
{
  BulkData::get(*this).set_relation_orientation(*this, meshObj, ord, orientation);
}

void Entity::reserve_relation(const unsigned num)
{
  BulkData::get(*this).reserve_relation(*this, num);
}

#endif // allow deprecated

// ---------------------------------------------------------------------


#endif // SIERRA_MIGRATION

BOOST_STATIC_ASSERT(( (int)MetaData::NODE_RANK == (int)Entity::NODE ));
BOOST_STATIC_ASSERT(( (int)MetaData::EDGE_RANK == (int)Entity::EDGE ));
BOOST_STATIC_ASSERT(( (int)MetaData::FACE_RANK == (int)Entity::FACE ));
BOOST_STATIC_ASSERT(( (int)MetaData::ELEMENT_RANK == (int)Entity::ELEMENT ));

} // namespace mesh
} // namespace stk
