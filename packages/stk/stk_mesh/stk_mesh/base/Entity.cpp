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

const unsigned int INVALID_LOCAL_ID = std::numeric_limits<unsigned int>::max();
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
  ThrowAssert(!impl::internal_is_handled_generically(rel_itr->getRelationType()));

  RelationVector& aux_rels = aux_relations();
  aux_rels.erase(aux_rels.begin() + (rel_itr - aux_rels.begin())); // Need to convert to non-const iterator

  if (aux_rels.empty()) {
    reserve_relation(0);
  }
}

void Entity::internal_verify_initialization_invariant()
{
#ifndef NDEBUG
  int my_global_id = global_id();
  EntityKey my_key = key();
#endif
  ThrowAssert ( !(my_global_id < 0 && my_key.id() == static_cast<EntityId>(my_global_id)) &&
                !(my_global_id > 0 && my_key.id() != static_cast<EntityId>(my_global_id)) );
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
  ThrowAssert(!impl::internal_is_handled_generically(rel->getRelationType()));

  const RelationType backRelType = back_relation_type(rel->getRelationType());

  Entity meshObj = rel->entity();
  Relation backRel_obj(entity_rank(), *this, backRelType, rel->getOrdinal(), rel->getOrientation());
  RelationIterator backRel_itr = meshObj.find_aux_relation(backRel_obj);

  ThrowRequire(backRel_itr != sierra::Fmwk::INVALID_RELATION_ITR);

  // Allow clients to make changes to orientation
  // Orientations do not affect Relation ordering, so this is safe.
  const_cast<Relation*>(&*rel)->setOrientation(orientation);
  const_cast<Relation*>(&*backRel_itr)->setOrientation(orientation);
}

void Entity::set_relation_orientation(Entity meshObj, ConnectivityOrdinal ord, unsigned orientation)
{
  BulkData& bulk      = BulkData::get(*this);

  Entity const*              fwd_rels  = bulk.begin(*this, meshObj.entity_rank());
  ConnectivityOrdinal const* fwd_ords  = bulk.begin_ordinals(*this, meshObj.entity_rank());
  Permutation *              fwd_perms = const_cast<Permutation*>(bulk.begin_permutations(*this, meshObj.entity_rank()));
  const int                  num_fwd   = bulk.num_connectivity(*this, meshObj.entity_rank());

  Entity const*              back_rels  = bulk.begin(meshObj, entity_rank());
  ConnectivityOrdinal const* back_ords  = bulk.begin_ordinals(meshObj, entity_rank());
  Permutation *              back_perms = const_cast<Permutation*>(bulk.begin_permutations(meshObj, entity_rank()));
  const int                  num_back   = bulk.num_connectivity(meshObj,entity_rank());

  // Find and change fwd connectivity
  for (int i = 0; i < num_fwd; ++i, ++fwd_rels, ++fwd_ords, ++fwd_perms) {
    // Allow clients to make changes to orientation
    // Orientations do not affect Relation ordering, so this is safe.
    if (*fwd_rels == meshObj && *fwd_ords == ord) {
      *fwd_perms = static_cast<Permutation>(orientation);
    }
  }

  // Find and change back connectivity
  for (int i = 0; i < num_back; ++i, ++back_rels, ++back_ords, ++back_perms) {
    // Allow clients to make changes to orientation
    // Orientations do not affect Relation ordering, so this is safe.
    if (*back_rels == *this && *back_ords == ord) {
      *back_perms = static_cast<Permutation>(orientation);
    }
  }
}

// ---------------------------------------------------------------------

void Entity::reserve_relation(const unsigned num)
{
  if (num == 0 && aux_relations().empty()) {
    RelationVector tmp;
    aux_relations().swap(tmp); // clear memory of m_relations.
  }
  else {
    aux_relations().reserve(num);
  }
}

RelationIterator Entity::find_aux_relation(const Relation& relation) const
{
  // Extremely hacky: It would be better to set up the < operator for relations so that lower_bound
  // can return the desired iterator, but any sane definition would probably force a change in
  // relation ordering and that's more than I'm willing to take on now.
  //
  // The current semantics for relation-searching is as follows:
  // Ordered based on derived_type, relation_type, and ordinal in descending precedence
  //   If multiple relations have the same derived_type, relation_type, and ordinal, a linear
  //   scan takes place looking for a matching meshobj. If no such meshobj was found, then
  //   we are left with an iterator pointing to the first relation with a different derived_type,
  //   relation_type, or ordinal. To sum up, the result of the search can either be equivalent to
  //   lower_bound OR upper_bound depending upon the state of the relations... YUCK!

  ThrowAssert(!impl::internal_is_handled_generically(relation.getRelationType()));

  for (RelationIterator rel = aux_relation_begin(); rel != aux_relation_end(); ++rel) {
    if (same_specification(*rel, relation) && rel->entity() != relation.entity()) {
      return rel;
    }
  }

  return sierra::Fmwk::INVALID_RELATION_ITR;
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
