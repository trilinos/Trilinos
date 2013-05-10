/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

// dan is newb

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

unsigned get_derived_type(const stk::mesh::Entity entity)
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return stk::mesh::BulkData::get(entity).entity_rank(entity);
#else
  ThrowErrorMsg("sierra::Fmwk::get_derived_type(const stk::mesh::Entity) has been deprecated.");
  return 0;
#endif
}

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

bool Entity::is_valid() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return (m_value != 0) && BulkData::get(*this).is_valid(*this);
#else
  ThrowErrorMsg("Entity::is_valid() has been deprecated");
  return false;
#endif
}

EntityState Entity::state() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).state(*this);
#else
  ThrowErrorMsg("Entity::state() has been deprecated");
  static EntityState *bad = 0;
  return *bad;
#endif
}

EntityRank Entity::entity_rank() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).entity_rank(*this);
#else
  ThrowErrorMsg("Entity::entity_rank() has been deprecated");
  static EntityRank *bad = 0;
  return *bad;
#endif
}

EntityId Entity::identifier() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).identifier(*this);
#else
  ThrowErrorMsg("Entity::identifier() has been deprecated");
  static EntityId *bad;
  return *bad;
#endif
}

const EntityKey Entity::key() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).entity_key(*this);
#else
  ThrowErrorMsg("Entity::key() has been deprecated");
  static EntityKey *bad = 0;
  return *bad;
#endif
}

Bucket & Entity::bucket() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).bucket(*this);
#else
  ThrowErrorMsg("Entity::bucket() has been deprecated");
  static Bucket *bad = 0;
  return *bad;
#endif
}

Bucket * Entity::bucket_ptr() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).bucket_ptr(*this);
#else
  ThrowErrorMsg("Entity::bucket_ptr() has been deprecated");
  static Bucket *bad = 0;
  return bad;
#endif
}

unsigned Entity::bucket_ordinal() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).bucket_ordinal(*this);
#else
  ThrowErrorMsg("Entity::bucket_ordinal() has been deprecated");
  return 0;
#endif
}

size_t Entity::synchronized_count() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).synchronized_count(*this);
#else
  ThrowErrorMsg("Entity::synchronize_count() has been deprecated");
  return 0;
#endif
}

PairIterRelation Entity::relations() const
{
  ThrowErrorMsg("Entity::relations() has been deprecated");
  PairIterRelation dummy;
  return dummy;
}

PairIterRelation Entity::relations( EntityRank type ) const
{
  ThrowErrorMsg("Entity::relations(EntityRank) has been deprecated");
  PairIterRelation dummy;
  return dummy;
}

PairIterRelation Entity::node_relations() const
{
  ThrowErrorMsg("Entity::node_relations() has been deprecated");
  PairIterRelation dummy;
  return dummy;
}

RelationIterator Entity::node_relation(unsigned ordinal) const
{
  ThrowErrorMsg("Entity::node_relation(unsigned) has been deprecated");
  RelationIterator dummy;
  return dummy;
}

int Entity::owner_rank() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).parallel_owner_rank(*this);
#else
  ThrowErrorMsg("Entity::owner_rank() has been deprecated");
  return 0;
#endif
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

void Entity::compress_relation_capacity()
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).compress_relation_capacity(*this);
#else
  ThrowErrorMsg("Entity::compress_relation_capacity() has been deprecated");
#endif
}

int Entity::global_id() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).global_id(*this);
#else
  ThrowErrorMsg("Entity::global_id() has been deprecated");
  return 0;
#endif
}

unsigned Entity::local_id() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).local_id(*this);
#else
  ThrowErrorMsg("Entity::local_id() has been deprecated");
  return 0;
#endif
}

void Entity::set_local_id(unsigned int l_id)
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  BulkData::get(*this).set_local_id(*this, l_id);
#else
  ThrowErrorMsg("Entity::set_local_id(int) has been deprecated");
#endif
}

int Entity::owner_processor_rank() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).parallel_owner_rank(*this);
#else
  ThrowErrorMsg("Entity::owner_processor_rank() has been deprecated");
  return 0;
#endif
}

void Entity::set_owner_processor_rank(int owner)
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  BulkData::get(*this).set_parallel_owner_rank(*this, owner);
#else
  ThrowErrorMsg("Entity::set_owner_processor_rank(.) has been deprecated");
#endif
}

void Entity::set_owner_rank(int owner)
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  BulkData::get(*this).set_parallel_owner_rank(*this, owner);
#else
  ThrowErrorMsg("Entity::set_owner_rank(.) has been deprecated");
#endif
}

void Entity::erase_and_clear_if_empty(RelationIterator rel_itr)
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  ThrowAssert(!impl::internal_is_handled_generically(rel_itr->getRelationType()));

  RelationVector& aux_rels = aux_relations();
  aux_rels.erase(aux_rels.begin() + (rel_itr - aux_rels.begin())); // Need to convert to non-const iterator

  if (aux_rels.empty()) {
    reserve_relation(0);
  }
#else
  ThrowErrorMsg("Entity::erase_and_clear_if_empty(..) has been deprecated");
#endif
}

void Entity::internal_verify_initialization_invariant()
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
#ifndef NDEBUG
  int my_global_id = global_id();
  EntityKey my_key = key();
#endif
  ThrowAssert ( !(my_global_id < 0 && my_key.id() == static_cast<EntityId>(my_global_id)) &&
                !(my_global_id > 0 && my_key.id() != static_cast<EntityId>(my_global_id)) );
#else
  ThrowErrorMsg("Entity::internal_verify_initialization_invariant() has been deprecated");
#endif
}

unsigned Entity::size_connection() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).get_connect_count(*this);
#else
  ThrowErrorMsg("Entity::size_connection() has been deprecated");
  return 0;
#endif
}

unsigned Entity::inc_connection()
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  BulkData& bulk_data = BulkData::get(*this);
  int count = bulk_data.get_connect_count(*this) + 1;
  bulk_data.set_connect_count(*this, count);
  return count;
#else
  ThrowErrorMsg("Entity::inc_connection() has been deprecated");
  return 0;
#endif
}

unsigned Entity::dec_connection()
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  BulkData& bulk_data = BulkData::get(*this);
  int count = bulk_data.get_connect_count(*this) - 1;
  bulk_data.set_connect_count(*this, count);
  return count;
#else
  ThrowErrorMsg("Entity::dec_connection() has been deprecated");
  return 0;
#endif
}

RelationIterator Entity::aux_relation_begin() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).aux_relations(*this).begin();
#else
  ThrowErrorMsg("Entity::aux_relation_begin() has been deprecated");
  RelationIterator dummy;
  return dummy;
#endif
}

RelationIterator Entity::aux_relation_end() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).aux_relations(*this).end();
#else
  ThrowErrorMsg("Entity::aux_relation_end() has been deprecated");
  RelationIterator dummy;
  return dummy;
#endif
}

RelationVector& Entity::aux_relations()
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).aux_relations(*this);
#else
  ThrowErrorMsg("Entity::aux_relations has been deprecated");
  static RelationVector *bad;
  return *bad;
#endif
}

const RelationVector& Entity::aux_relations() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).aux_relations(*this);
#else
  ThrowErrorMsg("Entity::aux_relations has been deprecated");
  static RelationVector *bad;
  return *bad;
#endif
}

RelationIterator Entity::internal_begin_relation(const RelationType relation_type) const
{
  ThrowErrorMsg("Entity::internal_begin_relation(.) has been deprecated");
  RelationIterator dummy;
  return dummy;
}

RelationIterator Entity::internal_end_relation(const RelationType relation_type) const
{
  ThrowErrorMsg("Entity::internal_end_relation(.) has been deprecated");
  RelationIterator dummy;
  return dummy;
}

void Entity::set_shared_attr(const void* attr)
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  BulkData::get(*this).set_shared_attr(*this, attr);
#else
  ThrowErrorMsg("Entity::set_shared_attr(.) has been deprecated");
#endif
}

const void* Entity::get_shared_attr() const
{
#ifdef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  return BulkData::get(*this).get_shared_attr(*this);
#else
  ThrowErrorMsg("Entity::get_shared_attr() has been deprecated;");
  return 0;
#endif
}

void Entity::set_relation_orientation(RelationIterator rel, unsigned orientation)
{
#ifndef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  ThrowErrorMsg("Entity::set_relation_orientation(RelationIterator, unsigned) needs to be re-factored away for POD stk::mesh::Entity");
#else
  if (!impl::internal_is_handled_generically(rel->getRelationType())) {
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
  else {
    Entity meshObj      = rel->entity();
    const unsigned ord  = rel->getOrdinal();
    BulkData& bulk      = BulkData::get(*this);

    Entity const*              fwd_rels  = bulk.begin_entities(*this, meshObj.entity_rank());
    ConnectivityOrdinal const* fwd_ords  = bulk.begin_ordinals(*this, meshObj.entity_rank());
    Permutation *              fwd_perms = const_cast<Permutation*>(bulk.begin_permutations(*this, meshObj.entity_rank()));
    const int                  num_fwd   = bulk.num_connectivity(*this, meshObj.entity_rank());

    Entity const*              back_rels  = bulk.begin_entities(meshObj, entity_rank());
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
#endif
}

// ---------------------------------------------------------------------

void Entity::reserve_relation(const unsigned num)
{
#ifndef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  ThrowErrorMsg("Entity::reserve_relation(unsigned) needs to be re-factored away for POD stk::mesh::Entity");
#else
  if (num == 0 && aux_relations().empty()) {
    RelationVector tmp;
    aux_relations().swap(tmp); // clear memory of m_relations.
  }
  else {
    aux_relations().reserve(num);
  }
#endif
}

RelationIterator Entity::find_aux_relation(const Relation& relation) const
{
#ifndef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  static RelationIterator dummy_rel;
  ThrowErrorMsg("Entity::find_relation has been deprecated for POD stk::mesh::Entity.");
  return rel;
#else
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
#endif
}

// ---------------------------------------------------------------------

void Entity::internal_verify_meshobj_invariant() const
{
#ifndef STK_MESH_ALLOW_DEPRECATED_ENTITY_FNS
  ThrowErrorMsg("Entity::internal_verify_meshobj_invariant() needs to be re-factored away for POD stk::mesh::Entity");
#else
  PairIterRelation stk_relations = relations();
  for ( ; !stk_relations.empty(); ++stk_relations ) {
    ThrowRequireMsg(stk_relations->entity().is_valid(), "Problem with: " << *stk_relations);
  }

  const RelationVector& aux_relations = this->aux_relations();
  for (RelationVector::const_iterator itr = aux_relations.begin(), end = aux_relations.end(); itr != end; ++itr) {
    ThrowRequireMsg(itr->entity().is_valid(), "Problem with: " << *itr);
  }
#endif
}

// ---------------------------------------------------------------------


#endif // SIERRA_MIGRATION

BOOST_STATIC_ASSERT(( (int)MetaData::NODE_RANK == (int)Entity::NODE ));
BOOST_STATIC_ASSERT(( (int)MetaData::EDGE_RANK == (int)Entity::EDGE ));
BOOST_STATIC_ASSERT(( (int)MetaData::FACE_RANK == (int)Entity::FACE ));
BOOST_STATIC_ASSERT(( (int)MetaData::ELEMENT_RANK == (int)Entity::ELEMENT ));

} // namespace mesh
} // namespace stk
