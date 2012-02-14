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

#ifdef SIERRA_MIGRATION
namespace sierra {
namespace Fmwk {

const unsigned int INVALID_LOCAL_ID = std::numeric_limits<unsigned int>::max();
const stk::mesh::RelationIterator INVALID_RELATION_ITR;

unsigned get_derived_type(const stk::mesh::Entity&);

}
}
#endif

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
  return print_entity_key(MetaData::get(entity),
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

#ifdef SIERRA_MIGRATION

void Entity::internal_verify_initialization_invariant()
{
  // If this MeshObj has a proper ID (fully initialized), then the id should match
  // the id in the entity-key; otherwise they should not match.
  ThrowRequire( !(m_global_id < 0 && key().id() == static_cast<uint64_t>(m_global_id)) &&
                !(m_global_id > 0 && key().id() != static_cast<uint64_t>(m_global_id)) );
}

// ---------------------------------------------------------------------

void Entity::internal_swap_in_real_entity(const int globalId)
{
  ThrowRequire(globalId > 0);
  m_global_id  = globalId;

  BulkData::get(*this).change_entity_id(globalId, *this);

  internal_verify_initialization_invariant();

  internal_verify_meshobj_invariant();
}

// ---------------------------------------------------------------------

void Entity::reserve_relation(const unsigned num)
{
  if (num == 0 && m_relations.empty()) {
    std::vector<Relation> tmp;
    m_relations.swap(tmp); // clear memory of m_relations.
  }
  else {
    m_relations.reserve(num);
  }
}

// ---------------------------------------------------------------------

namespace {

// In order to preserve relation order, we need to use fmwk-style relation
// comparing when dealing with a fmwk-managed relation; otherwise, use
// stk-style relation ordering.
struct relation_compare
{
  relation_compare(bool use_stk_compare) : m_use_stk_compare(use_stk_compare) {}

  bool operator()(const Relation& lhs, const Relation& rhs) const
  {
    if (m_use_stk_compare) {
      return lhs < rhs;
    }
    else {
      // Fmwk version of relation comparison
      if ( lhs.entity_rank() < rhs.entity_rank()) return true;
      if ( rhs.entity_rank() < lhs.entity_rank()) return false;

      if ( lhs.getRelationType() < rhs.getRelationType() ) return true;
      if ( rhs.getRelationType() < lhs.getRelationType() ) return false;

      return lhs.getOrdinal() < rhs.getOrdinal();
    }
  }

  bool m_use_stk_compare;
};

}

RelationIterator Entity::find_relation(const Relation& relation) const
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

  const Relation::RelationType relation_type = relation.getRelationType();

  RelationIterator rel = std::lower_bound(internal_begin_relation(relation_type),
                                          internal_end_relation(relation_type),
                                          relation,
                                          relation_compare(internal_is_handled_generically(relation_type)));

  while (rel != internal_end_relation(relation_type) &&
         rel->entity_rank()     == relation.entity_rank() &&
         rel->getRelationType() == relation.getRelationType() &&
         rel->getOrdinal()      == relation.getOrdinal() &&
         rel->getMeshObj()      != relation.getMeshObj())
    ++rel;

  return rel;
}

// ---------------------------------------------------------------------

bool Entity::update_relation(
  const RelationIterator        ir ,
  const bool                    back_rel_flag) const
{
  const Relation::RelationType relType = ir->getRelationType();
  ThrowAssert(verify_relation_ordering(internal_begin_relation(relType), internal_end_relation(relType)));

  if (!internal_is_handled_generically(relType)) {
    Entity & meshObj = *ir->getMeshObj() ;

    const bool real_back_rel_flag = back_rel_flag || (relType == Relation::USES || relType == Relation::USED_BY);

    const Relation::RelationType backRelType = back_relation_type(relType);

    ThrowAssert(verify_relation_ordering(meshObj.internal_begin_relation(backRelType), meshObj.internal_end_relation(backRelType)));

    // Create the corresponding back relation to ir
    Relation backRel_obj(const_cast<Entity*>(this), backRelType, ir->getOrdinal(), ir->getOrientation());
    RelationIterator backRel_itr = meshObj.find_relation(backRel_obj);

    const bool exists = backRel_itr != meshObj.internal_end_relation(backRelType) && *backRel_itr == backRel_obj;

    if (exists && !real_back_rel_flag) {
      // Remove back relation and increment the counter

      meshObj.erase_and_clear_if_empty(backRel_itr);

      //ThrowAssert(sierra::Fmwk::get_derived_type(meshObj) != Entity::ELEMENT);

      meshObj.inc_connection();
    }
    else if (!exists && real_back_rel_flag) {
      // Insert back relation

      const unsigned k = backRel_itr - meshObj.internal_begin_relation(backRelType) ;

      // 'relations' may change

      meshObj.reserve_relation(meshObj.m_relations.size() + 1);

      meshObj.m_relations.insert(meshObj.m_relations.begin() + k, backRel_obj);

      //ThrowAssert(sierra::Fmwk::get_derived_type(meshObj) != Entity::ELEMENT);

      meshObj.dec_connection();
    }

    ThrowAssert(verify_relation_ordering(meshObj.internal_begin_relation(relType), meshObj.internal_end_relation(relType)));
  }

  return true;
}

// ---------------------------------------------------------------------

void Entity::erase_and_clear_if_empty(RelationIterator rel_itr)
{
  ThrowRequire(!internal_is_handled_generically(rel_itr->getRelationType()));

  m_relations.erase(m_relations.begin() + (rel_itr - m_relations.begin())); // Need to convert to non-const iterator

  if (m_relations.empty()) {
    reserve_relation(0);
  }
}

// ---------------------------------------------------------------------

void Entity::internal_verify_meshobj_invariant() const
{
  PairIterRelation stk_relations = relations();
  for ( ; !stk_relations.empty(); ++stk_relations ) {
    ThrowRequireMsg(stk_relations->getMeshObj() != NULL, "Problem with: " << *stk_relations);
  }

  for (std::vector<Relation>::const_iterator itr = m_relations.begin(), end = m_relations.end(); itr != end; ++itr) {
    ThrowRequireMsg(itr->getMeshObj() != NULL, "Problem with: " << *itr);
  }
}

// ---------------------------------------------------------------------

void Entity::set_relation_orientation(RelationIterator rel, unsigned orientation)
{
  const Relation::RelationType  backRelType = back_relation_type(rel->getRelationType());

  Entity & meshObj = *rel->getMeshObj();
  Relation backRel_obj(const_cast<Entity*>(this), backRelType, rel->getOrdinal(), rel->getOrientation());
  RelationIterator backRel_itr = meshObj.find_relation(backRel_obj);

  const bool exists = backRel_itr != meshObj.internal_end_relation(backRelType) && *backRel_itr == backRel_obj;
  ThrowRequire(exists);

  // Allow clients to make changes to orientation
  // Orientations do not affect Relation ordering, so this is safe.
  const_cast<Relation*>(&*rel)->setOrientation(orientation);
  const_cast<Relation*>(&*backRel_itr)->setOrientation(orientation);
}

#endif

} // namespace mesh
} // namespace stk

