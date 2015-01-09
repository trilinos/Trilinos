/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Relation.hpp>
#include <stk_mesh/base/FieldData.hpp>

namespace stk_classic {
namespace mesh {

//----------------------------------------------------------------------

std::ostream &
operator << ( std::ostream & s , const Relation & rel )
{
  Entity * const e = rel.entity();

  if ( e ) {
    const MetaData & meta_data = MetaData::get(*e);
    s << "[" << rel.identifier() << "]->" ;
    print_entity_key( s , meta_data , e->key() );
  }
  else {
    s << "[" << rel.identifier() << "]->" << rel.entity_rank();
  }

  return s ;
}

//----------------------------------------------------------------------

Relation::Relation( Entity & entity , RelationIdentifier identifier )
  : m_raw_relation( Relation::raw_relation_id( entity.entity_rank() , identifier ) ),
    m_target_entity( & entity )
{
#ifdef SIERRA_MIGRATION
  setRelationType(INVALID);
#endif
}

bool Relation::operator < ( const Relation & rhs ) const
{
  bool result = false;

#ifdef SIERRA_MIGRATION
  if (entity_rank() != rhs.entity_rank()) {
    result = entity_rank() < rhs.entity_rank();
  }
  else if (getRelationType() != rhs.getRelationType()) {
    result = getRelationType() < rhs.getRelationType();
  }
  else if (identifier() != rhs.identifier()) {
    result = identifier() < rhs.identifier();
  }
#else
  if ( m_raw_relation.value != rhs.m_raw_relation.value ) {
    result = m_raw_relation.value < rhs.m_raw_relation.value ;
  }
#endif
  else {
    const EntityKey lhs_key = m_target_entity     ? m_target_entity->key()     : EntityKey();
    const EntityKey rhs_key = rhs.m_target_entity ? rhs.m_target_entity->key() : EntityKey();
    result = lhs_key < rhs_key ;
  }
  return result ;
}

//----------------------------------------------------------------------

#ifdef SIERRA_MIGRATION

Relation::Relation(Entity *obj, const unsigned relation_type, const unsigned ordinal, const unsigned orient)
  :
  m_raw_relation( Relation::raw_relation_id( obj->entity_rank(), ordinal )),
  m_attribute( (relation_type << fmwk_orientation_digits) | orient ),
  m_target_entity(obj)
{
  ThrowAssertMsg( orient <= fmwk_orientation_mask,
                  "orientation " << orient << " exceeds maximum allowed value");
}

void Relation::setMeshObj(Entity *object)
{
  if (object != NULL) {
    m_raw_relation = Relation::raw_relation_id( object->entity_rank(), identifier() );
  }
  m_target_entity = object;
}

#endif

namespace {

void get_entities_through_relations(
  PairIterRelation rel ,
  const std::vector<Entity*>::const_iterator i_beg ,
  const std::vector<Entity*>::const_iterator i_end ,
  std::vector<Entity*> & entities_related )
{
  for ( ; rel.first != rel.second ; ++rel.first ) {

    // Do all input entities have a relation to this entity ?

    Entity * const e = rel.first->entity();

    std::vector<Entity*>::const_iterator i = i_beg ;

    for ( ; i != i_end ; ++i ) {
      PairIterRelation r = (*i)->relations();
      while ( r.first != r.second && e != r.first->entity() ) {
        ++r.first ;
      }
      if ( r.first == r.second ) { break ; }
    }

    if ( i == i_end ) {
      entities_related.push_back( e );
    }
  }
}

inline
void insert_part_and_supersets(OrdinalVector& induced_parts,
                               Part& part,
                               bool include_supersets)
{
  insert_ordinal( induced_parts , part.mesh_meta_data_ordinal() );

  // In order to preserve superset/subset consistency we should add supersets of
  // induced parts to the induced part lists. Unfortunately, this opens up an ambiguity
  // where, when a relation is removed, we cannot know if an unranked superset
  // part should be removed.
  if (include_supersets) {
    const PartVector & supersets = part.supersets();
    for (PartVector::const_iterator itr = supersets.begin(), end = supersets.end(); itr != end; ++itr) {
      insert_ordinal( induced_parts, (*itr)->mesh_meta_data_ordinal() );
    }
  }
}

}

void get_entities_through_relations(
  const std::vector<Entity*> & entities ,
        std::vector<Entity*> & entities_related )
{
  entities_related.clear();

  if ( ! entities.empty() ) {
          std::vector<Entity*>::const_iterator i = entities.begin();
    const std::vector<Entity*>::const_iterator j = entities.end();

    PairIterRelation rel = (*i)->relations(); ++i ;

    get_entities_through_relations( rel , i , j , entities_related );
  }
}

void get_entities_through_relations(
  const std::vector<Entity*> & entities ,
        unsigned               entities_related_rank ,
        std::vector<Entity*> & entities_related )
{
  entities_related.clear();

  if ( ! entities.empty() ) {
          std::vector<Entity*>::const_iterator i = entities.begin();
    const std::vector<Entity*>::const_iterator j = entities.end();

    PairIterRelation rel = (*i)->relations( entities_related_rank ); ++i ;

    get_entities_through_relations( rel , i , j , entities_related );
  }
}

//----------------------------------------------------------------------

/** \brief  Query if a member entity of the given entity type
 *          has an induced membership.
 */
bool membership_is_induced( const Part & part , unsigned entity_rank )
{
  const MetaData & meta = MetaData::get(part);

  const bool induced_by_type =
     entity_rank < part.primary_entity_rank() &&
                   part.primary_entity_rank() < meta.entity_rank_count() ;

  const bool induced_by_stencil =
    ! part.relations().empty() &&
      part.relations().begin()->m_target == & part ;

  return induced_by_type || induced_by_stencil ;
}

//----------------------------------------------------------------------

void induced_part_membership( Part & part ,
                              unsigned entity_rank_from ,
                              unsigned entity_rank_to ,
                              RelationIdentifier relation_identifier ,
                              OrdinalVector & induced_parts,
                              bool include_supersets)
{
  if ( entity_rank_to < entity_rank_from &&
       part.primary_entity_rank() == entity_rank_from ) {

    // Direct relationship:

    insert_part_and_supersets( induced_parts , part, include_supersets );

    // Stencil relationship where 'part' is the root:
    // The 'target' should not have subsets or supersets.

    const std::vector<PartRelation> & part_rel = part.relations();

    for ( std::vector<PartRelation>::const_iterator
          j = part_rel.begin() ; j != part_rel.end() ; ++j ) {

      if ( & part == j->m_root &&
           0 <= (* j->m_function)( entity_rank_from , entity_rank_to ,
                                   relation_identifier ) ) {
        insert_part_and_supersets( induced_parts , * j->m_target, include_supersets );
      }
    }
  }
}

//----------------------------------------------------------------------
//  What are this entity's part memberships that can be deduced from
//  this entity's relationship.  Can only trust 'entity_from' to be
//  accurate if it is owned by the local process.

void induced_part_membership( const Entity           & entity_from ,
                              const OrdinalVector       & omit ,
                                    unsigned           entity_rank_to ,
                                    RelationIdentifier relation_identifier ,
                                    OrdinalVector       & induced_parts,
                                    bool include_supersets)
{
  const Bucket   & bucket_from    = entity_from.bucket();
  const BulkData & mesh           = BulkData::get(bucket_from);
  const unsigned local_proc_rank  = mesh.parallel_rank();
  const unsigned entity_rank_from = entity_from.entity_rank();

  // Only induce parts for normal (not back) relations. Can only trust
  // 'entity_from' to be accurate if it is owned by the local process.
  if ( entity_rank_to < entity_rank_from &&
       local_proc_rank == entity_from.owner_rank() ) {
    const PartVector & all_parts   = mesh.mesh_meta_data().get_parts();

    const std::pair<const unsigned *, const unsigned *>
      bucket_superset_ordinals = bucket_from.superset_part_ordinals();

    OrdinalVector::const_iterator omit_begin = omit.begin(),
                                  omit_end   = omit.end();

    // Contributions of the 'from' entity:
    for ( const unsigned * i = bucket_superset_ordinals.first ;
                           i != bucket_superset_ordinals.second ; ++i ) {
      ThrowAssertMsg( *i < all_parts.size(), "Index " << *i << " out of bounds" );
      Part & part = * all_parts[*i] ;

      if ( part.primary_entity_rank() == entity_rank_from && ! contains_ordinal( omit_begin, omit_end , *i )) {
        induced_part_membership( part,
                                 entity_rank_from ,
                                 entity_rank_to ,
                                 relation_identifier ,
                                 induced_parts,
                                 include_supersets);
      }
    }
  }
}

//----------------------------------------------------------------------

void induced_part_membership( const Entity     & entity ,
                              const OrdinalVector & omit ,
                                    OrdinalVector & induced_parts,
                                    bool include_supersets)
{
  for ( PairIterRelation
        rel = entity.relations() ; ! rel.empty() ; ++rel ) {

    induced_part_membership( * rel->entity() , omit ,
                             entity.entity_rank() ,
                             rel->identifier() ,
                             induced_parts,
                             include_supersets);
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk_classic
