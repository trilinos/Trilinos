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

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

std::ostream &
operator << ( std::ostream & s , const Relation & rel )
{
  Entity const e = rel.entity();

  s << "[" << rel.relation_ordinal() << "]->(" << rel.entity_rank()
    << ", " << e.local_offset() << ")";

  return s ;
}

namespace {

void get_entities_through_relations(
  const BulkData &mesh,
  Entity const *rels_begin,
  Entity const *rels_end,
  const std::vector<Entity>::const_iterator i_beg ,
  const std::vector<Entity>::const_iterator i_end ,
  std::vector<Entity> & entities_related )
{
  for (Entity const *rels_left = rels_begin ; rels_left != rels_end ; ++rels_left )
  {
    // Do all input entities have a relation to this entity ?

    Entity const e = *rels_left;
    EntityRank erank = mesh.entity_rank(e);

    std::vector<Entity>::const_iterator i = i_beg ;
    for ( ; i != i_end ; ++i )
    {
      Entity const *irels_j = mesh.begin_entities(*i, erank);
      Entity const *irels_end = mesh.end_entities(*i, erank);

      while ( irels_j != irels_end && e != *irels_j) {
        ++irels_j ;
      }
      if ( irels_j == irels_end ) {
        // Entity *i does not have a relation to Entity e.
        break ;
      }
    }

    if ( i == i_end ) {
      entities_related.push_back( e );
    }
  }
}

inline
void insert_part_and_supersets(OrdinalVector& induced_parts,
                               const Part& part,
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
  const BulkData &mesh,
  const std::vector<Entity> & entities ,
        std::vector<Entity> & entities_related )
{
  entities_related.clear();

  if ( ! entities.empty() ) {
          std::vector<Entity>::const_iterator i = entities.begin();
    const std::vector<Entity>::const_iterator j = entities.end();

    const Bucket &ibucket = mesh.bucket(*i);
    const Ordinal &ibordinal = mesh.bucket_ordinal(*i);
    const EntityRank end_rank = mesh.mesh_meta_data().entity_rank_count();

    std::vector<Entity>::const_iterator next_i = i + 1;
    for (EntityRank rank = stk::topology::BEGIN_RANK; rank < end_rank; ++rank)
    {
      Entity const *rels_begin = ibucket.begin_entities(ibordinal, rank);
      Entity const *rels_end = ibucket.end_entities(ibordinal, rank);
      get_entities_through_relations(mesh, rels_begin, rels_end, next_i, j,
                                     entities_related);
    }
  }
}

void get_entities_through_relations(
  const BulkData& mesh,
  const std::vector<Entity> & entities ,
        unsigned               entities_related_rank ,
        std::vector<Entity> & entities_related )
{
  entities_related.clear();

  if ( ! entities.empty() ) {
          std::vector<Entity>::const_iterator i = entities.begin();
    const std::vector<Entity>::const_iterator j = entities.end();

    Entity const *rels_begin = mesh.begin_entities(*i, entities_related_rank );
    Entity const *rels_end = mesh.end_entities(*i, entities_related_rank );
    ++i;
    get_entities_through_relations(mesh, rels_begin, rels_end, i, j, entities_related);
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

  const bool induced_by_stencil = false;

  return induced_by_type || induced_by_stencil ;
}

//----------------------------------------------------------------------

void induced_part_membership( const Part & part ,
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

  }
}

//----------------------------------------------------------------------
//  What are this entity's part memberships that can be deduced from
//  this entity's relationship.  Can only trust 'entity_from' to be
//  accurate if it is owned by the local process.

void induced_part_membership(const BulkData& mesh, const Entity entity_from ,
                              const OrdinalVector       & omit ,
                                    unsigned           entity_rank_to ,
                                    RelationIdentifier relation_identifier ,
                                    OrdinalVector       & induced_parts,
                                    bool include_supersets)
{
  const Bucket   & bucket_from    = mesh.bucket(entity_from);
  const int      local_proc_rank  = mesh.parallel_rank();
  const unsigned entity_rank_from = bucket_from.entity_rank();
  const bool dont_check_owner     = mesh.parallel_size() == 1; // critical for fmwk

  // Only induce parts for normal (not back) relations. Can only trust
  // 'entity_from' to be accurate if it is owned by the local process.
  if ( entity_rank_to < entity_rank_from &&
       ( dont_check_owner || local_proc_rank == mesh.parallel_owner_rank(entity_from) ) ) {
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

void induced_part_membership(const BulkData& mesh, const Entity entity ,
                              const OrdinalVector & omit ,
                                    OrdinalVector & induced_parts,
                                    bool include_supersets)
{
  const MeshIndex e_idx = mesh.mesh_index(entity);
  ThrowAssertMsg(e_idx.bucket, "BulkData at " << &mesh << " does not know Entity" << mesh.identifier(entity));
  const Bucket &e_bucket = *e_idx.bucket;
  const Ordinal e_ordinal = e_idx.bucket_ordinal;

  const EntityRank e_rank = mesh.entity_rank(entity);
  const EntityRank end_rank = mesh.mesh_meta_data().entity_rank_count();

  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    int num_rels = e_bucket.num_connectivity(e_ordinal, irank);
    Entity const *rels = e_bucket.begin_entities(e_ordinal, irank);
    ConnectivityOrdinal const *ords = e_bucket.begin_ordinals(e_ordinal, irank);
    for (int j = 0; j < num_rels; ++j)
    {
      induced_part_membership(mesh, rels[j] , omit , e_rank, ords[j], induced_parts,
                              include_supersets);
    }
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk
