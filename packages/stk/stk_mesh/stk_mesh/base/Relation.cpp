// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_mesh/base/Relation.hpp>
#include <stddef.h>                     // for NULL
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/Bucket.hpp>     // for Bucket
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, get_connectivity
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <utility>                      // for pair
#include "stk_mesh/base/ConnectivityMap.hpp"  // for ConnectivityMap
#include "stk_mesh/base/Part.hpp"       // for Part, insert_ordinal, etc
#include "stk_mesh/base/Types.hpp"      // for EntityRank, OrdinalVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowAssertMsg, etc


namespace stk {
namespace mesh {


Relation::RawRelationType & Relation::RawRelationType::operator =(const Relation::RawRelationType & rhs)
{
    value = rhs.value;
    return *this;
}

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
  EntityVector temp_entities;
  Entity const* irels_j = NULL;
  int num_conn = 0;
  for (Entity const *rels_left = rels_begin ; rels_left != rels_end ; ++rels_left )
  {
    // Do all input entities have a relation to this entity ?

    Entity const e = *rels_left;
    EntityRank erank = mesh.entity_rank(e);

    std::vector<Entity>::const_iterator i = i_beg ;
    for ( ; i != i_end ; ++i )
    {
      if (mesh.connectivity_map().valid(mesh.entity_rank(*i), erank)) {
        num_conn = mesh.num_connectivity(*i, erank);
        irels_j  = mesh.begin(*i, erank);
      }
      else {
        num_conn = get_connectivity(mesh, *i, erank, temp_entities);
        irels_j  = &*temp_entities.begin();
      }

      Entity const *irels_end = irels_j + num_conn;

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

} // namespace

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
    const EntityRank end_rank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());

    std::vector<Entity>::const_iterator next_i = i + 1;
    EntityVector temp_entities;
    Entity const* rels_begin = NULL;
    int num_conn = 0;
    for (EntityRank rank = stk::topology::BEGIN_RANK; rank < end_rank; ++rank)
    {
      if (mesh.connectivity_map().valid(ibucket.entity_rank(), rank)) {
        num_conn   = mesh.num_connectivity(ibucket[ibordinal], rank);
        rels_begin = mesh.begin(ibucket[ibordinal], rank);
      }
      else {
        num_conn   = get_connectivity(mesh, ibucket[ibordinal], rank, temp_entities);
        rels_begin = &*temp_entities.begin();
      }

      get_entities_through_relations(mesh, rels_begin, rels_begin + num_conn, next_i, j,
                                     entities_related);
    }
  }
}

void get_entities_through_relations(
  const BulkData& mesh,
  const std::vector<Entity> & entities ,
        EntityRank              entities_related_rank ,
        std::vector<Entity> & entities_related )
{
  entities_related.clear();

  if ( ! entities.empty() ) {

          std::vector<Entity>::const_iterator i = entities.begin();
    const std::vector<Entity>::const_iterator j = entities.end();

    EntityVector temp_entities;
    Entity const* rel_entities = NULL;
    int num_rels = 0;
    if (mesh.connectivity_map().valid(mesh.entity_rank(*i), entities_related_rank)) {
      num_rels     = mesh.num_connectivity(*i, entities_related_rank);
      rel_entities = mesh.begin(*i, entities_related_rank);
    }
    else {
      num_rels    = get_connectivity(mesh, *i, entities_related_rank, temp_entities);
      rel_entities = &*temp_entities.begin();
    }

    ++i;
    get_entities_through_relations(mesh, rel_entities, rel_entities + num_rels, i, j, entities_related);
  }
}

//----------------------------------------------------------------------

void induced_part_membership( const Part & part ,
                              EntityRank entity_rank_from ,
                              EntityRank entity_rank_to ,
                              OrdinalVector & induced_parts)
{
  if ( entity_rank_to < entity_rank_from && part.should_induce(entity_rank_from) ) {
    // Direct relationship:
    insert_ordinal( induced_parts , part.mesh_meta_data_ordinal() );
  }
}

//----------------------------------------------------------------------
//  What are this entity's part memberships that can be deduced from
//  this entity's relationship.  Can only trust 'entity_from' to be
//  accurate if it is owned by the local process.

void induced_part_membership(const BulkData& mesh,
                             const PartVector& all_parts,
                             const Entity entity_from ,
                             const OrdinalVector       & omit ,
                                   EntityRank            entity_rank_to ,
                                   OrdinalVector       & induced_parts)
{
  const Bucket   & bucket_from    = mesh.bucket(entity_from);
  const int      local_proc_rank  = mesh.parallel_rank();
  const EntityRank entity_rank_from = bucket_from.entity_rank();
  const bool dont_check_owner     = mesh.parallel_size() == 1; // critical for fmwk
  ThrowAssert(entity_rank_from > entity_rank_to);

  // Only induce parts for normal (not back) relations. Can only trust
  // 'entity_from' to be accurate if it is owned by the local process.
  if ( dont_check_owner || local_proc_rank == mesh.parallel_owner_rank(entity_from) ) {

    const std::pair<const unsigned *, const unsigned *>
      bucket_superset_ordinals = bucket_from.superset_part_ordinals();

    OrdinalVector::const_iterator omit_begin = omit.begin(),
                                  omit_end   = omit.end();

    // Contributions of the 'from' entity:
    for ( const unsigned * i = bucket_superset_ordinals.first ;
                           i != bucket_superset_ordinals.second ; ++i ) {
      ThrowAssertMsg( *i < all_parts.size(), "Index " << *i << " out of bounds" );
      Part & part = * all_parts[*i] ;

      if ( part.should_induce(entity_rank_from) && ! contains_ordinal( omit_begin, omit_end , *i )) {
        induced_part_membership( part,
                                 entity_rank_from ,
                                 entity_rank_to ,
                                 induced_parts);
      }
    }
  }
}

//----------------------------------------------------------------------

void induced_part_membership(const BulkData& mesh,
                             const Entity entity ,
                             const OrdinalVector & omit ,
                                   OrdinalVector & induced_parts)
{
  ThrowAssertMsg(mesh.is_valid(entity), "BulkData at " << &mesh << " does not know Entity" << entity.local_offset());

  const EntityRank e_rank = mesh.entity_rank(entity);
  const EntityRank end_rank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());
  const PartVector& all_parts = mesh.mesh_meta_data().get_parts();

  EntityVector temp_entities;
  Entity const* rels = NULL;
  int num_rels = 0;
  for (EntityRank irank = static_cast<EntityRank>(e_rank + 1); irank < end_rank; ++irank)
  {
    if (mesh.connectivity_map().valid(e_rank, irank)) {
      num_rels = mesh.num_connectivity(entity, irank);
      rels     = mesh.begin(entity, irank);
    }
    else {
      num_rels = get_connectivity(mesh, entity, irank, temp_entities);
      rels     = &*temp_entities.begin();
    }

    for (int j = 0; j < num_rels; ++j)
    {
      induced_part_membership(mesh, all_parts, rels[j], omit, e_rank, induced_parts);
    }
  }
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk
