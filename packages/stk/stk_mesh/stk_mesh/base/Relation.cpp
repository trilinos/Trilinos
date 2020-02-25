// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
#include "stk_mesh/base/Types.hpp"      // for EntityRank, OrdinalVector, etc
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg, etc


namespace stk {
namespace mesh {


Relation::RawRelationType & Relation::RawRelationType::operator =(const Relation::RawRelationType & rhs)
{
    value = rhs.value;
    return *this;
}

//----------------------------------------------------------------------

namespace {

void get_entities_through_relations(
  const BulkData &mesh,
  Entity const *rels_begin,
  Entity const *rels_end,
  const Entity* i_beg ,
  const Entity* i_end ,
  std::vector<Entity> & entities_related )
{
  for (Entity const *rels_left = rels_begin ; rels_left != rels_end ; ++rels_left )
  {
    // Do all input entities have a relation to this entity ?

    Entity const e = *rels_left;
    EntityRank erank = mesh.entity_rank(e);

    const Entity* i = i_beg ;
    for ( ; i != i_end ; ++i )
    {
      const MeshIndex& meshIndex = mesh.mesh_index(*i);
      const Bucket* bucket = meshIndex.bucket;
      unsigned bucketOrd = meshIndex.bucket_ordinal;
      int num_conn = bucket->num_connectivity(bucketOrd, erank);
      const Entity* irels_j  = bucket->begin(bucketOrd, erank);
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
    const Entity* i = entities.data();
    const Entity* j = i+entities.size();

    const Bucket &ibucket = mesh.bucket(*i);
    const Ordinal &ibordinal = mesh.bucket_ordinal(*i);
    const EntityRank end_rank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());

    const Entity* next_i = i + 1;
    for (EntityRank rank = stk::topology::BEGIN_RANK; rank < end_rank; ++rank)
    {
      int num_conn   = mesh.num_connectivity(ibucket[ibordinal], rank);
      const Entity* rels_begin = mesh.begin(ibucket[ibordinal], rank);

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
  impl::find_entities_these_nodes_have_in_common(mesh, entities_related_rank,
                                                 entities.size(), entities.data(),
                                                 entities_related);
}

//----------------------------------------------------------------------

void induced_part_membership(const BulkData& mesh,
                             const Entity entity ,
                                   OrdinalVector & induced_parts)
{
  ThrowAssertMsg(mesh.is_valid(entity), "BulkData at " << &mesh << " does not know Entity" << entity.local_offset());

  const EntityRank e_rank = mesh.entity_rank(entity);
  const EntityRank end_rank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());

  for (EntityRank irank = static_cast<EntityRank>(e_rank + 1); irank < end_rank; ++irank)
  {
    int num_rels = mesh.num_connectivity(entity, irank);
    Entity const* rels     = mesh.begin(entity, irank);

    for (int j = 0; j < num_rels; ++j)
    {
      impl::get_part_ordinals_to_induce_on_lower_ranks(mesh, rels[j], e_rank, induced_parts);
    }
  }
}

//----------------------------------------------------------------------

stk::mesh::Entity find_by_ordinal(stk::mesh::Entity mesh_obj, stk::mesh::EntityRank rank, size_t ordinal, const stk::mesh::BulkData& mesh)
{
  stk::mesh::ConnectivityOrdinal const* ordinals = mesh.begin_ordinals(mesh_obj, rank);
  stk::mesh::Entity const* conn = mesh.begin(mesh_obj, rank);
  for (int i = 0, ie = mesh.num_connectivity(mesh_obj, rank); i < ie; ++i) {
    if (ordinals[i] == ordinal) {
      return conn[i];
    }
  }
  return stk::mesh::Entity();
}

} // namespace mesh
} // namespace stk
