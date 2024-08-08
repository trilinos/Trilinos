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
#include "stk_mesh/base/Types.hpp"      // for EntityRank, OrdinalVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg, etc


namespace stk {
namespace mesh {


Relation::RawRelationType & Relation::RawRelationType::operator =(const Relation::RawRelationType & rhs)
{
    value = rhs.value;
    return *this;
}

void induced_part_membership(const BulkData& mesh,
                             const Entity entity ,
                                   OrdinalVector & induced_parts)
{
  STK_ThrowAssertMsg(mesh.is_valid(entity), "BulkData at " << &mesh << " does not know Entity" << entity.local_offset());

  const EntityRank e_rank = mesh.entity_rank(entity);
  const EntityRank end_rank = static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());

  for (EntityRank irank = static_cast<EntityRank>(e_rank + 1); irank < end_rank; ++irank)
  {
    const int num_rels = mesh.num_connectivity(entity, irank);
    if (num_rels > 0) {
      const Entity * rels     = mesh.begin(entity, irank);

      const Bucket* prevBucketPtr = nullptr;
      for (int j = 0; j < num_rels; ++j)
      {
        const Bucket* curBucketPtr = mesh.bucket_ptr(rels[j]);
        if (prevBucketPtr != curBucketPtr) {
          prevBucketPtr = curBucketPtr;
          impl::get_part_ordinals_to_induce_on_lower_ranks(mesh, *curBucketPtr, e_rank, induced_parts);
        }
      }
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
