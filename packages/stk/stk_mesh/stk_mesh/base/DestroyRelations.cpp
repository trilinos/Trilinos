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

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/DestroyRelations.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace stk { namespace mesh {

void destroy_relations(stk::mesh::BulkData &bulk,
                       stk::mesh::Entity entity,
                       stk::mesh::EntityRank connectedRank)
{
  const bool downward = bulk.entity_rank(entity) > connectedRank;
  const int numConn = bulk.num_connectivity(entity, connectedRank);
  for(int j = numConn - 1; j>= 0; --j) {
    const MeshIndex& meshIdx = bulk.mesh_index(entity);
    const Entity* conn = meshIdx.bucket->begin(meshIdx.bucket_ordinal, connectedRank);
    const ConnectivityOrdinal* ords = meshIdx.bucket->begin_ordinals(meshIdx.bucket_ordinal, connectedRank);
    if (bulk.is_valid(conn[j])) {
      if (downward) {
        STK_ThrowRequireMsg(bulk.destroy_relation(entity, conn[j], ords[j]), "Failed to destroy relation "<<bulk.entity_key(entity)<<" -> "<<bulk.entity_key(conn[j]));
      }
      else {
        STK_ThrowRequireMsg(bulk.destroy_relation(conn[j], entity, ords[j]), "Failed to destroy relation "<<bulk.entity_key(conn[j])<<" -> "<<bulk.entity_key(entity));
      }
    }
  }
}

}} // namespace stk::mesh

