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

#include <stk_mesh/base/DestroyRelations.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk { namespace mesh {

void destroy_relations(stk::mesh::BulkData &bulk,
                       stk::mesh::Entity entity,
                       stk::mesh::EntityRank connectedRank)
{
  const int numConn = bulk.num_connectivity(entity, connectedRank);
  const Entity* conn = bulk.begin(entity, connectedRank);
  const ConnectivityOrdinal* ords = bulk.begin_ordinals(entity, connectedRank);

  switch(numConn) {
  case 0:
    return; break;
  case 1:
    if (bulk.entity_rank(entity) > connectedRank) {
      bulk.destroy_relation(entity, conn[0], ords[0]);
    }
    else {
      bulk.destroy_relation(conn[0], entity, ords[0]);
    }
    break;
  case 2:
    {
      Entity tmpEntity = conn[1];
      ConnectivityOrdinal tmpOrd = ords[1];
      if (bulk.entity_rank(entity) > connectedRank) {
        bulk.destroy_relation(entity, conn[0], ords[0]);
        bulk.destroy_relation(entity, tmpEntity, tmpOrd);
      }
      else {
        bulk.destroy_relation(conn[0],   entity, ords[0]);
        bulk.destroy_relation(tmpEntity, entity, tmpOrd);
      }
      break;
    }
  default:
    {
      stk::mesh::EntityVector connv(conn, bulk.end(entity, connectedRank));
      std::vector<ConnectivityOrdinal> ordv(ords, bulk.end_ordinals(entity, connectedRank));
      if (bulk.entity_rank(entity) > connectedRank) {
        for(int i=0; i<numConn; ++i) {
          bulk.destroy_relation(entity, connv[i], ordv[i]);
        }
      }
      else {
        for(int i=0; i<numConn; ++i) {
          bulk.destroy_relation(connv[i], entity, ordv[i]);
        }
      }
      break;
    }
  };
}

}} // namespace stk::mesh

