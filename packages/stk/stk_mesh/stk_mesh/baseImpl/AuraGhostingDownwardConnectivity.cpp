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

#include <stk_mesh/baseImpl/AuraGhostingDownwardConnectivity.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/Visitors.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/EntityProcMapping.hpp>

namespace stk {
namespace mesh {
namespace impl {

AuraGhostingDownwardConnectivity::AuraGhostingDownwardConnectivity()
 : AuraGhosting()
{
}

AuraGhostingDownwardConnectivity::~AuraGhostingDownwardConnectivity()
{
}

void AuraGhostingDownwardConnectivity::generate_aura(BulkData& bulkData)
{
  EntityProcMapping sendAuraEntityProcs(bulkData.get_size_of_entity_index_space());
  const MetaData& meta = bulkData.mesh_meta_data();

  std::vector<int> sharingProcs;
  for(EntityRank rank : {stk::topology::ELEM_RANK, stk::topology::CONSTRAINT_RANK}) {
    impl::for_each_selected_entity_run_no_threads(bulkData, rank, meta.locally_owned_part(),
      [&sendAuraEntityProcs, &sharingProcs](const BulkData& bulk, const MeshIndex& meshIndex) {
        Entity entity = (*meshIndex.bucket)[meshIndex.bucket_ordinal];
        for(EntityRank lowerRank : {stk::topology::NODE_RANK, stk::topology::EDGE_RANK, stk::topology::FACE_RANK}) {
          const unsigned numConn = bulk.num_connectivity(entity, lowerRank);
          const Entity* conn = bulk.begin(entity, lowerRank);

          for(unsigned i=0; i<numConn; ++i) {
            if (bulk.bucket(conn[i]).shared()) {
              bulk.comm_shared_procs(conn[i], sharingProcs);

              for(int p : sharingProcs) {
                sendAuraEntityProcs.addEntityProc(entity, p);
              }
            }
          }
        }
      }
    );
  }

  change_ghosting(bulkData, sendAuraEntityProcs);
}

}}} // end namepsace stk mesh impl

