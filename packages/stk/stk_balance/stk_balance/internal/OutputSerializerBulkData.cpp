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

#include "OutputSerializerBulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/GetEntities.hpp"

namespace stk {
namespace balance {

OutputSerializerBulkData::OutputSerializerBulkData(ParallelMachine parallel)
  : BulkData(std::make_shared<stk::mesh::MetaData>(), parallel, stk::mesh::BulkData::NO_AUTO_AURA
#ifdef SIERRA_MIGRATION
        , true
#endif
             )
{
}

OutputSerializerBulkData::OutputSerializerBulkData(unsigned spatialDim, ParallelMachine parallel)
  : BulkData(std::make_shared<stk::mesh::MetaData>(), parallel, stk::mesh::BulkData::NO_AUTO_AURA
#ifdef SIERRA_MIGRATION
        , true
#endif
             )
{
}

void
OutputSerializerBulkData::switch_to_serial_mesh()
{
  modification_begin();

  stk::mesh::EntityVector nodesToUnshare;
  stk::mesh::get_entities(*this, stk::topology::NODE_RANK, mesh_meta_data().globally_shared_part(), nodesToUnshare);

  for (const stk::mesh::Entity& node : nodesToUnshare) {
    internal_set_owner(node, 0);
    remove_entity_comm(node);
    entity_comm_map_clear(entity_key(node));
  }

  destroy_all_ghosting();

  stk::mesh::PartVector addParts{&mesh_meta_data().locally_owned_part()};
  stk::mesh::PartVector removeParts{&mesh_meta_data().globally_shared_part()};
  internal_verify_and_change_entity_parts(nodesToUnshare, addParts, removeParts);

  m_parallel = Parallel(MPI_COMM_SELF);

  modification_end();
}

}
}
