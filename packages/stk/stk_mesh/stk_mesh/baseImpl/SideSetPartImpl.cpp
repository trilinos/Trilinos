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

#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_mesh/baseImpl/SideSetPartImpl.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {
namespace impl {

bool part_is_connected_to_shell_block(const BulkData& bulk, const stk::mesh::Part &part)
{
  bool connected = false;
  const MetaData& meta = bulk.mesh_meta_data();
  std::vector<const stk::mesh::Part*> touchingBlocks = meta.get_blocks_touching_surface(&part);

  for(const stk::mesh::Part* touchingBlock : touchingBlocks) {
    connected |= meta.get_topology(*touchingBlock).is_shell();
  }
  return connected;
}

void check_sideset_part_constraints(const BulkData& bulk, const stk::mesh::Part &part)
{
  const MetaData& meta = bulk.mesh_meta_data();
  if(part.primary_entity_rank() != meta.side_rank() &&
     !part_is_connected_to_shell_block(bulk, part)) {
    stk::RuntimeWarning() << "create_sideset: part " << part.name()
                                                     << " has rank " << part.primary_entity_rank();
  }
  if((part.id() == stk::mesh::Part::INVALID_ID) &&
     (part.name() != "universal_sideset") && !part.subsets().empty()) {
    stk::RuntimeWarning() << "create_sideset: part " << part.name()
                                                       << " has invalid id ";
  }

  for(const stk::mesh::Part* subsetPart : part.subsets()) {
    if(subsetPart->primary_entity_rank() == meta.side_rank()) {
      if(subsetPart->id() != part.id()) {
        stk::RuntimeWarning() << "create_sideset: part " << part.name()
                                                         << " with id " << part.id()
                                                         << "; subset sideblock part " << subsetPart->name()
                                                         << " has different id " << subsetPart->id();
      }
    }
  }
}

} // namespace impl
} // namespace mesh
} // namespace stk

