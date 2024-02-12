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

#ifndef stk_mesh_base_MeshUtils_hpp
#define stk_mesh_base_MeshUtils_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

namespace stk {
namespace mesh {

// Function to promote ghosts to shared:
// This function is a transition aid until an application can consistently and
// correctly call BulkData::add_node_sharing because this function is not
// expected to perform well or to scale well at large processor counts.
void fixup_ghosted_to_shared_nodes(BulkData& bulk);

// Helper functions:
void find_ghosted_nodes_that_need_to_be_shared(const BulkData & bulk, stk::mesh::EntityVector& ghosted_nodes_that_are_now_shared);

void cleanup_all_downward_orphan_entities(BulkData& bulk, stk::mesh::EntityRank maxOrphanRank = stk::topology::FACE_RANK);
void cleanup_orphan_nodes(BulkData& bulk);
void cleanup_orphan_edges(BulkData& bulk);
void cleanup_orphan_faces(BulkData& bulk);

} // namespace mesh
} // namespace stk

#endif // stk_mesh_base_MeshUtils_hpp
