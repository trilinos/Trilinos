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

#include "stk_util/stk_config.h"
#include "stk_tools/mesh_tools/DetectHinges.hpp"
#include "stk_tools/mesh_tools/DetectHingesImpl.hpp"
#include "stk_tools/mesh_tools/DisconnectBlocksImpl.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Bucket.hpp"
#include <vector>

namespace stk {
namespace tools {

stk::mesh::EntityVector get_common_elements(const stk::mesh::BulkData& bulk, stk::mesh::Entity node1, stk::mesh::Entity node2)
{
  stk::mesh::EntityVector commonElements;
  stk::mesh::Entity nodes[] = {node1, node2};
  stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::ELEM_RANK, 2, nodes, commonElements);

  return commonElements;
}

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes)
{
  std::vector<std::string> blocksToDetect;
  fill_mesh_hinges( bulk,  blocksToDetect, hingeNodes);
}

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocksToDetect, HingeNodeVector& hingeNodes)
{
  hingeNodes = impl::get_hinge_nodes(bulk, blocksToDetect);
}

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, HingeNodeVector& hingeNodes, HingeEdgeVector& hingeEdges, bool onlyIfConnectedToSolidElements)
{
  std::vector<std::string> blocksToDetect;
  fill_mesh_hinges( bulk,  blocksToDetect, hingeNodes, hingeEdges, onlyIfConnectedToSolidElements);
}

void fill_mesh_hinges(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blocksToDetect, HingeNodeVector& hingeNodes, HingeEdgeVector& hingeEdges, bool onlyIfConnectedToSolidElements)
{
  hingeNodes = impl::get_hinge_nodes(bulk, blocksToDetect, onlyIfConnectedToSolidElements);

  if(hingeNodes.size() != 0) {
    hingeEdges = impl::get_hinge_edges(bulk, hingeNodes);
  }

  impl::prune_hinge_nodes(bulk, hingeNodes, hingeEdges);
}

}} // namespace stk::tools

