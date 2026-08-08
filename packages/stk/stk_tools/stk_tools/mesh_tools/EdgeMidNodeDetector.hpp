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

#ifndef STK_STK_TOOLS_STK_TOOLS_MESH_TOOLS_EDGEMIDNODEDETECTOR_HPP_
#define STK_STK_TOOLS_STK_TOOLS_MESH_TOOLS_EDGEMIDNODEDETECTOR_HPP_

#include "stk_mesh/base/Types.hpp"

#include <string>
#include <algorithm>
#include <map>
#include <utility>

namespace stk { namespace mesh { class BulkData; } }

namespace stk {
namespace tools {

using EdgeVertices = std::pair<stk::mesh::Entity, stk::mesh::Entity>;

class EdgeMidNodeDetector
{
public:

  EdgeMidNodeDetector(const stk::mesh::BulkData& bulk);

  const stk::mesh::BulkData& get_bulk_data() const { return m_bulk; }

  stk::mesh::BulkData& get_bulk_data() { return const_cast<stk::mesh::BulkData&> (m_bulk); }

  void fill_all_edge_mid_nodes(stk::mesh::EntityVector& edgeMidNodes) const;

  stk::mesh::EntityVector get_all_edge_mid_nodes() const;

  void fill_all_edges(std::vector<EdgeVertices>& edges) const;

  std::vector<EdgeVertices> get_all_edges() const;

  const stk::mesh::EntityVector& get_mid_nodes(const EdgeVertices& edge) const;

  void fill_mid_nodes(const EdgeVertices& edge, stk::mesh::EntityVector& nodes) const;

  bool is_mid_edge_node(stk::mesh::Entity node) const;

  bool is_mid_edge_node(stk::mesh::EntityId nodeId) const;

  const std::vector<EdgeVertices>& get_edge_vertices(stk::mesh::Entity midEdgeNode) const;

  void fill_edge_vertices(stk::mesh::Entity midEdgeNode, std::vector<EdgeVertices>& edges) const;

  size_t num_mid_nodes() const { return m_midEdgeNodeToEdgeMap.size(); }

  size_t num_edges() const { return m_edgeToMidEdgeNodesMap.size(); }

private:
  const stk::mesh::BulkData& m_bulk;
  stk::mesh::EntityVector m_scratchElementNodes;
  stk::mesh::EntityVector m_scratchEdgeMidNodes;

  const std::vector<EdgeVertices> m_emptyEdges;
  const stk::mesh::EntityVector m_emptyNodes;

  std::map<stk::mesh::Entity, std::vector<EdgeVertices>> m_midEdgeNodeToEdgeMap;
  std::map<EdgeVertices, stk::mesh::EntityVector> m_edgeToMidEdgeNodesMap;

private:
  void initialize();

  bool fill_edge(const stk::mesh::EntityVector& allEdgeNodes, EdgeVertices& edge);

  void fill_mid_nodes(const stk::mesh::EntityVector& allEdgeNodes, stk::mesh::EntityVector& edgeMidNodes);

  void update_maps();

  void process_edge_mid_nodes(stk::mesh::Entity element, stk::topology elemTopology);

  void process_non_vertex_mid_nodes(stk::mesh::Entity element, stk::topology elemTopology);
};

}
}

#endif /* STK_STK_TOOLS_STK_TOOLS_MESH_TOOLS_EDGEMIDNODEDETECTOR_HPP_ */
