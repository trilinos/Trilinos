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

#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>
#include <stk_tools/mesh_tools/EdgeMidNodeDetector.hpp>

namespace stk {
namespace tools {

  EdgeMidNodeDetector::EdgeMidNodeDetector(const stk::mesh::BulkData& bulk)
  : m_bulk(bulk)
  {
    initialize();
  }

  void EdgeMidNodeDetector::fill_all_edge_mid_nodes(stk::mesh::EntityVector& edgeMidNodes) const
  {
    edgeMidNodes.clear();
    edgeMidNodes.reserve(m_midEdgeNodeToEdgeMap.size());

    for(const auto& entry : m_midEdgeNodeToEdgeMap) {
      edgeMidNodes.push_back(entry.first);
    }
  }

  stk::mesh::EntityVector EdgeMidNodeDetector::get_all_edge_mid_nodes() const
  {
    stk::mesh::EntityVector edgeMidNodes;
    fill_all_edge_mid_nodes(edgeMidNodes);
    return edgeMidNodes;
  }

  void EdgeMidNodeDetector::fill_all_edges(std::vector<EdgeVertices>& edges) const
  {
    edges.clear();
    edges.reserve(m_edgeToMidEdgeNodesMap.size());

    for(const auto& entry : m_edgeToMidEdgeNodesMap) {
      edges.push_back(entry.first);
    }
  }

  std::vector<EdgeVertices> EdgeMidNodeDetector::get_all_edges() const
  {
    std::vector<EdgeVertices> edges;
    fill_all_edges(edges);
    return edges;
  }

  const stk::mesh::EntityVector& EdgeMidNodeDetector::get_mid_nodes(const EdgeVertices& edge) const
  {
    auto iter = m_edgeToMidEdgeNodesMap.find(edge);
    if(iter == m_edgeToMidEdgeNodesMap.end()) return m_emptyNodes;

    return iter->second;
  }

  void EdgeMidNodeDetector::fill_mid_nodes(const EdgeVertices& edge, stk::mesh::EntityVector& nodes) const
  {
    nodes.clear();

    auto iter = m_edgeToMidEdgeNodesMap.find(edge);
    if(iter == m_edgeToMidEdgeNodesMap.end()) return;

    nodes = iter->second;
  }

  bool EdgeMidNodeDetector::is_mid_edge_node(stk::mesh::Entity node) const
  {
    return (m_midEdgeNodeToEdgeMap.find(node) != m_midEdgeNodeToEdgeMap.end());
  }

  bool EdgeMidNodeDetector::is_mid_edge_node(stk::mesh::EntityId nodeId) const
  {
    stk::mesh::Entity node = m_bulk.get_entity(stk::topology::NODE_RANK, nodeId);
    if(!m_bulk.is_valid(node)) return false;

    return is_mid_edge_node(node);
  }

  const std::vector<EdgeVertices>& EdgeMidNodeDetector::get_edge_vertices(stk::mesh::Entity midEdgeNode) const
  {
    auto iter = m_midEdgeNodeToEdgeMap.find(midEdgeNode);
    if(iter == m_midEdgeNodeToEdgeMap.end()) return m_emptyEdges;

    return iter->second;
  }

  void EdgeMidNodeDetector::fill_edge_vertices(stk::mesh::Entity midEdgeNode, std::vector<EdgeVertices>& edges) const
  {
    edges.clear();

    auto iter = m_midEdgeNodeToEdgeMap.find(midEdgeNode);
    if(iter == m_midEdgeNodeToEdgeMap.end()) return;

    edges = iter->second;
  }

  void EdgeMidNodeDetector::initialize()
  {
    const stk::mesh::MetaData& meta = m_bulk.mesh_meta_data();
    const stk::mesh::BucketVector &buckets = m_bulk.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());

    m_scratchElementNodes.clear();
    m_scratchEdgeMidNodes.clear();
    m_midEdgeNodeToEdgeMap.clear();
    m_edgeToMidEdgeNodesMap.clear();

    for(const stk::mesh::Bucket* bucket : buckets) {
      stk::topology elemTopology = bucket->topology();

      unsigned numSubTopology = elemTopology.num_sub_topology(stk::topology::EDGE_RANK);
      if(numSubTopology != 0) {
        for(stk::mesh::Entity element : *bucket) {
          process_edge_mid_nodes(element, elemTopology);
        }
      } else {
        for(stk::mesh::Entity element : *bucket) {
          process_non_vertex_mid_nodes(element, elemTopology);
        }
      }
    }
  }

  bool EdgeMidNodeDetector::fill_edge(const stk::mesh::EntityVector& allEdgeNodes, EdgeVertices& edge)
  {
    if(allEdgeNodes.size() < 2) return false;

    if(m_bulk.identifier(allEdgeNodes[0]) < m_bulk.identifier(allEdgeNodes[1])) {
      edge.first = allEdgeNodes[0];
      edge.second = allEdgeNodes[1];
    } else {
      edge.first = allEdgeNodes[1];
      edge.second = allEdgeNodes[0];
    }

    return true;
  }

  void EdgeMidNodeDetector::fill_mid_nodes(const stk::mesh::EntityVector& allEdgeNodes, stk::mesh::EntityVector& edgeMidNodes)
  {
    edgeMidNodes.clear();

    for(size_t i=2; i<allEdgeNodes.size(); i++) {
      edgeMidNodes.push_back(allEdgeNodes[i]);
    }
  }

  void EdgeMidNodeDetector::update_maps()
  {
    EdgeVertices edge;
    if(fill_edge(m_scratchElementNodes, edge)) {
      if(m_edgeToMidEdgeNodesMap.find(edge) != m_edgeToMidEdgeNodesMap.end()) return;

      fill_mid_nodes(m_scratchElementNodes, m_scratchEdgeMidNodes);

      m_edgeToMidEdgeNodesMap[edge] = m_scratchEdgeMidNodes;

      for(stk::mesh::Entity node : m_scratchEdgeMidNodes) {
        auto iter = m_midEdgeNodeToEdgeMap.find(node);
        if(iter == m_midEdgeNodeToEdgeMap.end()) {
          m_midEdgeNodeToEdgeMap[node] = std::vector<EdgeVertices>{edge};
        } else {
          iter->second.push_back(edge);
        }
      }
    }
  }

  void EdgeMidNodeDetector::process_edge_mid_nodes(stk::mesh::Entity element, stk::topology elemTopology)
  {
    const stk::mesh::Entity* elemNodes = m_bulk.begin_nodes(element);
    unsigned numSubTopology = elemTopology.num_sub_topology(stk::topology::EDGE_RANK);

    for(unsigned ordinal=0; ordinal<numSubTopology; ordinal++) {
      stk::topology subTopology = elemTopology.sub_topology(stk::topology::EDGE_RANK, ordinal);
      unsigned numSubTopologyNodes = subTopology.num_nodes();
      m_scratchElementNodes.resize(numSubTopologyNodes);
      elemTopology.sub_topology_nodes(elemNodes, stk::topology::EDGE_RANK, ordinal, m_scratchElementNodes.data());

      update_maps();
    }
  }

  void EdgeMidNodeDetector::process_non_vertex_mid_nodes(stk::mesh::Entity element, stk::topology elemTopology)
  {
    const stk::mesh::Entity* elemNodes = m_bulk.begin_nodes(element);
    unsigned numNodes = m_bulk.num_nodes(element);
    unsigned numVertices = elemTopology.num_vertices();

    STK_ThrowRequireMsg(numVertices <= 2, "Expected less than 2 vertices for topology: " << elemTopology);

    m_scratchElementNodes.clear();
    for(unsigned i=0; i<numNodes; i++) {
      m_scratchElementNodes.push_back(elemNodes[i]);
    }

    update_maps();
  }
}
}

