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

#ifndef STK_MESH_CONNECTED_TOPOLOGY_NODES_HPP
#define STK_MESH_CONNECTED_TOPOLOGY_NODES_HPP

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>      // for EntityRank, ConnectivityOrdinal, EntityVector
#include <stk_topology/topology.hpp>
#include <vector>

namespace stk {
namespace mesh {

class ConnectedTopologyNodes {
public:
  ConnectedTopologyNodes(stk::topology fromTopology, const EntityVector& fromEntities)
  : m_connectedNodes(fromTopology.num_nodes()*fromEntities.size(), Entity()),
    m_ordinals(fromTopology.num_nodes()),
    m_numNodesPerEntity(fromTopology.num_nodes())
  {
    for (unsigned n=0; n<m_numNodesPerEntity; ++n) {
      m_ordinals[n] = n;
    }
  }

  ~ConnectedTopologyNodes(){}

  unsigned num_nodes_per_entity(unsigned fromEntity) const { return m_numNodesPerEntity; }

  const Entity* begin_connectivity(unsigned fromEntity) const
  { return m_connectedNodes.data() + fromEntity*m_numNodesPerEntity; }

  const Entity* end_connectivity(unsigned fromEntity) const
  { return m_connectedNodes.data() + ((fromEntity+1)*m_numNodesPerEntity); }

  const ConnectivityOrdinal* begin_ordinals(unsigned fromEntity) const
  { return m_ordinals.data(); }

  const ConnectivityOrdinal* end_ordinals(unsigned fromEntity) const
  { return m_ordinals.data() + m_numNodesPerEntity; }

  bool set_connected_node(unsigned fromEntity, Entity connectedNode, ConnectivityOrdinal ordinal)
  {
     if (m_connectedNodes[fromEntity*m_numNodesPerEntity + ordinal] != connectedNode) {
       m_connectedNodes[fromEntity*m_numNodesPerEntity + ordinal] = connectedNode;
       return true;
     }
     return false;
  }

  void set_connected_nodes(unsigned fromEntity, unsigned numNodes, const Entity* connectedNodes)
  {
    const unsigned start = fromEntity*m_numNodesPerEntity;
    for(unsigned i=0; i<numNodes; ++i) {
      m_connectedNodes[start+i] = connectedNodes[i];
    }
  }

private:
  EntityVector m_connectedNodes;
  std::vector<ConnectivityOrdinal> m_ordinals;
  unsigned m_numNodesPerEntity;
};

} // namespace mesh
} // namespace stk

#endif //STK_MESH_CONNECTED_TOPOLOGY_NODES_HPP

