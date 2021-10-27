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

#ifndef STK_MESH_CONNECTED_SPARSE_NODES_HPP
#define STK_MESH_CONNECTED_SPARSE_NODES_HPP

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/SparseConnectivity.hpp>
#include <stk_topology/topology.hpp>
#include <vector>

namespace stk {
namespace mesh {

class ConnectedSparseNodes {
public:
  ConnectedSparseNodes(SparseConnectivity& sparseConnectivity,
                       EntityRank fromRank,
                       const EntityVector& fromEntities)
  : m_sparseConnectivity(sparseConnectivity),
    m_fromEntities(fromEntities),
    m_fromRank(fromRank)
  {}

  ~ConnectedSparseNodes(){}

  unsigned num_nodes_per_entity(unsigned fromEntity) const
  { return m_sparseConnectivity.num_connectivity(m_fromEntities[fromEntity], stk::topology::NODE_RANK); }

  const Entity* begin_connectivity(unsigned fromEntity) const
  { return m_sparseConnectivity.begin_connectivity(m_fromEntities[fromEntity], stk::topology::NODE_RANK); }

  const Entity* end_connectivity(unsigned fromEntity) const
  { return m_sparseConnectivity.end_connectivity(m_fromEntities[fromEntity], stk::topology::NODE_RANK); }

  const ConnectivityOrdinal* begin_ordinals(unsigned fromEntity) const
  { return m_sparseConnectivity.begin_ordinals(m_fromEntities[fromEntity], stk::topology::NODE_RANK); }

  const ConnectivityOrdinal* end_ordinals(unsigned fromEntity) const
  { return m_sparseConnectivity.end_ordinals(m_fromEntities[fromEntity], stk::topology::NODE_RANK); }

  bool set_connected_node(unsigned fromEntity, Entity connectedNode, ConnectivityOrdinal ordinal)
  {
    if (connectedNode.local_offset() > 0) {
      return m_sparseConnectivity.add_connectivity(m_fromRank, m_fromEntities[fromEntity],
                                                   stk::topology::NODE_RANK, connectedNode, ordinal,
                                                   INVALID_PERMUTATION);
    }
    else {
      const unsigned numBefore = m_sparseConnectivity.num_connectivity(m_fromEntities[fromEntity], stk::topology::NODE_RANK);
      bool returnVal = m_sparseConnectivity.remove_connectivity(m_fromRank, m_fromEntities[fromEntity],
                                                                stk::topology::NODE_RANK, connectedNode, ordinal);
      const unsigned numAfter = m_sparseConnectivity.num_connectivity(m_fromEntities[fromEntity], stk::topology::NODE_RANK);
      ThrowRequireMsg((numBefore==0 && numAfter==0) || (numBefore>0 && numAfter==numBefore-1),"remove didn't work");
      return returnVal;
    }
  }

  bool reset_connected_node(unsigned fromEntity, Entity connectedNode, ConnectivityOrdinal ordinal)
  {
     return m_sparseConnectivity.replace_or_add_connectivity(m_fromRank, m_fromEntities[fromEntity],
                                                             stk::topology::NODE_RANK, connectedNode, ordinal);
  }

  void reset_connected_nodes(unsigned fromEntity, unsigned numNodes, const Entity* connectedNodes)
  {
    for(unsigned i=0; i<numNodes; ++i) {
      m_sparseConnectivity.replace_or_add_connectivity(m_fromRank, m_fromEntities[fromEntity+i],
                                                       stk::topology::NODE_RANK, connectedNodes[i], i);
    }
  }

private:
  SparseConnectivity& m_sparseConnectivity;
  const EntityVector& m_fromEntities;
  EntityRank m_fromRank;
};

} // namespace mesh
} // namespace stk

#endif //STK_MESH_CONNECTED_SPARSE_NODES_HPP

