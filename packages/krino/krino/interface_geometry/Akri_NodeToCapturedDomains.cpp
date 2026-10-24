// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshHelpers.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/CommSparse.hpp>

namespace krino {

static void pack_node_captured_domains_for_sharing_or_ghosting_procs(const stk::mesh::BulkData & mesh,
    const stk::mesh::Entity node,
    const std::vector<int> & nodeCapturedDomains,
    std::vector<int> procsThatNeedToKnowAboutNode,
    stk::CommSparse &commSparse)
{
  fill_procs_owning_or_sharing_or_ghosting_node(mesh, node, procsThatNeedToKnowAboutNode);
  for (int procId : procsThatNeedToKnowAboutNode)
  {
    if (procId != mesh.parallel_rank())
    {
      commSparse.send_buffer(procId).pack(mesh.identifier(node));
      commSparse.send_buffer(procId).pack(nodeCapturedDomains.size());
      for (int domain : nodeCapturedDomains)
        commSparse.send_buffer(procId).pack(domain);
    }
  }
}

static
void pack_all_node_captured_domains_for_sharing_or_ghosting_procs(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    stk::CommSparse &commSparse)
{
  std::vector<int> scratchVec;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto && entry : nodesToCapturedDomains)
    {
      stk::mesh::Entity node = entry.first;
      if (mesh.parallel_owner_rank(node) == mesh.parallel_rank())
      {
        const auto & nodeCapturedDomains = entry.second;
        pack_node_captured_domains_for_sharing_or_ghosting_procs(mesh, node, nodeCapturedDomains, scratchVec, commSparse);
      }
    }
  });
}

static
void pack_node_captured_domains_of_given_nodes_for_sharing_or_ghosting_procs(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & nodes,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    stk::CommSparse &commSparse)
{
  std::vector<int> scratchVec;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (auto node : nodes)
    {
      if (mesh.parallel_owner_rank(node) == mesh.parallel_rank())
      {
        const auto nodeCapturedDomains = nodesToCapturedDomains.at(node);
        pack_node_captured_domains_for_sharing_or_ghosting_procs(mesh, node, nodeCapturedDomains, scratchVec, commSparse);
      }
    }
  });
}

static
void unpack_snap_node_domains(const stk::mesh::BulkData & mesh,
    NodeToCapturedDomainsMap & nodesToSnappedDomains,
    stk::CommSparse &commSparse)
{
  std::vector<int> nodeSnappedDomains;

  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId nodeId;
      commSparse.recv_buffer(procId).unpack(nodeId);
      size_t numNodes;
      commSparse.recv_buffer(procId).unpack(numNodes);
      nodeSnappedDomains.resize(numNodes);
      for (auto && domain : nodeSnappedDomains)
        commSparse.recv_buffer(procId).unpack(domain);

      stk::mesh::Entity snapNode = mesh.get_entity(stk::topology::NODE_RANK, nodeId);
      nodesToSnappedDomains[snapNode] = nodeSnappedDomains;
    }
  });
}

void communicate_node_captured_domains_for_given_nodes(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & nodes,
    NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_node_captured_domains_of_given_nodes_for_sharing_or_ghosting_procs(mesh, nodes, nodesToCapturedDomains, commSparse);
  unpack_snap_node_domains(mesh, nodesToCapturedDomains, commSparse);
}

void communicate_node_captured_domains_for_all_nodes(const stk::mesh::BulkData & mesh,
    NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_all_node_captured_domains_for_sharing_or_ghosting_procs(mesh, nodesToCapturedDomains, commSparse);
  unpack_snap_node_domains(mesh, nodesToCapturedDomains, commSparse);
}

}


