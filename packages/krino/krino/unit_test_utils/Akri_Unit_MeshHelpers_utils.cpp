// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <set>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element

#include <Akri_Unit_MeshHelpers.hpp>

namespace krino {

void build_mesh(stk::mesh::BulkData & mesh, const std::vector<stk::mesh::EntityIdVector> & elem_nodes, const std::vector<int> & elem_procs, std::vector<stk::mesh::PartVector> & elem_parts)
{
  if (elem_nodes.empty()) return;

  const size_t numElem = elem_nodes.size();
  EXPECT_TRUE(elem_procs.size() == numElem);
  EXPECT_TRUE(elem_parts.size() == numElem);
  const size_t nodesPerElem = elem_nodes[0].size();
  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();

  mesh.modification_begin();

  // Create nodes and elements
  std::set<stk::mesh::EntityId> nodes_I_declared;
  for (size_t e=0; e<numElem; ++e)
  {
    if (parallel_size == 1 || elem_procs[e] == parallel_rank)
    {
      for(size_t n=0; n<nodesPerElem; ++n)
      {
        stk::mesh::EntityId nodeGlobalId = elem_nodes[e][n];

        stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeGlobalId);
        if (!mesh.is_valid(node))
        {
          mesh.declare_node(nodeGlobalId);
          nodes_I_declared.insert(nodeGlobalId);
        }
      }
      stk::mesh::EntityId elemGlobalId = e+1;
      stk::mesh::declare_element(mesh, elem_parts[e], elemGlobalId, elem_nodes[e]);
    }
  }

  // declare node sharing
  if (parallel_size > 1)
  {
    for (size_t e=0; e<numElem; ++e)
    {
      if (elem_procs[e] != parallel_rank)
      {
        for(size_t n=0; n<nodesPerElem; ++n)
        {
          stk::mesh::EntityId nodeGlobalId = elem_nodes[e][n];
          if (nodes_I_declared.find(nodeGlobalId) != nodes_I_declared.end())
          {
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeGlobalId);
            mesh.add_node_sharing(node, elem_procs[e]);
          }
        }
      }
    }
  }

  mesh.modification_end();
}

}
