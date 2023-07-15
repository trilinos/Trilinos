/*
 * Akri_Edge.cpp
 *
 *  Created on: Sep 23, 2022
 *      Author: drnoble
 */
#include "Akri_Edge.hpp"

#include <type_traits>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "Akri_MeshHelpers.hpp"

namespace krino {

static Edge edge_from_edge_node_offsets(stk::mesh::Entity::entity_value_type edgeNodeOffset0, stk::mesh::Entity::entity_value_type edgeNodeOffset1)
{
  static_assert(std::is_same<stk::mesh::Entity::entity_value_type, uint32_t>::value, "stk::mesh::Entity must be 32 bit.");
  const uint64_t edgeValue = (static_cast<uint64_t>(edgeNodeOffset1) << 32) + edgeNodeOffset0;
  return Edge(edgeValue);
}

Edge edge_from_edge_nodes(const stk::mesh::BulkData & mesh, stk::mesh::Entity edgeNode0, stk::mesh::Entity edgeNode1)
{
  return (mesh.identifier(edgeNode0) < mesh.identifier(edgeNode1)) ?
      edge_from_edge_node_offsets(edgeNode0.local_offset(), edgeNode1.local_offset()) :
      edge_from_edge_node_offsets(edgeNode1.local_offset(), edgeNode0.local_offset());
}

std::array<stk::mesh::Entity,2> get_edge_nodes(const Edge edge)
{
  static_assert(std::is_same<stk::mesh::Entity::entity_value_type, uint32_t>::value, "stk::mesh::Entity must be 32 bit.");
  return std::array<stk::mesh::Entity, 2>{stk::mesh::Entity(edge.value() & 0xFFFFFFFF), stk::mesh::Entity(edge.value() >> 32)};
}

void fill_edge_nodes(const Edge edge, std::vector<stk::mesh::Entity> & edgeNodes)
{
  static_assert(std::is_same<stk::mesh::Entity::entity_value_type, uint32_t>::value, "stk::mesh::Entity must be 32 bit.");
  edgeNodes.clear();
  edgeNodes.reserve(2);
  edgeNodes.emplace_back(edge.value() & 0xFFFFFFFF);
  edgeNodes.emplace_back(edge.value() >> 32);
}

void fill_entity_edges(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, std::vector<Edge> & entityEdges)
{
  const stk::topology topology = mesh.bucket(entity).topology();
  const unsigned numEdges = topology.num_edges();

  entityEdges.clear();
  entityEdges.reserve(numEdges);
  const stk::mesh::Entity * entityNodes = mesh.begin_nodes(entity);

  for (unsigned iEdge = 0; iEdge < numEdges; ++iEdge)
  {
    const unsigned * edgeNodeOrdinals = get_edge_node_ordinals(topology, iEdge);
    entityEdges.push_back(edge_from_edge_nodes(mesh, entityNodes[edgeNodeOrdinals[0]], entityNodes[edgeNodeOrdinals[1]]));
  }
}

int get_edge_parallel_owner_rank(const stk::mesh::BulkData & mesh, const Edge edge)
{
  const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
  return std::min(mesh.parallel_owner_rank(edgeNodes[0]), mesh.parallel_owner_rank(edgeNodes[1]));
}

std::string debug_edge(const stk::mesh::BulkData & mesh, const Edge edge)
{
  const std::array<stk::mesh::Entity,2> & edgeNodes = get_edge_nodes(edge);
  std::ostringstream out;
  out << "Edge with nodes " << mesh.identifier(edgeNodes[0]) << " and " << mesh.identifier(edgeNodes[1]);
  return out.str();
}

}
