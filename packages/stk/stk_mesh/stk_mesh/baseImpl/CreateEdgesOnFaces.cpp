#include "stk_mesh/baseImpl/CreateEdgesOnFaces.hpp"

#include <vector>
#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/EntityLess.hpp"


namespace stk {
namespace mesh {
namespace impl {

Part* EdgesOnFacesCreator::create_edges()
{
  m_edgePart = create_edge_part_if_needed();
  if (m_edgePart)
  {
    check_edge_part_topology();
    create_surface_edges();
    attach_edges_to_non_owned_faces();
  }

  return m_edgePart;
}

void EdgesOnFacesCreator::create_surface_edges()
{
  std::vector<stk::mesh::EntityId> newEntityIds;
  m_bulk.generate_new_ids(stk::topology::EDGE_RANK, compute_upper_bound_on_num_edges(), newEntityIds);

  std::vector<stk::mesh::Part*> addParts = {m_edgePart};
  std::vector<stk::mesh::Entity> edgeNodes;
  size_t newEntityIdx = 0;
  for (stk::mesh::Bucket* bucket : m_bulk.get_buckets(stk::topology::FACE_RANK, m_surfaceAndOwned))
  {
    stk::topology faceTopo = bucket->topology();
    check_face_topology(faceTopo);

    for (stk::mesh::Entity face : *bucket)
    {
      stk::mesh::Entity const* faceNodes = m_bulk.begin_nodes(face);
      for (unsigned int i=0; i < faceTopo.num_edges(); ++i)
      {
        edgeNodes.clear();
        faceTopo.edge_nodes(faceNodes, i, std::back_inserter(edgeNodes));

        stk::mesh::Entity edge = get_common_edge(edgeNodes[0], edgeNodes[1]);
        if (!m_bulk.is_valid(edge))
        {
          stk::mesh::EntityId entityId = newEntityIds.at(newEntityIdx++);
          edge = m_bulk.declare_edge(entityId, addParts);

          sort_nodes_for_global_consistency(edgeNodes);
          for (size_t j=0; j < edgeNodes.size(); ++j)
          {
            m_bulk.declare_relation(edge, edgeNodes[j], j);
          }
        }

        m_bulk.declare_relation(face, edge, i);
      }
    }
  }
}

void EdgesOnFacesCreator::attach_edges_to_non_owned_faces()
{
  stk::mesh::Selector surfaceAndNotOwned = m_surfaceSelector &
    (m_bulk.mesh_meta_data().universal_part() - m_bulk.mesh_meta_data().locally_owned_part());

  std::vector<stk::mesh::Entity> faceNodes, edgeNodes;
  for (stk::mesh::Bucket* bucket : m_bulk.get_buckets(stk::topology::FACE_RANK, surfaceAndNotOwned))
  {
    stk::topology faceTopo = bucket->topology();
    check_face_topology(faceTopo);
    for (stk::mesh::Entity face : *bucket)
    {
      faceNodes.assign(m_bulk.begin_nodes(face), m_bulk.end_nodes(face));
      for (unsigned int i=0; i < faceTopo.num_edges(); ++i)
      {
        edgeNodes.clear();
        faceTopo.edge_nodes(faceNodes, i, std::back_inserter(edgeNodes));
        stk::mesh::Entity edge = get_common_edge(edgeNodes[0], edgeNodes[1]);
        if (m_bulk.is_valid(edge))
        {
          m_bulk.declare_relation(face, edge, i);
        }
      }
    }
  }
}

stk::topology EdgesOnFacesCreator::get_edge_topology() const
{
  stk::topology edgeTopo = stk::topology::INVALID_TOPOLOGY;
  for (stk::mesh::Bucket* bucket : m_bulk.get_buckets(stk::topology::FACE_RANK, m_surfaceAndOwned))
  {
    edgeTopo = bucket->topology().edge_topology();
    break;
  }

  int edgeTopoLocal = static_cast<int>(edgeTopo);
  int edgeTopoGlobal = 0;
  MPI_Allreduce(&edgeTopoLocal, &edgeTopoGlobal, 1, MPI_INT, MPI_MAX, m_bulk.parallel());

  return static_cast<stk::topology::topology_t>(edgeTopoGlobal);
}

stk::mesh::Part* EdgesOnFacesCreator::create_edge_part_if_needed() const
{
  stk::mesh::Part* edgePart = m_edgePart;
  if (!m_edgePart)
  {
    stk::topology edgeTopo = get_edge_topology();
    if (edgeTopo == stk::topology::INVALID_TOPOLOGY)
    {
      return nullptr;
    }

    edgePart = &(m_bulk.mesh_meta_data().declare_part_with_topology("face_edges_part", edgeTopo));
  }

  return edgePart;
}

void EdgesOnFacesCreator::check_edge_part_topology()
{
  stk::topology topo = m_edgePart->topology();
  STK_ThrowRequireMsg(topo.rank() == stk::topology::EDGE_RANK, "edge part topology must have EDGE_RANK");
  STK_ThrowRequireMsg(topo.num_vertices() <= 3, "edge part topology must have 3 or fewer vertices");
}

int EdgesOnFacesCreator::compute_upper_bound_on_num_edges() const
{
  std::vector<size_t> localEntityCounts;
  stk::mesh::count_entities(m_surfaceAndOwned, m_bulk, localEntityCounts);
  const int maxEdgesPerFace = 4;

  return maxEdgesPerFace * localEntityCounts[stk::topology::FACE_RANK] - localEntityCounts[stk::topology::EDGE_RANK];
}

void EdgesOnFacesCreator::check_face_topology(stk::topology faceTopo)
{
  for (unsigned int i=0; i < faceTopo.num_edges(); ++i)
  {
    STK_ThrowRequireMsg(faceTopo.edge_topology(i) == m_edgePart->topology(), "edge part topology does not match face edge topology");
  }
}

[[nodiscard]] stk::mesh::Entity EdgesOnFacesCreator::get_common_edge(stk::mesh::Entity entity1, stk::mesh::Entity entity2) const
{
  stk::mesh::Entity const* endEdge1 = m_bulk.end_edges(entity1);
  stk::mesh::Entity const* endEdge2 = m_bulk.end_edges(entity2);
  for (stk::mesh::Entity const * edge1 = m_bulk.begin_edges(entity1); edge1 != endEdge1; ++edge1)
  {
    for (stk::mesh::Entity const * edge2 = m_bulk.begin_edges(entity2); edge2 != endEdge2; ++edge2)
    {
      if (*edge1 == *edge2)
      {
        return *edge1;
      }
    }
  }

  return stk::mesh::Entity();
}

void EdgesOnFacesCreator::sort_nodes_for_global_consistency(std::vector<Entity>& edgeNodes) const
{
  stk::mesh::EntityLess entityLess(m_bulk);
  if (!entityLess(edgeNodes[0], edgeNodes[1]))
  {
    stk::mesh::Entity tmp = edgeNodes[0];
    edgeNodes[0] = edgeNodes[1];
    edgeNodes[1] = tmp;
  }
}

}
}
}